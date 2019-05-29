/*!
  * variants v2.1.0
  * Copyright 2018 IUCT-O
  * Author Frederic Escudie
  * Licensed under GNU General Public License
  */


function getURLFromVariantId(variant_id){
    let url = null
    if(variant_id.startsWith("COSM")){
        url = `https://cancer.sanger.ac.uk/cosmic/search?q=${variant_id}`
    } else if(variant_id.startsWith("rs")){
        url = `https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=${variant_id}`
    }
    return url
}


function getVariantLink(variant_id){
    const url = getURLFromVariantId(variant_id)
    let variant_link = variant_id
    if(url !== null){
        variant_link = `<a target="_blank" href="${url}">${variant_id}</a>`
    }
    return variant_link
}


class VariantAnnotSubject {
    constructor(feature, feature_type, symbol){
        this.feature = feature // "ENST00000369535"
        this.feature_type = feature_type // "Transcript"
        this.symbol = symbol // "NRAS"
    }

    static fromJSON(hash){
        return new VariantAnnotSubject(hash.feature, hash.feature_type, hash.symbol)
    }
}

class VariantAnnot {
    constructor(subject, changes={}, conseq="", pathogenicity=[]){
        this.subject = subject
        this.changes = changes // { "HGVSc": "ENST00000369535.4:c.318G>T", "HGVSp": "ENSP00000358548.4:p.Ser106%3D" }  OR  { "HGVSc":"NM_022965.3:c.931-528delG", "HGVSp":null }
        this.conseq = conseq // "missense_variant"
        this.pathogenicity = pathogenicity // {"sift":"deleterious(0)", polyphen:"probably_damaging(0.959)"}
    }

    get FinalHGVS(){
        const reverse_priority = ["HGVSg", "HGVSc", "HGVSp"]
        let HGVS = null
        reverse_priority.forEach( function( tag ){
            if( this.changes.hasOwnproperty(tag) && this.changes[tag] !== null ){
                HGVS = this.changes[tag]
            }
        });
        return HGVS
    }

    static fromJSON(hash){
        const subject = VariantAnnotSubject.fromJSON(hash.subject)
        let obj = new VariantAnnot(subject)
        Object.keys(hash).forEach(function (key) {
            if(key != "subject"){
                obj[key] = hash[key]
            }
        })
        return obj
    }
}

class VariantPopulation {
    constructor(AF, name, source){
        this.AF = AF // 0.000118
        this.name = name // "EAS"
        this.source = source // "ExAC"
    }

    static fromJSON(hash){
        return new VariantPopulation(hash.AF, hash.name, hash.source)
    }
}

class VariantSupport {
    constructor(qual=null, filters=["PASS"], libraries=[], source="unknown"){
        this.filters = filters
        this.qual = qual
        this.libraries = libraries
        this.source = source
    }

    get AF(){
        let alt_depth = 0
        let total_depth = 0
        this.libraries.forEach(function(curr_lib) {
            total_depth += curr_lib.depth
            alt_depth += curr_lib.alt_depth
        });
        return alt_depth/total_depth
    }

    get depth(){
        let total_depth = 0
        this.libraries.forEach(function(curr_lib) {
            total_depth += curr_lib.depth
        })
        return total_depth
    }

    get alt_depth(){
        let alt_depth = 0
        this.libraries.forEach(function(curr_lib) {
            alt_depth += curr_lib.alt_depth
        })
        return alt_depth
    }

    static fromJSON(hash){
        let obj = new VariantSupport()
        Object.keys(hash).forEach(function (key) {
            if(hash[key] !== null){
                obj[key] = hash[key]
            }
        })
        return obj
    }
}

class VariantCoord {
    constructor(ref, alt, region, pos, assembly){
        if( region.toLowerCase().startsWith("chr") ){
            region = region.substr(3)
        }
        this.alt = alt  // str not an array
        this.pos = pos // 15340455
        this.ref = ref
        this.region = region // "chr1"
        this.assembly = assembly // "GRCh38.p12"
    }

    get cpltPos() {
        return this.region + ":" + this.pos + "-" + (this.pos + Math.max(0, this.ref.length -1))
    }

  getName(){
    return this.cpltPos + "=" + this.ref + "/" + this.alt
  }

    refStart() {
        let start = this.pos
        if( [".", "-"].indexOf(this.ref) != -1 ){
            start -= 0.5
        }
        return start
    }

    refEnd() {
        let end = this.pos
        if( [".", "-"].indexOf(this.ref) != -1 ){
            end -= 0.5
        } else {
            end += this.ref.length - 1
        }
        return end
    }

    static fromJSON(hash){
        let obj = new VariantCoord(hash.ref, hash.alt, hash.region, hash.pos)
        if(hash.hasOwnProperty("assembly")){
            obj.assembly = hash.assembly
        }
        return obj
    }
}

class VariantCuration {
    constructor(status, author=null, description=null, date=Date.now()){
        this.status = status
        this.date = date
        this.author = author
        this.description = description
    }

    static fromJSON(hash){
        let obj = new VariantCuration(hash.status)
        Object.keys(hash).forEach(function (key) {
            obj[key] = hash[key]
        })
        return obj
    }
}


class Variant {
    constructor( coord, supports, annot, pop_AF, xref, curations ){
        this.coord = coord
        this.supports = supports
        this.pop_AF = pop_AF // []
        this.annot = annot
        this.xref = xref
        this.curations = curations
    }

    getName() {
          return this.coord.getName()
      }

    get MAF_pop() {
        let MAF_pop = null
        if( this.pop_AF.length > 0 ){
            MAF_pop = this.pop_AF[0]
            this.pop_AF.forEach(function(curr_pop) {
                if( curr_pop.AF > MAF_pop.AF ){
                    MAF_pop = curr_pop
                }
            })
        }
        return MAF_pop
    }

    get filters(){
        const self = this
        let filter_tags = []
        this.sources.forEach(function(curr_src){
            self.supports[curr_src].filters.forEach(function(curr_filter_tag){
                filter_tags.push(curr_filter_tag)
            })
        })
        filter_tags = Array.from(new Set(filter_tags))
        if(filter_tags.length > 1){  // Remove PASS coming from a support if others support tag as not pass
            filter_tags = filter_tags.filter(function(elt){ return elt != "PASS" })
        }
        return filter_tags
    }

    get sources(){
        let sources = Object.keys(this.supports)
        sources.sort()
        return sources
    }

    tagValid( filter ){
        return filter.eval(this.filters)
    }

    static fromJSON(variant){
        const coord = VariantCoord.fromJSON(variant.coord)
        // Supports
        const supports = {}
        if(variant.hasOwnProperty("supports")){
            variant.supports.forEach(function(curr_support){
                if(supports.hasOwnProperty(curr_support.source)){
                    throw `The source already exist '{curr_support.source}'`
                }
                const support = VariantSupport.fromJSON(curr_support)
                supports[support.source] = support
            })
        }
        // Annotations
        const annot = []
        if(variant.hasOwnProperty("annot")){
            variant.annot.forEach(function(curr_annot){
                annot.push(VariantAnnot.fromJSON(curr_annot))
            })
        }
        // Populations
        const pop_AF = []
        if(variant.hasOwnProperty("pop_AF")){
            variant.pop_AF.forEach(function(curr_pop){
                pop_AF.push(VariantPopulation.fromJSON(curr_pop))
            })
        }
        // Curations
        const curations = []
        if(variant.hasOwnProperty("curations")){
            variant.curations.forEach(function(curr_curation){
                curations.push(VariantCuration.fromJSON(curr_curation))
            })
        }
        // Return
        return new Variant(
            coord, supports, annot, pop_AF, variant.xref, curations
        )
    }
}
