
function displayKnownVariants(variant){
    let display = getVariantLink(variant.id)
    if(variant.HGVSp != ""){
        display += " " + variant.gene + ":" + variant.HGVSp
    } else if(variant.HGVSc != ""){
        display += " " + variant.gene + ":" + variant.HGVSc
    } else if(variant.gene != ""){
        display += " " + variant.gene
    }
    return display
}

function displayAnnotations(annot){
    let display = ""
    if(annot.SYMBOL != ""){
        display = annot.SYMBOL
        if(annot.Feature_type == "Transcript"){
            let start_type = "exon"
            let start_pos = null
            if(annot.start_EXON === null){
                start_type = "intron"
                start_pos = annot.start_INTRON
            } else {
                start_pos = annot.start_EXON
            }
            let end_type = "exon"
            let end_pos = null
            if(annot.end_EXON === null){
                end_type = "intron"
                end_pos = annot.end_INTRON
            } else {
                end_pos = annot.end_EXON
            }
            if(start_type == end_type && start_pos == end_pos){
                display += `: ${annot.Feature} on ${start_type} ${start_pos}`
            } else {
                display += `: ${annot.Feature} from ${start_type} ${start_pos} to ${end_type} ${end_pos}`
            }

        }
    }
    return display
}

Vue.component('shallows-analysis-table', {
    props: {
        analysis: {
            type: Object,
            required: true
        },
        columns: {
            type: Array,
            default: function(){
                return [
                    {"title": "Region", "value": function(entry, col){ return entry.reference }},
                    {"title": "Start", "value": function(entry, col){ return entry.start }},
                    {"title": "End", "value": function(entry, col){ return entry.end }},
                    {"title": "Length", "value": function(entry, col){ return entry.end - entry.start + 1 }},
                    {
                        "title": "Overlapped regions",
                        "value": function(entry, col){
                            return entry.annotations.map(displayAnnotations).join("<br />")
                        },
                        "is_html": true
                    },
                    {
                        "title": "Overlapped known variants",
                        "value": function(entry, col){
                            return entry.known_variants.map(displayKnownVariants).join("<br />")
                        },
                        "is_html": true
                    }
                ]
            }
        }
    },
    computed: {
        title: function(){
			return `Regions with depth < ${this.analysis.parameters.min_depth} ${this.analysis.parameters.depth_mode}s`
		}
    },
	template:
		`<dynamic-table
            :data="analysis.results"
            :header="columns"
            :title="title">
        </dynamic-table>`
})
