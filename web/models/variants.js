/*!
  * shallowsAnalysisTable component v1.0.0
  * Copyright 2019 IUCT-O
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
