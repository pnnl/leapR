test_enrichment_comp <- function(){

  protdata.enrichment.svl <- leapR::leapR(
    geneset = ncipid,
    enrichment_method = "enrichment_comparison",
    eset = pset,
    assay_name = "proteomics",
    primary_columns = shortlist,
    secondary_columns = longlist
  )

}
