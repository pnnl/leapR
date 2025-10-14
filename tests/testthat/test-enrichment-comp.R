test_enrichment_comp <- function(){

  protdata.enrichment.svl <- leapR::leapR(
    geneset = ncipid,
    enrichment_method = "enrichment_comparison",
    eset = pset,
    assay_name = "proteomics",
    primary_columns = shortlist,
    secondary_columns = longlist
  )

  sigs <- subset(protdata.enrichment.svl,BH_pvalue < 0.05)

  expect_equal(nrow(sigs), 50)
}


