test_enrichment_path <- function(){

  protdata.enrichment.svl <- leapR::leapR(
    geneset = ncipid,
    enrichment_method = "enrichment_in_pathway",
    eset = pset,
    assay_name = "proteomics",
    primary_columns = shortlist,
    secondary_columns = longlist
  )

  sigs <- subset(protdata.enrichment.svl,BH_pvalue < 0.05)

  expect_equal(nrow(sigs), 123)
}


