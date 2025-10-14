### using enrichment_wrapper function
protdata.enrichment.correlation <- leapR::leapR(
  geneset = ncipid,
  enrichment_method = "correlation_enrichment",
  assay_name = "proteomics",
  eset = pset
)

sigs <- subset(protdata.enrichment.correlation,BH_pvalue < 0.05)

expect_equal(nrow(sigs), 71)
