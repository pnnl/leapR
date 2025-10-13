### using enrichment_wrapper function
protdata.enrichment.correlation <- leapR::leapR(
  geneset = ncipid,
  enrichment_method = "correlation_enrichment",
  assay_name = "proteomics",
  eset = pset
)
