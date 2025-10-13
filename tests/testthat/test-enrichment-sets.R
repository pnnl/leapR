genelist <- rownames(pset)[which(SummarizedExperiment::assay(pset, "proteomics")[, 1] > 0.5)]

protdata.enrichment.sets.test <- leapR::leapR(
  geneset = ncipid,
  enrichment_method = "enrichment_in_sets",
  eset = pset,
  assay_name = "proteomics",
  targets = genelist
)
