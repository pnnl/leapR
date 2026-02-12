genelist <- rownames(pset)[which(SummarizedExperiment::assay(pset, "proteomics")[, 1] > 0.5)]

protdata.enrichment.sets.test <- leapR::leapR(
  geneset = ncipid,
  enrichment_method = "enrichment_in_sets",
  eset = pset,
  assay_name = "proteomics",
  targets = genelist
)

nsigs <- subset(protdata.enrichment.sets.test,BH_pvalue < 0.5)


expect_equal(nrow(nsigs), 6)

