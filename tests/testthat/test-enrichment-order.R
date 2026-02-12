protdata.enrichment.order <- leapR::leapR(
  geneset = ncipid,
  enrichment_method = "enrichment_in_order",
  eset = pset, # threshold=.5,
  assay_name = "proteomics",
  primary_columns = "TCGA-13-1484"
)


sigs <- subset(protdata.enrichment.order,BH_pvalue < 0.5)

expect_equal(nrow(sigs), 5)

