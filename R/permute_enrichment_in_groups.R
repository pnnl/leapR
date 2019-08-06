#'permute_enrichment_in_groups
#'
#'
#'

permute_enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5,
                                         mapping_column=NULL, abundance_column=NULL, ntimes=100) {
  # Permute Fisher's exact test and provide an FDR for the real thing

  real = enrichment_in_groups(geneset, targets=targets, background=background, minsize=minsize,
                              mapping_column=mapping_column, abundance_column=abundance_column,
                              method=method)

  # all we need to do is keep track of the number of times we get a better
  #    p-value out of the permutations
  results = rep(0, nrow(real))
  names(results) = rownames(real)

  for (i in 1:ntimes) {
    # permute labels
    unreal = enrichment_in_groups(geneset, targets=targets, background=background, minsize=minsize,
                                  mapping_column=mapping_column, abundance_column=abundance_column,
                                  method=method, randomize=T)

    results[!is.na(unreal[,"pvalue"])] = results[!is.na(unreal[,"pvalue"])] +
      (unreal[which(!is.na(unreal[,"pvalue"])),"pvalue"]<real[which(!is.na(unreal[,"pvalue"])),"pvalue"])
  }
  # return the percentage of times that each functional group
  #        permutation returns a p-value better than the real
  #        functional group
  return(results/ntimes)
}
