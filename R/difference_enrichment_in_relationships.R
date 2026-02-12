#' difference_enrichment_in_relationships
#'
#' difference_enrichment_in_relationships function description calculates
#' the difference in two matrices, designed to be symmetric, that each
#' represent a relationship between features. e.g. correlation
#' a t-test is assessed between the two sets of values for each gene in a gene
#' set
#'
#' @param geneset is a list of pathways
#' @param relationships1 is a matrix of pairwise relationships between genes
#' @param relationships2 is a second matrix of pairwise relationships between
#' genes
#'
#' @noRd

difference_enrichment_in_relationships <- function(geneset,relationships1,
                                                   relationships2, idmap = NA) {
  stopifnot(is(geneset,'geneset_data'))
  results2 <- vapply(seq_len(length(geneset$names)), function(i){
    thisname <- geneset$names[i]
    thissize <- geneset$size[i]
    thisdesc <- geneset$desc[i]
    grouplist <- geneset$matrix[i, seq_len(thissize)]
    if (!all(is.na(idmap))) {
      grouplist <- names(idmap)[which(idmap %in% grouplist)]
    }

    # do these separately since the matrices could be different
    ingroup1_ids <- grouplist[which(grouplist %in% rownames(relationships1))]
    ingroup2_ids <- grouplist[which(grouplist %in% rownames(relationships2))]

    ingroup1 <- relationships1[ingroup1_ids, ingroup1_ids, drop = FALSE]
    ingroup1 <- ingroup1[upper.tri(relationships1[ingroup1_ids,
                                                   ingroup1_ids])]
    ingroup2 <- relationships2[ingroup2_ids, ingroup2_ids, drop = FALSE]
    ingroup2 <- ingroup2[upper.tri(relationships2[ingroup2_ids,
                                                    ingroup2_ids])]
    in_mean1 <- mean(unlist(ingroup1), na.rm = TRUE)
    in_mean2 <- mean(unlist(ingroup2), na.rm = TRUE)
    pvalue <- NA
    if (length(ingroup1) > 1 && length(ingroup2) > 1) {
      pvalue <- try(t.test(ingroup1, ingroup2)$p.value, silent = TRUE)
      if (is(pvalue, "try-error")) pvalue <- NA
    }

    res <- c(ingroup_n = length(ingroup1), ingroupnames = "",
             ingroup_mean = in_mean1,
             outgroup_n = length(ingroup2), outgroup_mean = in_mean2,
             zscore = NA,
             oddsratio = in_mean1 - in_mean2, pvalue = pvalue,
             BH_pvalue = NA, SignedBH_pvalue = NA,
             background_n = NA, background_mean = NA)

    return(res)

  },c(ingroup_n = numeric(1), ingroupnames = "", ingroup_mean = numeric(1),
      outgroup_n = numeric(1), outgroup_mean = numeric(1), zscore = numeric(1),
      oddsratio = numeric(1), pvalue = numeric(1),
      BH_pvalue = numeric(1), SignedBH_pvalue = numeric(1),
      background_n = numeric(1), background_mean = numeric(1)))

  results <- t(results2)
  rownames(results) <- geneset$names
  results <- as.data.frame(results)
  return(results)
}
