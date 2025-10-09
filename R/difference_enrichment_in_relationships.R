#' difference_enrichment_in_relationships
#'
#' difference_enrichment_in_relationships function description is...
#'
#' @param geneset is...
#' @param relationships1 is...
#' @param relationships2 is...
#' @param idmap defaults to an NA
#' @param tag defaults to an NA
#'
#'
#' @noRd

difference_enrichment_in_relationships <- function(geneset, relationships1, relationships2,
                                                   idmap = NA, tag = NA,
                                                   mode = "original") {
  # for each category in geneset calculates enrichment of
  #     relationships in matrix 1 versus those in matrix 2, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference

  results2 <- vapply(seq_along(1:length(geneset$names)), function(i){
    thisname <- geneset$names[i]
    thissize <- geneset$size[i]
    thisdesc <- geneset$desc[i]
    grouplist <- geneset$matrix[i, 1:thissize]

    if (!is.na(tag)) grouplist <- vapply(grouplist, function(n) paste(tag, n, sep = "_"), "")

    if (!all(is.na(idmap))) {
      grouplist <- names(idmap)[which(idmap %in% grouplist)]
    }

    # do these separately since the matrices could be different
    ingroup1_ids <- grouplist[which(grouplist %in% rownames(relationships1))]
    ingroup2_ids <- grouplist[which(grouplist %in% rownames(relationships2))]

    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup1 <- relationships1[ingroup1_ids, ingroup1_ids, drop = FALSE][upper.tri(relationships1[ingroup1_ids, ingroup1_ids])]
      ingroup2 <- relationships2[ingroup2_ids, ingroup2_ids, drop = FALSE][upper.tri(relationships2[ingroup2_ids, ingroup2_ids])]
    } else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits <- strsplit(mode, "-")[[1]]
      tag1 <- bits[1]
      tag2 <- bits[2]

      # filter to give a matrix where rows are tag1 and cols are tag2
      relationships1a <- relationships1[grep(tag1, rownames(relationships1)), grep(tag2, rownames(relationships1)), drop = FALSE]
      relationships2a <- relationships2[grep(tag1, rownames(relationships2)), grep(tag2, rownames(relationships2)), drop = FALSE]

      relationships1a <- relationships1a[
        which(rownames(relationships1a) %in% ingroup1_ids),
        which(colnames(relationships1a) %in% ingroup1_ids)
      ]
      relationships2a <- relationships2a[
        which(rownames(relationships2a) %in% ingroup2_ids),
        which(colnames(relationships2a) %in% ingroup2_ids)
      ]

      # now match up cognates
      match1 <- detag(rownames(relationships1a))[detag(rownames(relationships1a)) %in% detag(colnames(relationships1a))]
      match2 <- detag(rownames(relationships2a))[detag(rownames(relationships2a)) %in% detag(colnames(relationships2a))]

      ingroup1 <- vapply(match1, function(p) relationships1a[paste(tag1, p, sep = "_"), paste(tag2, p, sep = "_")], "")
      ingroup2 <- vapply(match2, function(p) relationships2a[paste(tag1, p, sep = "_"), paste(tag2, p, sep = "_")], "")
      # browser()
    }

    in_mean1 <- mean(unlist(ingroup1), na.rm = TRUE)
    in_mean2 <- mean(unlist(ingroup2), na.rm = TRUE)

    pvalue <- NA

    if (length(ingroup1) > 1 && length(ingroup2) > 1) {
      pvalue <- try(t.test(ingroup1, ingroup2)$p.value, silent = TRUE)
      if (is(pvalue, "try-error")) pvalue <- NA
    }

    res <- c(ingroup_n = length(ingroup1), ingroupnames = "", ingroup_mean = in_mean1,
             outgroup_n = length(ingroup2), outgroup_mean = in_mean2, zscore = NA,
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
  results[,'pvalue'] <- as.numeric(results[,'pvalue'])
  results[,'ingroup_n'] <- as.numeric(results[,'ingroup_n'])
  results[,'outgroup_n'] <- as.numeric(results[,'outgroup_n'])
  results[,'oddsratio'] <- as.numeric(results[,'oddsratio'])
  results[,'ingroup_mean'] <- as.numeric(results[,'ingroup_mean'])
  results[,'outgroup_mean'] <- as.numeric(results[,'outgroup_mean'])
  results[,'BH_pvalue'] <- p.adjust(results[,'pvalue'], method = "BH")
  results[,'SignedBH_pvalue'] <- results[,'BH_pvalue'] * sign(as.numeric(results[,'ingroup_mean']))
  return(results)
}
