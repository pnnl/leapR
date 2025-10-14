#' enrichment_in_relationships
#'
#' enrichment_in_relationships function description is a general way to determine if a pathway
#' is enriched in relationships (interactions, correlation) between its members
#' # access through leapr wrapper
#'
#' @param geneset List of pathways in gmt format
#' @param relationships table of relationship information, e.g. correlation
#' @param idmap list of identifiers to use for mapping, the names of the items should agree with names of features in matrix
#' @param mode allowable values are 'original' or anything else
#' @param silence_try_errors boolean to silence errors
#' @return table of enrichment statistics
#' @importFrom stats sd
#' @importFrom stats p.adjust
#'
enrichment_in_relationships <- function(geneset, relationships, idmap = NA,
                                        # tag = NA,
                                        mode = "original",
                                        silence_try_errors = TRUE) {
  # for each category in geneset calculates enrichment of within-group
  #     relationships relative to between-group relationships, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference

  results2 <- vapply(seq_along(1:length(geneset$names)), function(i){
    thisname <- geneset$names[i]
    thissize <- geneset$size[i]
    thisdesc <- geneset$desc[i]
    grouplist <- geneset$matrix[i, seq_along(1:thissize)]

    if (!all(is.na(idmap))) {
      grouplist <- names(idmap)[which(idmap %in% grouplist)]
    }


    ingroup_ids <- grouplist[which(grouplist %in% rownames(relationships))]
    outgroup_ids <- rownames(relationships)[which(!rownames(relationships) %in% grouplist)]

    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup <- relationships[ingroup_ids, ingroup_ids, drop =FALSE][upper.tri(relationships[ingroup_ids, ingroup_ids])]

      # pathway members versus all other components
      outgroup <- relationships[ingroup_ids, outgroup_ids]

      # non-pathway members versus non-pathway members
      outgroup_2 <- relationships[outgroup_ids, outgroup_ids, drop =FALSE][upper.tri(relationships[outgroup_ids, outgroup_ids])]
    } else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits <- strsplit(mode, "-")[[1]]
      tag1 <- bits[1]
      tag2 <- bits[2]

      # filter to give a matrix where rows are tag1 and cols are tag2
      relationshipsa <- relationships[grep(tag1, rownames(relationships)), grep(tag2, rownames(relationships))]

      relationshipsa <- relationshipsa[
        which(rownames(relationshipsa) %in% ingroup_ids),
        which(colnames(relationshipsa) %in% ingroup_ids)
      ]

      # now match up cognates
      match1 <- detag(rownames(relationshipsa))[detag(rownames(relationshipsa)) %in% detag(colnames(relationshipsa))]

      ingroup <- vapply(match1, function(p) relationshipsa[paste(tag1, p, sep = "_"), paste(tag2, p, sep = "_")], "")

      # in the cognate comparison mode this no longer makes complete sense - that is,
      #    the comparison between groups isn't the same
      outgroup <- relationships[ingroup_ids, outgroup_ids, drop =FALSE]

      # non-pathway members versus non-pathway members
      outgroup_2 <- relationships[outgroup_ids, outgroup_ids, drop =FALSE][upper.tri(relationships[outgroup_ids, outgroup_ids])]
      # browser()
    }

    in_mean <- mean(unlist(ingroup), na.rm = TRUE)
    out_mean <- mean(unlist(outgroup), na.rm = TRUE)

    out_mean_2 <- mean(unlist(outgroup_2), na.rm = TRUE)

    pvalue <- NA
    pvalue_2 <- NA
    if (length(ingroup) > 1) {
      pvalue <- try(t.test(ingroup, outgroup)$p.value, silent = silence_try_errors)
      if (is(pvalue, "try-error")) pvalue <- NA
      pvalue_2 <- try(t.test(ingroup, outgroup_2)$p.value, silent = silence_try_errors)
      if (is(pvalue_2, "try-error")) pvalue_2 <- NA
      if(is.na(pvalue)) pvalue <- pvalue_2
    }

    delta <- in_mean - out_mean

    res <- c(ingroup_n = length(ingroup), ingroupnames = paste(ingroup_ids, collapse = ','),  ingroup_mean = in_mean,
             outgroup_n =  length(outgroup),  outgroup_mean = out_mean, zscore = NA,
             oddratio = delta, pvalue = pvalue,
             BH_pvalue = NA, SignedBH_pvalue = NA,
             background_n = NA, background_mean = NA)
    return(res)

  }, c(ingroup_n = numeric(1), ingroupnames = "", ingroup_mean = numeric(1),
       outgroup_n = numeric(1), outgroup_mean = numeric(1), zscore = numeric(1),
       oddsratio = numeric(1), pvalue = numeric(1),
       BH_pvalue = numeric(1), SignedBH_pvalue = numeric(1),
       background_n = numeric(1), background_mean = numeric(1)))

  results <- t(results2)
  rownames(results) <- geneset$names

  return(results)
}
