#' difference_enrichment_in_relationships
#'
#' # access through leapr wrapper
#'

difference_enrichment_in_relationships <- function(geneset, relationships1, relationships2, idmap=NA, tag=NA,
                                                   mode="original") {
  # for each category in geneset calculates enrichment of
  #     relationships in matrix 1 versus those in matrix 2, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference
  results = matrix(nrow=length(geneset$names), ncol=6)
  rownames(results) = geneset$names
  colnames(results) = c("group_mean_1", "group_mean_2", "ingroup1_n",
                        "ingroup2_n", "pvalue", "BH_pvalue")

  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]

    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    #cat(thisname, thissize, "\n")
    grouplist = geneset$matrix[i,1:thissize]

    if (!is.na(tag)) grouplist = sapply(grouplist, function (n) paste(tag, n, sep="_"))

    if (!is.na(idmap)) {
      grouplist = rownames(idmap)[which(idmap[,1] %in% grouplist)]
    }

    # do these separately since the matrices could be different
    ingroup1_ids = grouplist[which(grouplist %in% rownames(relationships1))]
    ingroup2_ids = grouplist[which(grouplist %in% rownames(relationships2))]

    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup1 = relationships1[ingroup1_ids, ingroup1_ids][upper.tri(relationships1[ingroup1_ids, ingroup1_ids])]
      ingroup2 = relationships2[ingroup2_ids, ingroup2_ids][upper.tri(relationships2[ingroup2_ids, ingroup2_ids])]
    }
    else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits = strsplit(mode, "-")[[1]]
      tag1 = bits[1]
      tag2 = bits[2]

      # filter to give a matrix where rows are tag1 and cols are tag2
      relationships1a = relationships1[grep(tag1, rownames(relationships1)), grep(tag2, rownames(relationships1))]
      relationships2a = relationships2[grep(tag1, rownames(relationships2)), grep(tag2, rownames(relationships2))]

      relationships1a = relationships1a[which(rownames(relationships1a) %in% ingroup1_ids),
                                        which(colnames(relationships1a) %in% ingroup1_ids)]
      relationships2a = relationships2a[which(rownames(relationships2a) %in% ingroup2_ids),
                                        which(colnames(relationships2a) %in% ingroup2_ids)]

      # now match up cognates
      match1 = detag(rownames(relationships1a))[detag(rownames(relationships1a)) %in% detag(colnames(relationships1a))]
      match2 = detag(rownames(relationships2a))[detag(rownames(relationships2a)) %in% detag(colnames(relationships2a))]

      ingroup1 = sapply(match1, function (p) relationships1a[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])
      ingroup2 = sapply(match2, function (p) relationships2a[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])
      #browser()
    }

    in_mean1 = mean(unlist(ingroup1), na.rm=T)
    in_mean2 = mean(unlist(ingroup2), na.rm=T)

    pvalue = NA

    if (length(ingroup1)>1 && length(ingroup2)>1) {
      pvalue = try(t.test(ingroup1, ingroup2)$p.value, silent=T);
      if (class(pvalue)=="try-error") pvalue = NA;
    }

    this = c(in_mean1, in_mean2,
             length(ingroup1), length(ingroup2),
             pvalue, NA)

    results[thisname,] = this
  }

  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  return(results)
}
