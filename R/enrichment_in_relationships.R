#' enrichment_in_relationships
#'
#' enrichment_in_relationships function description is a general way to determine if a pathway
#' is enriched in relationships (interactions, correlation) between its members
#' # access through leapr wrapper
#'
#'

enrichment_in_relationships <- function(geneset, relationships, idmap=NA, tag=NA, mode="original") {
  # for each category in geneset calculates enrichment of within-group
  #     relationships relative to between-group relationships, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference
  results = data.frame(row.names = geneset$names,
                       ingroup_n=rep(NA_real_, length(geneset$names)), ingroupnames=rep(NA_character_, length(geneset$names)), 
                       ingroup_mean=rep(NA_real_, length(geneset$names)), outgroup_n=rep(NA_real_, length(geneset$names)), 
                       zscore=rep(NA_real_, length(geneset$names)), oddsratio=rep(NA_real_, length(geneset$names)), 
                       pvalue=rep(NA_real_, length(geneset$names)), BH_pvalue=rep(NA_real_, length(geneset$names)), 
                       SignedBH_pvalue=rep(NA_real_, length(geneset$names)), background_n=rep(NA_real_, length(geneset$names)),
                       bacground_mean=rep(NA_real_, length(geneset$names)), stringsAsFactors = F)
  
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


    ingroup_ids = grouplist[which(grouplist %in% rownames(relationships))]
    outgroup_ids = rownames(relationships)[which(!rownames(relationships) %in% grouplist)]

    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup = relationships[ingroup_ids, ingroup_ids][upper.tri(relationships[ingroup_ids, ingroup_ids])]

      # pathway members versus all other components
      outgroup = relationships[ingroup_ids, outgroup_ids]

      # non-pathway members versus non-pathway members
      outgroup_2 = relationships[outgroup_ids, outgroup_ids][upper.tri(relationships[outgroup_ids, outgroup_ids])]

    }
    else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits = strsplit(mode, "-")[[1]]
      tag1 = bits[1]
      tag2 = bits[2]

      # filter to give a matrix where rows are tag1 and cols are tag2
      relationshipsa = relationships[grep(tag1, rownames(relationships)), grep(tag2, rownames(relationships))]

      relationshipsa = relationshipsa[which(rownames(relationshipsa) %in% ingroup_ids),
                                      which(colnames(relationshipsa) %in% ingroup_ids)]

      # now match up cognates
      match1 = detag(rownames(relationshipsa))[detag(rownames(relationshipsa)) %in% detag(colnames(relationshipsa))]

      ingroup = sapply(match1, function (p) relationshipsa[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])

      # in the cognate comparison mode this no longer makes complete sense - that is,
      #    the comparison between groups isn't the same
      outgroup = relationships[ingroup_ids, outgroup_ids]

      # non-pathway members versus non-pathway members



      outgroup_2 = relationships[outgroup_ids, outgroup_ids][upper.tri(relationships[outgroup_ids, outgroup_ids])]
      #browser()
    }

    in_mean = mean(unlist(ingroup), na.rm=T)
    out_mean = mean(unlist(outgroup), na.rm=T)

    out_mean_2 = mean(unlist(outgroup_2), na.rm=T)

    pvalue = NA
    pvalue_2 = NA
    if (length(ingroup)>1) {
      pvalue = try(t.test(ingroup, outgroup)$p.value, silent=T);
      if (class(pvalue)=="try-error") pvalue = NA;
      pvalue_2 = try(t.test(ingroup, outgroup_2)$p.value, silent=T)
      if (class(pvalue_2)=="try-error") pvalue_2 = NA
    }

    delta = in_mean - out_mean
    
    results[thisname,"ingroup_n"] = length(ingroup)
    #results[thisname,"ingroupnames"] = in_group_name
    results[thisname,"ingroup_mean"] = in_mean
    results[thisname,"outgroup_n"] = length(outgroup)
    results[thisname,"outgroup_mean"] = out_mean
    #results[thisname,"zscore"] = zscore
    results[thisname,"pvalue"] = pvalue
    #results[thisname,"oddsratio"] = foldx
    
  }

  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  #results[,"BH_pvalue_background"] = p.adjust(results[,"pvalue_background"], method="BH")
  return(results)
}
