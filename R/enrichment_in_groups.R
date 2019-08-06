#'enrichment_in_groups
#'
#'
#'


enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5,
                                 mapping_column=NULL, abundance_column=NULL, randomize=F) {

  resultp = c()
  resultf = c()
  results = data.frame(row.names = geneset$names,
                       in_path=rep(NA_real_, length(geneset$names)), in_path_names=rep(NA_character_, length(geneset$names)), out_path=rep(NA_real_, length(geneset$names)),
                       in_back=rep(NA_real_, length(geneset$names)), out_back=rep(NA_real_, length(geneset$names)), foldx=rep(NA_real_, length(geneset$names)),
                       pvalue=rep(NA_real_, length(geneset$names)), Adjusted_pvalue=rep(NA_real_, length(geneset$names)), Signed_AdjP=rep(NA_real_, length(geneset$names)),
                       stringsAsFactors = F)

  if (method == "ks") colnames(results)[c(3,5)] = c("MeanPath", "Zscore")

  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    grouplist = geneset$matrix[i,1:thissize]
    if (randomize) {
      # choose a random set of genes as this grouplist
      # A disadvantage is that we resample for each functional group rather than
      #   running one set of analyses on a fully scrambled set of functions.
      #   I don't think this should be a huge problem though.
      grouplist = sample(unlist(geneset$matrix), length(grouplist))
    }
    in_back = length(background)

    if (method == "fishers") {
      enr = enrichment_by_fishers(targets, background, grouplist)
      p = enr$fisher$p.value
      f = enr$foldx
      mat = enr$mat
      names = enr$in_path_names

      results[thisname, ] = list(mat[1,1], names, mat[1,2], mat[2,1], mat[2,2], f, p, NA, NA)
    }
    else if (method == "ks") {    #Kolmogorov-Smirnov test
      # in this case "background" must be the continuous variable from which grouplist can be drawn
      backlist = background

      if (is.null(mapping_column)) {
        in_group = background[grouplist[which(grouplist %in% rownames(background))],abundance_column]
        in_group_name = paste(intersect(grouplist, rownames(background)), collapse = ", ")
        backlist = background[,abundance_column]
      }
      else {
        # mapping_column adds the ability to use phospho-type data where the gene name (non-unique) is in the
        #       first column and the rownames are peptide ids
        # unfortunately this means that "background" has to be the whole matrix and abundance_column
        #       has to be specified, which is a bit ugly
        in_group = background[which(background[,mapping_column] %in% grouplist),abundance_column]
        in_group_name = paste(intersect(background[,mapping_column], grouplist), collapse = ", ")
        backlist = background[,abundance_column]
      }

      in_path = length(in_group)


      if (in_path > minsize) {
        in_back = length(backlist)

        enr = try(ks.test(in_group, backlist))
        if (class(enr) == "try-error") {
          enr = NA
          p.value = NA
        }
        else {
          p.value = enr$p.value
        }

        # this expression of foldx might be subject to some weird pathological conditions
        # e.g. one sample has a background that is always negative, another that's positive
        # may pertain to zscore too (although not sure it should)
        #foldx = mean(in_group, na.rm=T)/mean(background, na.rm=T)

        # rank from largest to smallest
        if (is.null(mapping_column)) in_rank = rank(backlist)[grouplist[which(grouplist %in% names(background))]]
        else in_rank = rank(backlist)[which(background[,mapping_column] %in% grouplist)]

        foldx = mean(in_rank, na.rm=T)/length(backlist)

        zscore = (mean(in_group, na.rm=T)-mean(backlist, na.rm=T))/sd(in_group, na.rm=T)

        #padj = enr$p.value*length(geneset$names)
        # c("in_path", "MeanPath", "in_back", "Zscore", "foldx", "pvalue", "Adjusted_pvalue")
        results[thisname, ] = list(in_path, in_group_name, mean(in_group, na.rm=T), in_back, zscore, foldx, p.value, NA, NA)
      }
    }
  }
  results[,"Adjusted_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  results[,"Signed_AdjP"] = results[,"Adjusted_pvalue"]*sign(results[,"in_path"])
  return(results)
}
