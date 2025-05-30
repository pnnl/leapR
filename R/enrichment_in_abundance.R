#' enrichment_in_abundance
#'
#' Enrichment in abundance calculates enrichment in pathways by the difference in abundance of the pathway members.
# # access through leapr wrapper
#' @import stats
#' @param geneset Gene set to calculate enrichmnet
#' @param abundance Molecular abundance matrix
#' @param mapping_column Column to use to map identifiers
#' @param abundance_column  Columns to use to quantify abundance
#' @param fdr number of times to sample for FDR value
#' @param matchset Name of a set to use for enrichment
#' @param sample_comparison list of samples to use as comparison
#' @param min_p_threshold  Only include p-values lower than this
#' @param tag tag to use for name of group
#' @param sample_n size of sample to use
#' @param silence_try_errors set to true to silence try errors
#' @return data frame of enrichment result

enrichment_in_abundance <-
  function(geneset,
           abundance,
           mapping_column = NULL,
           abundance_column = NULL,
           fdr = 0,
           matchset = NULL,
           sample_comparison = NULL,
           min_p_threshold = NULL,
           tag = NA,
           sample_n = NULL,
           silence_try_errors = TRUE) {
    # for each category in geneset calculates the abundance level
    #     of the genes/proteins in the category versus those
    #     not in the category and calculate a pvalue based on a
    #     two sided t test
    # If the sample_comparison variable is set then do a comparison between this
    #     abundance and sample comparison set which is vector of valid column (sample) ids
    results = data.frame(
      row.names = geneset$names,
      ingroup_n = rep(NA_real_, length(geneset$names)),
      ingroupnames = rep(NA_character_, length(geneset$names)),
      ingroup_mean = rep(NA_real_, length(geneset$names)),
      outgroup_n = rep(NA_real_, length(geneset$names)),
      outgroup_mean = rep(NA_real_, length(geneset$names)),
      zscore = rep(NA_real_, length(geneset$names)),
      oddsratio = rep(NA_real_, length(geneset$names)),
      pvalue = rep(NA_real_, length(geneset$names)),
      BH_pvalue = rep(NA_real_, length(geneset$names)),
      SignedBH_pvalue = rep(NA_real_, length(geneset$names)),
      background_n = rep(NA_real_, length(geneset$names)),
      background_mean = rep(NA_real_, length(geneset$names)),
      stringsAsFactors = FALSE
    )
    
    
    if (!is.null(mapping_column))
      groupnames = unique(abundance[, mapping_column])
    
    for (i in 1:length(geneset$names)) {
      thisname = geneset$names[i]
      if (!is.null(matchset) && matchset != thisname)
        next
      
      thissize = geneset$size[i]
      thisdesc = geneset$desc[i]
      #cat(thisname, thissize, "\n")
      grouplist = geneset$matrix[i, 1:thissize]
      if (!is.na(tag))
        grouplist = sapply(grouplist, function (n)
          paste(tag, n, sep = "_"))
      
      if (!is.null(mapping_column)) {
        ingroupnames = grouplist[which(grouplist %in% groupnames)]
        outgroupnames = groupnames[which(!groupnames %in% grouplist)]
        
        if (!is.null(sample_n)) {
          if (sample_n > length(ingroupnames) ||
              sample_n > length(outgroupnames)) {
            next
          }
          #cat(sample_n, length(ingroupnames), length(outgroupnames), "\n")
          ingroupnames = sample(ingroupnames, sample_n)
          outgroupnames = sample(outgroupnames, sample_n)
        }
        
        ingroup = abundance[which(abundance[, mapping_column] %in% ingroupnames),
                            abundance_column[which(abundance_column %in% colnames(abundance[, 2:ncol(abundance)]))]]
        
        if (!is.null(sample_comparison))
          outgroup = abundance[which(abundance[, mapping_column] %in% ingroupnames),
                               sample_comparison[which(sample_comparison %in% colnames(abundance[, 2:ncol(abundance)]))]]
        else
          outgroup = abundance[which(abundance[, mapping_column] %in% outgroupnames), abundance_column]
      }
      else {
        # I changed this to make it easier to use- may break old code!
        # ingroup = abundance[grouplist[which(grouplist %in% names(abundance))]]
        # outgroup = abundance[which(!names(abundance) %in% grouplist)]
        ingroupnames = grouplist[which(grouplist %in% rownames(abundance))]
        ingroup = unlist(abundance[ingroupnames,
                                   abundance_column[which(abundance_column %in% colnames(abundance))]])
        
        if (!is.null(sample_comparison))
          outgroup = abundance[ingroupnames,
                               sample_comparison[which(sample_comparison %in% colnames(abundance))]]
        else
          outgroup = abundance[which(!rownames(abundance) %in% grouplist),
                               abundance_column[which(abundance_column %in% colnames(abundance))]]
      }
      #cat(length(ingroup), length(outgroup), "\n")
      in_mean = mean(unlist(ingroup), na.rm = TRUE)
      out_mean = mean(unlist(outgroup), na.rm = TRUE)
      pvalue = NA
      if (length(ingroup) > 1) {
        pvalue = try(t.test(unlist(ingroup), unlist(outgroup))$p.value, silent =
                       silence_try_errors)
        ;
        if (is(pvalue,"try-error"))
          pvalue = NA
        
        
        # we step through all components to calculate the number
        #    that are above/below the threshold
        # if (!is.null(min_p_threshold)) {
        #   count_above = 0
        #   count_below = 0
        #   # this works with sample_comparison
        #   for (this_bit in rownames(ingroup)) {
        #     petevalue = try(t.test(ingroup[this_bit,], outgroup[this_bit,])$p.value, silent=F);
        #     if (class(petevalue)=="try-error" || is.na(petevalue)) {
        #       petevalue = NA;
        #     }
        #     else {
        #       if (petevalue > min_p_threshold) count_above = count_above + 1
        #       if (petevalue <= min_p_threshold) count_below = count_below + 1
        #     }
        #   }
        #   results[thisname,"count_above"] = count_above
        #   results[thisname,"count_below"] = count_below
        # }
      }
      
      delta = in_mean - out_mean
      if (fdr) {
        background = c()
        abundances = c(ingroup, outgroup)
        for (i in 1:fdr) {
          # randomly sample genes for fdr times
          ingroup = sample(1:length(abundances), length(ingroup))
          outgroup = which(!1:length(abundances) %in% ingroup)
          ingroup = abundances[ingroup]
          outgroup = abundances[outgroup]
          in_mean = mean(ingroup, na.rm = TRUE)
          out_mean = mean(outgroup, na.rm = TRUE)
          delta_r = in_mean - out_mean
          background = c(background, delta_r)
        }
        pvalue = sum(abs(background) > abs(delta)) / length(background)
      }
      
      # changing the order of output columns to match with the results from the other enrichment methods
      
      ingroupnames = paste(ingroupnames, collapse = ", ")
      results[thisname, "ingroup_n"] = length(unlist(ingroup))
      results[thisname, "ingroupnames"] = ingroupnames
      results[thisname, "ingroup_mean"] = in_mean
      results[thisname, "outgroup_n"] = length(unlist(outgroup))
      results[thisname, "outgroup_mean"] = out_mean
      results[thisname, "pvalue"] = pvalue
      results[thisname, "oddsratio"] = delta
      
      #question : do we want to calculate an oddsratio for this too?
      # answer: yes, but for now we'll use the mean of the ingroup compared with the
      #        distribution of the background as a zscore
      zscore = (out_mean - in_mean) / sd(unlist(outgroup), na.rm = TRUE)
      #browser()
      results[thisname, "zscore"] = zscore
    }
    #update
    results[, "BH_pvalue"] = p.adjust(results[, "pvalue"], method = "BH")
    
    if (!is.null(min_p_threshold)) {
      results = results[results$pvalue < min_p_threshold,]
      return(results)
    } else{
      return(results)
    }
  }
