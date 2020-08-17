#' enrichment_in_abundance
#'
#' Enrichment in abundance calculates enrichment in pathways by the difference in abundance of the pathway members.
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param Is a \emph{mxn} matrix of abundance data (protein, gene, etc.), with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param mapping_column Is an optional character string, a column name of \code{abundance}, for specifying gene mapping (e.g. for phosphoproteomics data). Defaults to NULL
#' @param abundance_column Is a character vector composed of column names from \code{abundance}, that identifies which conditions should be tested.
#' @param fdr A numerical value which specifies how many times to randomly sample genes, defaults to 0.
#' @param matchset defaults to NULL
#' @param sample_comparison Is a character vector composed of column names from \code{abundance}, that identifies which conditions should be compared. Defaults to NULL
#' @param min_p_threshold Is a numeric value, a lower p-value threshold for returning results, defaults to NULL.
#' @param a prefix tag (optional) for specifying source of data default is NA
#' @param sample_n If set specifies the number of randomizations to perform to calculate an FDR [warning: time intensive], default is NA
#' 
#' @details Case 1: Single condition compatison (sample_comparison = NULL). This usage applies a T test to 
#' @details assess the significance of the abundance of pathway members relative to
#' @details genes/proteins for a condition or conditions that are not in the pathway for a single condition
#' @details Case 2: Case/control comparison (sample_comparison is specified. Compare the abundance of
#' @details the members in the pathway for one set of conditions (abundance_column - which can be a vector) versus
#' @details another set of conditions (sample_comparison - also a vector). 
#' @details 
#' @details These are calculated for every pathway in geneset.
#' @details 
#' @details The significance of this is based on the hypothesis that pathway activity is related to abdundance,
#' @details so that if the members of a praticular pathway are higher abundance than everything else, then it
#' @details can be considered more 'active' (Case 1). Or in Case 2, if the abundance in pathway members is higher
#' @details in one set of conditions (case) relative to in another set of conditions (controls) then it is more
#' @details 'active' in that case.
#' @details 
#' @details Note that the use of the T test to statistically test the difference between the two sets of abundance
#' @details measures is very simple and may be subject to caveats (e.g. low number of comparisons, not-normally 
#' @details distributed abundance values, etc.)
#' @details 
#' @details Providing a number for \code{sample_n} specifies the number of times the pathway membership will be
#' @details randomized to calculate a fasle discovery rate for the enrichment. This is a time-intensive option
#' @details since it requires running the enrichment sample_n+1 times.
#' @details 
#' @details The interpretation of results from Case 1 to Case 2 should vary considerably.
#'
#' @examples
#' dontrun{
#'         library(LEAP)
#'
#'         # read in the example abundance data
#'         data("protdata")
#'
#'         # read in the pathways
#'         data("ncipid")
#'
#'         # read in the patient groups
#'         data("short_list")
#'         data("longlist")
#'
#'         #in this example we lump a bunch of patients together (the 'short survivors') and compare them to another group (the 'long survivors')
#'         protdata.enrichment.svl = enrichment_in_abundance(ncipid, protdata, abundance_column=shortlist, sample_comparison=longlist)
#'
#'         #another application is to compare just one patient against another (this would be the equivalent of comparing one time point to another)
#'         protdata.enrichment.svl.ovo = enrichment_in_abundance(ncipid, protdata, abundance_column=shortlist[1], sample_comparison=longlist[1])
#'
#' }
#'
#'
#'

enrichment_in_abundance <- function(geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                    fdr=0, matchset=NULL, longform=F, sample_comparison=NULL,
                                    min_p_threshold=NULL, tag=NA, sample_n=NULL) {
  # for each category in geneset calculates the abundance level
  #     of the genes/proteins in the category versus those
  #     not in the category and calculate a pvalue based on a
  #     two sided t test
  # If the sample_comparison variable is set then do a comparison between this
  #     abundance and sample comparison set which is vector of valid column (sample) ids
  results = data.frame(row.names = geneset$names,
                       ingroup_n=rep(NA_real_, length(geneset$names)), ingroupnames=rep(NA_character_, length(geneset$names)), 
                       ingroup_mean=rep(NA_real_, length(geneset$names)), outgroup_n=rep(NA_real_, length(geneset$names)), 
                       zscore=rep(NA_real_, length(geneset$names)), oddsratio=rep(NA_real_, length(geneset$names)), 
                       pvalue=rep(NA_real_, length(geneset$names)), BH_pvalue=rep(NA_real_, length(geneset$names)), 
                       SignedBH_pvalue=rep(NA_real_, length(geneset$names)), background_n=rep(NA_real_, length(geneset$names)),
                       bacground_mean=rep(NA_real_, length(geneset$names)), stringsAsFactors = F)
  
  
  if (!is.null(mapping_column)) groupnames = unique(abundance[,mapping_column])

  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    if (!is.null(matchset) && matchset != thisname) next

    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    #cat(thisname, thissize, "\n")
    grouplist = geneset$matrix[i,1:thissize]
    if (!is.na(tag)) grouplist = sapply(grouplist, function (n) paste(tag, n, sep="_"))

    if (!is.null(mapping_column)) {
      ingroupnames = grouplist[which(grouplist %in% groupnames)]
      outgroupnames = groupnames[which(!groupnames %in% grouplist)]

      if (!is.null(sample_n)) {
        if (sample_n > length(ingroupnames) || sample_n > length(outgroupnames)) {
          next
        }
        #cat(sample_n, length(ingroupnames), length(outgroupnames), "\n")
        ingroupnames = sample(ingroupnames, sample_n)
        outgroupnames = sample(outgroupnames, sample_n)
      }

      ingroup = abundance[which(abundance[,mapping_column] %in% ingroupnames),
                          abundance_column[which(abundance_column %in% colnames(abundance[,2:ncol(abundance)]))]]

      if (!is.null(sample_comparison)) outgroup = abundance[which(abundance[,mapping_column] %in% ingroupnames),
                                                            sample_comparison[which(sample_comparison %in% colnames(abundance[,2:ncol(abundance)]))]]
      else outgroup = abundance[which(abundance[,mapping_column] %in% outgroupnames), abundance_column]
    }
    else {
      # I changed this to make it easier to use- may break old code!
      # ingroup = abundance[grouplist[which(grouplist %in% names(abundance))]]
      # outgroup = abundance[which(!names(abundance) %in% grouplist)]
      ingroupnames = grouplist[which(grouplist %in% rownames(abundance))]
      ingroup = unlist(abundance[ingroupnames,
                                 abundance_column[which(abundance_column %in% colnames(abundance))]])

      if (!is.null(sample_comparison)) outgroup = abundance[which(rownames(abundance) %in% grouplist),
                                                            sample_comparison[which(sample_comparison %in% colnames(abundance))]]
      else outgroup = abundance[which(!rownames(abundance) %in% grouplist),
                                abundance_column[which(abundance_column %in% colnames(abundance))]]
    }
    #cat(length(ingroup), length(outgroup), "\n")
    in_mean = mean(unlist(ingroup), na.rm=T)
    out_mean = mean(unlist(outgroup), na.rm=T)
    pvalue = NA
    if (length(ingroup)>1) {
      pvalue = try(t.test(unlist(ingroup), unlist(outgroup))$p.value, silent=T);
      if (class(pvalue)=="try-error") pvalue = NA;

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
        in_mean = mean(ingroup, na.rm=T)
        out_mean = mean(outgroup, na.rm=T)
        delta_r = in_mean - out_mean
        background = c(background, delta_r)
      }
      pvalue = sum(abs(background)>abs(delta))/length(background)
    }
    
    # changing the order of output columns to match with the results from the other enrichment methods
    
    ingroupnames = paste(ingroupnames, collapse = ", ")
    results[thisname,"ingroup_n"] = length(unlist(ingroup))
    results[thisname,"ingroupnames"] = ingroupnames
    results[thisname,"ingroup_mean"] = in_mean
    results[thisname,"outgroup_n"] = length(unlist(outgroup))
    results[thisname,"outgroup_mean"] = out_mean
    results[thisname,"pvalue"] = pvalue
    
    #question : do we want to calculate an oddsratio for this too?
  }
  #update
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  
  #if (!is.null(matchset)) {
  #  results = results[matchset,]
  #  if (longform==T) results = list(results, ingroup, outgroup)
  #}
  return(results)
}
