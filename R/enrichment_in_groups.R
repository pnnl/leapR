#' enrichment_in_groups
#'
#' Calculate the enrichment in pathways using Fisher's exact or Kolgmorov-Smirnov test
#'
#' @param genesets is a GeneSet object for pathway annotation
#' @param targets defaults to NULL
#' @param background Is a \emph{mxn} matrix of gene expression data, with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param method A character string specifying which enrichment method to use fishers or ks for Kolgmorov-Smirnov, defaults to 'fishers'.
#' @param minsize Is the minimum size of a gene set that will be considered, defaults to 5
#' @param mapping_column Is an optional character string, a column name of \code{abundance}, for specifying gene mapping (e.g. for phosphoproteomics data). Defaults to NULL
#' @param abundance_column Is a character vector composed of column names from \code{background}, that ...??.
#' @param randomize is a logical, defaults to FALSE
#' 
#' @details  
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#' @details
#'
#' @examples
#' dontrun{
#'        library(LEAP)
#'
#'        # read in the example abundance data
#'        data("protdata")
#'
#'        # read in the pathways
#'        data("ncipid")
#'
#'        #for this example we will construct a list of genes from the expression data to emulate what you might be inputting
#'        genelist = rownames(protdata)[which(protdata[,1]>0.5)]
#'        background = rownames(protdata)
#'
#'        prodata.enrichment.fishers = enrichment_in_groups(ncipid, targets=genelist, background=background)
#'
#'        #in this example we construct some modules from the heirarchical clustering of the data
#'        protdata_naf = as.matrix(protdata)
#'
#'        #hierarchical clustering is not too happy with lots of missing values so we'll do a zero fill on this to get the modules
#'        protdata_naf[which(is.na(protdata_naf))] = 0
#'
#'        #construct the hierarchical clustering using the 'wardD' method, which seems to give more even sized modules
#'        protdata_hc = hclust(dist(protdata_naf), method="ward")
#'
#'        #arbitrarily we'll chop the clusters into 5 modules
#'        modules = cutree(protdata_hc, k=5)
#'
#'        #modules is a named list of values where each value is a module number and the name is the gene name
#'
#'        #To do enrichment for one module (module 1 in this case) do this
#'        protdata.enrichment.fishers.module_1 = enrichment_in_groups(ncipid, targets=names(modules[which(modules==1)]),background=names(modules))
#'
#' }
#'
#'

enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5,
                                 mapping_column=NULL, abundance_column=NULL, randomize=F) {

  resultp = c()
  resultf = c()
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
      
      results[thisname,"ingroup_n"] = mat[1,1]
      results[thisname,"ingroupnames"] = names
      results[thisname,"outgroup_n"] = mat[1,2]
      results[thisname,"pvalue"] = p
      results[thisname,"oddsratio"] = f
      results[thisname,"background_n"] = mat[2,2]
      # putting this here is a bit of a kludge (since it's not actually a mean, it's a count)
      results[thisname,"background_mean"] = mat[2,1]
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

        results[thisname,"ingroup_n"] = in_path
        results[thisname,"ingroupnames"] = in_group_name
        results[thisname,"ingroup_mean"] = mean(in_group, na.rm=T)
        results[thisname,"outgroup_n"] = in_back
        results[thisname,"zscore"] = zscore
        results[thisname,"pvalue"] = p.value
        results[thisname,"oddsratio"] = foldx
      }
    }
  }
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  results[,"SignedBH_pvalue"] = results[,"BH_pvalue"]*sign(results[,"ingroup_mean"])
  return(results)
}
