#' cluster_enrichment
#'
#' Cluster enrichment Run enrichment (Fisher's exact) on clusters (lists of identifier groups)
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param clusters is a list of clusters (gene lists) to calculate enrichment on
#' @param background is a list of genes to serve as the background for enrichment
#' @param sigfilter minimum significance threshold default is .05 
#'
#' @details This function will calculate enrichment (Fisher's exact test for membership overlap) on
#' @details a series of lists of genes, such as from a set of clusters. The results are returned as
#' @details a list of results matrices in the order of the input clusters.
#'
#' @examples
#'dontrun{
#'         library(leapr)
#'
#'         # read in the example transcriptomic data
#'         data("transdata")
#'
#'         # read in the pathways
#'         data("ncipid")
#'         
#'         # for the example we will limit the number of transcripts considered - arbitrarily in this case
#'         transdata = transdata[1:3000,]
#'
#'         # perform heirarchical clustering on the  data
#'         transdata.hc = hclust(dist(transdata), method="ward.D2")
#'         
#'         transdata.hc.clusters = cutree(transdata.hc, k=5)
#'         
#'         #calculates enrichment for each of the clusters individually and returns a list
#'         #   of enrichment results
#'         transdata.hc.enrichment = cluster_enrichment(geneset=ncipid, clusters=transdata.hc.clusters, background=rownames(transdata))
#'         
#' }
#'
#' @export
#'

cluster_enrichment <- function(geneset, clusters, background=NA, sigfilter=0.05) {
  x = length(clusters)
  
  # The default background is all the genes that are in all the clusters
  # the user can submit their own background list if they want to do something different
  if (!all(is.na(background)))  background = sapply(1:x, function (n) unlist(clusters[[n]]))
  
  #this = sapply(1:x, function (i) list(enrichment_in_groups(geneset, clusters[[i]], background)))
  this = sapply(1:x, function (i) list(leapR(geneset=geneset, enrichment_method="enrichment_in_sets", 
                                             targets=clusters[[i]], background=background)))
  
  # if the sigfilter is set we'll only return those functions that have a p-value
  #   lower than the threshold
  if (is.na(sigfilter)) return(this)
  
  outlist = list()
  for (i in 1:x) {
    these = this[[i]]
    these = these[which(these[,"BH_pvalue"]<sigfilter),]
    outlist = c(outlist, list(these))
  }
  return(outlist)
}
