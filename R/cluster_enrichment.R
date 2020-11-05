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
  x = length(clusters)-1
  if (is.na(background))  background = clusters[[x+1]]
  this = sapply(1:x, function (i) list(enrichment_in_groups(geneset, clusters[[i]], background)))
  outlist = list()
  for (i in 1:x) {
    these = this[[i]][which(this[[i]][,7]<sigfilter),]
    if (length(which(this[[i]][,7]<sigfilter)) == 1) {
      these = as.matrix(t(these))
      rownames(these)[[1]] = rownames(this[[i]])[which(this[[i]][,7]<sigfilter)]
    }
    outlist = c(outlist, list(these))
  }
  return(outlist)
}
