#' cluster_enrichment
#'
#' Cluster enrichment Run enrichment (Fisher's exact) on clusters (lists of identifier groups)
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param clusters is a list of clusters (gene lists) to calculate enrichment on, generally the result of the `cuttree` function
#' @param background is a list of genes to serve as the background for enrichment
#' @param sigfilter minimum significance threshold default is .05 
#'
#' @details This function will calculate enrichment (Fisher's exact test for membership overlap) on
#' @details a series of lists of genes, such as from a set of clusters. The results are returned as
#' @details a list of results matrices in the order of the input clusters.
#' @importFrom utils read.table
#' @return data frame with enrichment results
#' @export
#' @examples
#'         library(leapR)
#'
#'         # read in the example transcriptomic data
#'         datadir='https://github.com/pnnl/leapR/raw/refs/heads/bioc-submission/csv/'
#'         transdata<-readr::read_csv(paste0(datadir,'transdata.csv.gz'))|>
#'            tibble::column_to_rownames('...1')|>
#'              as.matrix()
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
#'
#'

cluster_enrichment <- function(geneset, clusters, background=NA, sigfilter=0.05) {
  x = length(clusters)
  
  # The default background is all the genes that are in all the clusters
  # the user can submit their own background list if they want to do something different
  if (!all(is.na(background)))  background = sapply(1:x, function (n) unlist(clusters[[n]]))
  
  #this = sapply(1:x, function (i) list(enrichment_in_groups(geneset, clusters[[i]], background)))
  this = sapply(1:x, function (i) list(leapR(geneset=geneset, enrichment_method="enrichment_in_sets", 
                                             targets=clusters[[i]], background=background)))
  
  #print(lapply(this,function(i) nrow(subset(i,as.numeric(BH_pvalue)<sigfilter))))
  # if the sigfilter is set we'll only return those functions that have a p-value
  #   lower than the threshold
  if (is.na(sigfilter)) return(this)
  #print(sigfilter)
  outlist = list()
  for (i in 1:x) {
    these <- this[[i]]
    sigs <- which(these$BH_pvalue<sigfilter)
    #print(sigs)
    these <- these[sigs,]
    #print(these)
    outlist = c(outlist, list(these))
  }
  return(outlist)
}
