#' cluster_enrichment
#'
#' Cluster enrichment Run enrichment (Fisher's exact) on clusters (lists of
#' identifier groups)
#' @param eset is an SummarizedExperiment containing data that is clustered
#' @param assay_name is the name of the assay
#' @param geneset is a GeneSet object for pathway annotation
#' @param clusters is a list of clusters (gene lists) to calculate enrichment
#'  on, generally the result of the `cutree` function
#' @param sigfilter minimum significance threshold default is .05
#'
#' @details This function will calculate enrichment (Fisher's exact test for
#'  membership overlap) on
#' @details a series of lists of genes, such as from a set of clusters. The
#' results are returned as
#' @details a list of results matrices in the order of the input clusters.
#' @importFrom utils read.table
#' @return data frame with enrichment results
#' @export
#' @examples
#'         library(leapR)
#'         library(BiocFileCache)
#'         
#'         path <- tools::R_user_dir("leapR", which = "cache")
#'         bfc <- BiocFileCache(path, ask = FALSE)
#'         
#'         url <- "https://api.figshare.com/v2/file/download/56536214"
#'         tc <- bfcadd(bfc, "tdat", fpath = url)
#'         load(tc)
#'
#'         # read in the pathways
#'         data("ncipid")
#'
#'         # for the example we will limit the number of transcripts considered
#'         #- arbitrarily in this case
#'         transdata <- SummarizedExperiment::assay(tset,'transcriptomics')
#'         transdata[which(is.na(transdata),arr.ind=TRUE)]<-0.0
#'         # perform heirarchical clustering on the  data
#'         transdata.hc <- hclust(dist(transdata), method="ward.D2")
#'
#'         transdata.hc.clusters <- cutree(transdata.hc, k=5)
#'         clust.list <- lapply(seq_len(5), function(x) {
#'            return(names(which(transdata.hc.clusters==x)))})
#'         #calculates enrichment for each of the clusters individually a
#'         #and returns a list of enrichment results
#'         transdata.hc.enrichment <- leapR::cluster_enrichment(eset=tset,
#'                 assay_name='transcriptomics',
#'                 geneset=ncipid,
#'                 clusters=clust.list)
#'
#'
#'

cluster_enrichment <- function(eset, assay_name, geneset,
                               clusters, sigfilter=0.05) {

  stopifnot(is(geneset,'geneset_data'),
            is(eset,'SummarizedExperiment'),
            length(clusters) > 0)

  x <- length(clusters)

  this <- vapply(seq_len(x), function(i) list(leapR(eset = eset,
                                                      assay_name = assay_name,
                                                      geneset = geneset,
                                                      enrichment_method =
                                                        'enrichment_in_sets',
                                                      targets = clusters[[i]])),
                 data.frame(1))

  # if the sigfilter is set we'll only return those functions that
  # have a p-value  lower than the threshold
  if (is.na(sigfilter)) return(this)

  outlist <- vapply(this, function(these){
    sigs <- which(these$BH_pvalue < sigfilter)
    these <- these[sigs,]
    return(list(these))
    }, list(list))

  return(outlist)
}
