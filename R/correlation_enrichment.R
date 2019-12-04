#' correlation_enrichment
#'
#' correlation_enrichment function description is...
#'
#' @param geneset is a list of four vectors, gene names, gene descriptions, gene sizes and a matrix...??
#' @param abundance Is a \emph{mxn} matrix of gene expression data, with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param mapping_column Is a character string, a column name of \code{abundance}, that...??.
#' @param tag default is NA
#'
#' @examples
#' dontrun{
#'         library(mcdeR)
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
#'         #this is using all the patients regardless of survival time
#'         protdata.enrichment.correlation = correlation_enrichment(ncipid, protdata)
#'
#'         #now we'll do two more- one for the short survivors and one for the long survivors
#'         protdata.enrichment.correlation.short = correlation_enrichment(ncipid, protdata[,shortlist])
#'         protdata.enrichment.correlation.long = correlation_enrichment(ncipid, protdata[,longlist])
#' }
#'
#' @export
#'

correlation_enrichment <- function(geneset, abundance, mapping_column=NA, tag=NA) {
  allgenes = unique(unlist(as.list(geneset$matrix)))

  if (!is.na(tag)) allgenes = sapply(allgenes, function (n) paste(tag, n, sep="_"))

  # fixme: add support for an id column
  ids = rownames(abundance)
  cols = 1:ncol(abundance)
  map = NA
  if (!is.na(mapping_column)) {
    allgenes = rownames(abundance)[which(abundance[,1] %in% allgenes)]

    # NOTE: assumes that the mapping column is 1 and everything else
    #       is valid data- which may not be the case
    cols = 2:ncol(abundance)
    map = abundance
  }
  allgenes_present = allgenes[which(allgenes %in% ids)]
  allgenes_cor = cor(t(abundance[allgenes_present,cols]), use="p")

  return(list(enrichment=enrichment_in_relationships(geneset, allgenes_cor, idmap=map, tag=tag),
              corrmat=allgenes_cor))
}
