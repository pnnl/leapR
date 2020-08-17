#' correlation_enrichment
#'
#' correlation_enrichment function description is...
#'
#' @param geneset is GeneSet object for pathway annotation
#' @param abundance Is a \emph{mxn} matrix of abundance data, with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param mapping_column Is a character string, a column name of \code{abundance}, that contains gene names (e.g. for phosphoproteomics). Default is NA
#' @param tag is an optional prefix that specifies the data type. Default is NA
#'
#' @description This function calculates the correlation of pathway members, for each pathway, and calculates enrichment based
#' @description on the distribution of this correlation. This is based on a hypothesis that pathways that are more coordinated
#' @description in activity (i.e. members are more correlated across different conditions) might be more important - or more
#' @description active.
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
