#' correlation_enrichment
#'
#' correlation_enrichment function description is...
#'
#' @param geneset is...
#' @param abundance is...
#' @param mapping_column default is NA
#' @param tag default is NA
#'
#' @examples
#' dontrun{
#'         library(readr)
#'         data("protdata")
#'         data("shortlist")
#'         data("longlist")
#'
#'         #we have to update the colnames slightly to match what we have in the patient groups
#'         colnames(protdata) = sapply(colnames(protdata), function (r) {a=strsplit(r,"\\.");paste("TCGA",a[[1]][2], a[[1]][3], sep="-")})
#'
#'         #read in the pathways
#'         ncipid = read_gene_sets(gsfile = "/example/NCI_PID_genesymbol_corrected.gmt")
#'
#'         #this is using all the patients regardless of survival time
#'         protdata.enrichment.correlation = correlation_enrichment(ncipid, protdata)
#'
#'         #now we'll do two more- one for the short survivors and one for the long survivors
#'         protdata.enrichment.correlation.short = correlation_enrichment(ncipid, protdata[,shortlist])
#'         protdata.enrichment.correlation.long = correlation_enrichment(ncipid, protdata[,longlist])
#'
#'
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
