#' correlation_enrichment
#'
#' # calculate enrichment in correlation between pathway members
#' # access through leapr wrapper
#' @importFrom stats sd
#' @importFrom stats p.adjust
#' @import SummarizedExperiment
#' @param geneset Geneset list
#' @param eset a SummarizedExperiment object
#' @param assay_name name of assay
#' @param mapping_column Column to use to map identifiers, if not rownames
#' @return list of enrichment statistic table and correlation matrix

correlation_enrichment <- function(geneset, eset, assay_name, mapping_column=NA) {#}, tag=NA) {
  ##which genes are in geneset
  allgenes = unique(unlist(as.list(geneset$matrix)))

#  if (!is.na(tag)) allgenes = sapply(allgenes, function (n) paste(tag, n, sep="_"))

  # fixme: add support for an id column
  ids = rownames(eset)
  cols = 1:ncol(eset)
  if (!is.na(mapping_column)) {
    map <- SummarizedExperiment::rowData(eset)[,mapping_column]
    names(map) <- rownames(eset)
  }else{
    map = rownames(eset)
    names(map) <- rownames(eset)
  }
#  if (is.na(mapping_column)) {
  allgenes_present = names(map)[which(map%in%allgenes)]

  # NOTE: assumes that the mapping column is 1 and everything else
  #       is valid data- which may not be the case
  cols = 1:ncol(eset)
  #map = eset
  
  #allgenes_present = allgenes[which(allgenes %in% )]
  allgenes_cor = cor(t(SummarizedExperiment::assay(eset, assay_name)[allgenes_present,cols]), use = "p")

  return(list(enrichment = enrichment_in_relationships(geneset, allgenes_cor, idmap = map),
              corrmat = allgenes_cor))
}
