#' correlation_enrichment
#'
#' # calculate enrichment in correlation between pathway members
#' # access through leapr wrapper
#' @import stats
#' @import Biobase
#' @param geneset Geneset list
#' @param eset an ExpressionSet object
#' @param mapping_column Column to use to map identifiers, if not rownames
#' @return list of enrichment statistic table and correlation matrix

correlation_enrichment <- function(geneset, eset, mapping_column=NA) {#}, tag=NA) {
  ##which genes are in geneset
  allgenes = unique(unlist(as.list(geneset$matrix)))

#  if (!is.na(tag)) allgenes = sapply(allgenes, function (n) paste(tag, n, sep="_"))

  # fixme: add support for an id column
  ids = rownames(eset)
  cols = 1:ncol(eset)
  if (!is.na(mapping_column)) {
    map <- Biobase::fData(eset)[,mapping_column]
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
  allgenes_cor = cor(t(Biobase::exprs(eset)[allgenes_present,cols]), use = "p")

  return(list(enrichment = enrichment_in_relationships(geneset, allgenes_cor, idmap = map),
              corrmat = allgenes_cor))
}
