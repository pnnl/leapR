#' correlation_enrichment
#'
#' # calculate enrichment in correlation between pathway members
#' # access through leapr wrapper
#'
#'@export

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
