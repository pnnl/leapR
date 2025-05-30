#' correlation_comparison_enrichment
#'
#' # internal function to calculate enrichment in differences in correlation between two groups
#' # access through the leapr wrapper
#'
#' @import stats
#' @param geneset pathway to use for enrichment
#' @param abundance abundance matrix
#' @param set1 first set to use
#' @param set2 second set to use
#' @param mapping_column Column to use for id mapping
#' @param tag Tag to add
#' @param mode to use, default is 'original'
#' @return data frame with enrichment results
#' 
correlation_comparison_enrichment <- function(geneset, abundance, set1, set2, mapping_column=NA, tag=NA, mode="original") {
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
  cols1 = colnames(abundance)[which(colnames(abundance) %in% set1)]
  cols2 = colnames(abundance)[which(colnames(abundance) %in% set2)]
  
  allgenes_present = allgenes[which(allgenes %in% ids)]
  abcols1 = abundance[allgenes_present,cols1]
  abcols2 = abundance[allgenes_present,cols2]
  
  abcols1 = abcols1[which(rowSums(!is.na(abcols1))>0),]
  abcols2 = abcols2[which(rowSums(!is.na(abcols2))>0),]
  
  allgenes_cor1 = cor(t(abcols1), use="p")
  allgenes_cor2 = cor(t(abcols2), use="p")

  return(difference_enrichment_in_relationships(geneset, allgenes_cor1, allgenes_cor2, idmap=map, tag=tag, mode=mode))
}
