#' correlation_comparison_enrichment
#'
#' # internal function to calculate enrichment in differences in correlation between two groups
#' # access through the leapr wrapper
#'
#' @importFrom stats sd
#' @importFrom stats p.adjust
#' @param geneset pathway to use for enrichment
#' @param eset  SummarizedExperiment with abundance matrix
#' @param assay_name name of assay
#' @param set1 first set to use
#' @param set2 second set to use
#' @param mapping_column Column to use for id mapping within rowData
#' @param tag Tag to add
#' @param mode to use, default is 'original'
#' @return data frame with enrichment results
#' 
correlation_comparison_enrichment <- function(geneset, eset, assay_name,
                                              set1, set2, 
                                              mapping_column=NA, tag=NA, 
                                              mode="original") {
  
  allgenes = unique(unlist(as.list(geneset$matrix)))

  if (!is.na(tag)) allgenes = sapply(allgenes, function(n) paste(tag, n, sep = "_"))

  # fixme: add support for an id column
  ids = rownames(eset)
  names(ids) <- rownames(eset)
  cols = 1:ncol(eset)
  map = NA
  if (!is.na(mapping_column)) {
    #allgenes = rownames(eset)[which(Biobase::fData(eset)[,mapping_column] %in% allgenes)]

    # NOTE: assumes that the mapping column is 1 and everything else
    #       is valid data- which may not be the case
    cols = 1:ncol(eset)
    ids = SummarizedExperiment::rowData(eset)[,mapping_column]
    names(ids) <- rownames(eset)

  }
  allgenes_present <- names(ids[which(ids %in% intersect(allgenes,ids))])
  
  if (is.na(mapping_column)) {
    ids <- NA ##set this back so next command can use rownames
  }
  cols1 = colnames(eset)[which(colnames(eset) %in% set1)]
  cols2 = colnames(eset)[which(colnames(eset) %in% set2)]
  
  abcols1 = SummarizedExperiment::assay(eset,assay_name)[allgenes_present,cols1]
  abcols2 = SummarizedExperiment::assay(eset,assay_name)[allgenes_present,cols2]
  
  abcols1 = abcols1[which(rowSums(!is.na(abcols1)) > 0),]
  abcols2 = abcols2[which(rowSums(!is.na(abcols2)) > 0),]
  
  allgenes_cor1 = cor(t(abcols1), use = "p")
  allgenes_cor2 = cor(t(abcols2), use = "p")

  return(difference_enrichment_in_relationships(geneset, allgenes_cor1, 
                                                allgenes_cor2, 
                                                idmap = ids, tag = tag, mode = mode))
}
