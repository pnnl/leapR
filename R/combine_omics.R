#' combine_omics
#'
#' Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.
#'
#' @param proteomics is a matrix of protein abundance values where rownames are ids and columns are conditions
#' @param proteomics is a matrix of transcript values where rownames are ids and columns are conditions
#' @param methylation is a matrix of methylation values where rownames are ids and columns are conditions
#' @param CNV is a matrix of copy number variant(CNV) values where rownames are ids and columns are conditions
#' @param proteomics_tag is a text prefix to be added to protein ids
#' @param transcriptomics_tag is a text prefix to be added to transcript ids
#' @param methylation_tag is a text prefix to be added to methylation ids
#' @param cnv_tag is a text prefix to be added to cnv ids
#' 
#' @details This combines matrices of different omics types together and adds prefix tags to the ids.
#'
#' @examples
#' dontrun{
#'   multi_omics = comibine_omics(protmatrix, transcriptmatrix)
#'
#' }
#'
#' @export


combine_omics = function(proteomics=NA, transcriptomics=NA, methylation=NA, cnv=NA, proteomics_tag="prot_", 
                         transcriptomics_tag="txn_", methylation_tag="meth_", cnv_tag="cnv_") {
  
  # find the common subset of colnames
  common_conditions = NA
  for (this in list(proteomics, transcriptomics, methylation, cnv)) {
    if (!is.na(this)) {
      if (is.na(common_conditions)) common_conditions = colnames(this)
      that = colnames(this)
      common_conditions = common_conditions[which(common_conditions %in% that)]
    }
  }
  if (!common_conditions) stop("No common conditions found")
  
  result = NA
  for (i in 1:4) {
    this = list(proteomics, transcriptomics, methylation, cnv)[i]
    tag = list(proteomics_tag, transcriptomics_tag, methylation_tag, cnv_tag)[i]
    
    this = this[common_conditions]
    
    if (!is.na(this)) {
      # tag all the ids appropriately
      rownames(this) = sapply(rownames(this), function (n) (paste(tag, n, sep="")))
      if (is.na(result)) result = this
      else result = rbind(result, this)
    }
    
    return(result)
  }
}
