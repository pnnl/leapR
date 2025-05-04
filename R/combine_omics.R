#' combine_omics
#'
#' Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.
#'
#' @param proteomics is a matrix of protein abundance values where rownames are ids and columns are conditions
#' @param proteomics is a matrix of transcript values where rownames are ids and columns are conditions
#' @param methylation is a matrix of methylation values where rownames are ids and columns are conditions
#' @param CNV is a matrix of copy number variant(CNV) values where rownames are ids and columns are conditions
#' @param phospho is a dataframe of phosphorylation data 
#' @param proteomics_tag is a text prefix to be added to protein ids
#' @param transcriptomics_tag is a text prefix to be added to transcript ids
#' @param methylation_tag is a text prefix to be added to methylation ids
#' @param cnv_tag is a text prefix to be added to cnv ids
#' @param phospho_tag is a text prefix to be added to phospho ids
#' @param id_column is an optional column number for identifiers for the phospho data
#' 
#' @details This combines matrices of different omics types together and adds prefix tags to the ids.
#'
#' @examples
#' dontrun{
#'         library(leapr)
#'
#'         # read in the example protein data
#'         data("protdata")
#'         # read in the example transcriptomics data
#'         data("transdata")
#'         
#'         # merge the two datasets by rows and add prefix tags for different omics types
#'         multi_omics = combine_omics(proteomics=protdata, transcriptomics=transdata)
#'
#' }
#'
#' @export


combine_omics = function(proteomics=NA, transcriptomics=NA, methylation=NA, cnv=NA, phospho=NA, proteomics_tag="prot_", 
                         transcriptomics_tag="txn_", methylation_tag="meth_", cnv_tag="cnv_", phospho_tag="phospho_", 
                         id_column=NA) {
  
  # find the common subset of colnames
  common_conditions = NA
  for (this in list(proteomics, transcriptomics, methylation, cnv, phospho)) {
    if (!all(is.na(this))) {
      if (all(is.na(common_conditions))) common_conditions = colnames(this)
      that = colnames(this)
      common_conditions = common_conditions[which(common_conditions %in% that)]
    }
  }
  if (length(common_conditions)==0) stop("No common conditions found")
  
  result = NA
  for (i in 1:5) {
    this = list(proteomics, transcriptomics, methylation, cnv, phospho)[[i]]
    tag = c(proteomics_tag, transcriptomics_tag, methylation_tag, cnv_tag, phospho_tag)[i]
    
    if (!all(is.na(this))) {
      this = this[,common_conditions]
      if (!is.na(id_column)) {
        # we need to add an id_column or use one that's here
        if (tag == phospho_tag) {
          # add the idcolumn from the input phospho data
          if (!all(is.na(phospho))) this = cbind(phospho[,id_column], this)
        }
        else {
          # add an idcolumn that is the rownames
          this = cbind(rownames(this), this)
        }
        colnames(this)[1] = "id"
      }
      
      # tag all the ids appropriately
      rownames(this) = sapply(rownames(this), function (n) paste(tag, n, sep=""))
      
      # we will also add tags to the idcolumn if necessary
      if (!is.na(id_column)) {
        this[,1] = sapply(this[,1], function (n) paste(tag, n, sep=""))
      }
      if (all(is.na(result))) result = this
      else result = rbind(result, this)
          }
  }
  return(result)
}
