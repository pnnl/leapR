#' combine_omics
#'
#' Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.
#'
#' @import Biobase
#' @param proteomics is a matrix of protein abundance values where rownames are ids and columns are conditions
#' @param transcriptomics is a matrix of gene expression values  where rownames are ids and columns are conditions
#' @param proteomics is a matrix of transcript values where rownames are ids and columns are conditions
#' @param methylation is a matrix of methylation values where rownames are ids and columns are conditions
#' @param cnv is a matrix of copy number variant(CNV) values where rownames are ids and columns are conditions
#' @param phospho is a dataframe of phosphorylation data 
#' @param proteomics_tag is a text prefix to be added to protein ids
#' @param transcriptomics_tag is a text prefix to be added to transcript ids
#' @param methylation_tag is a text prefix to be added to methylation ids
#' @param cnv_tag is a text prefix to be added to cnv ids
#' @param phospho_tag is a text prefix to be added to phospho ids
#' @param id_column is an optional column number for identifiers for the phospho data
#' 
#' @return combined omics table 
#' @details This combines matrices of different omics types together and adds prefix tags to the ids.
#'
#' @examples
#'         library(leapR)
#'
#'         # read in the example protein data
#'         pdata <- download.file('https://figshare.com/ndownloader/files/55158767',method='libcurl',destfile='protdata')
#'         protdata<-read.csv("protdata",check.names=FALSE,row.names=1)|>
#'             as.matrix()
#'         file.remove("protdata")
#'
#'           
#'         tdata <- download.file("https://figshare.com/ndownloader/files/55158764",method='libcurl',destfile='transdata')
#'         transdata<-readr::read_csv('transdata')|>
#'             tibble::column_to_rownames('...1')|>
#'             as.matrix()
#'         file.remove("transdata")
#'
#'         # merge the two datasets by rows and add prefix tags for different omics types
#'         multi_omics = combine_omics(proteomics=protdata, transcriptomics=transdata)
#'
#'
#' @export


combine_omics <- function(omics_list, id_list = rep(NA,length(omics_list))){
  
  # find the common subset of colnames
  common_conditions = NA
  for (this in omics_list) {
    #if (!all(is.na(this))) {
      if (all(is.na(common_conditions))) common_conditions = colnames(this)
      that = colnames(this)
      common_conditions = common_conditions[which(common_conditions %in% that)]
    #}
  }
  
  if (length(common_conditions) == 0) stop("No common conditions found")
  
  result <- NA
  for (i in 1:length(omics_list)) {
    this = omics_list[[i]]
    id_column = id_list[[i]]
    data = this@annotation
    
    tag = paste0(substr(data,0,3),'_')
    #tag = c(proteomics_tag, transcriptomics_tag, methylation_tag, cnv_tag, phospho_tag)[i]
    
    #if (!all(is.na(this))) {
    this = this[,common_conditions]
      
    if (!is.na(id_column)) {
        # we need to add an id_column or use one that's here
          # add the idcolumn from the input phospho data
          this = cbind(Biobase::featureData(this)[[id_column]], Biobase::exprs(this))
    }else {
        this = cbind(rownames(this),Biobase::exprs(this))
      }
        
        #else {
        # add an idcolumn that is the rownames
        #this = cbind(rownames(this), this)
        
     colnames(this)[1] = "id"
      
      # tag all the ids appropriately
      rownames(this) = sapply(rownames(this), function (n) paste(tag, n, sep=""))
      
      # we will also add tags to the idcolumn if necessary
    #  if (!is.na(id_column)) {
      this[,1] = sapply(this[,1], function (n) paste(tag, n, sep = ""))
    #  }
      
      #  print(dim(this))
      #    print(dim(result))
      if (all(is.na(result))) result = this
      else result = rbind(result, this)
      
  }
  
  nres <- apply(result[,common_conditions],2,as.numeric)
  rownames(nres) <- rownames(result)
  print(dim(result))
  if ('id' %in% colnames(result))
    ids = result[,'id']
  else
    ids = rownames(result)
  
  return(data.frame(id=ids,nres,check.names=FALSE))
}

old_combine_omics = function(proteomics=NA, transcriptomics=NA, methylation=NA, cnv=NA, phospho=NA, proteomics_tag="prot_", 
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
        if(!all(is.na(phospho))) {}
        if (tag == phospho_tag) {
          # add the idcolumn from the input phospho data
          this = cbind(phospho[,id_column], this)
        }else {
          this = cbind(rownames(this),this)
        }
      
      #else {
          # add an idcolumn that is the rownames
          #this = cbind(rownames(this), this)
      
      colnames(this)[1] = "id"
      }
      # tag all the ids appropriately
      rownames(this) = sapply(rownames(this), function (n) paste(tag, n, sep=""))
      
      # we will also add tags to the idcolumn if necessary
      if (!is.na(id_column)) {
        this[,1] = sapply(this[,1], function (n) paste(tag, n, sep=""))
      }
 
      #  print(dim(this))
    #    print(dim(result))
      if (all(is.na(result))) result = this
      else result = rbind(result, this)
    }
  }
  
  ##SG: added in second check to mak esure all values are numeric
  nres <- apply(result[,common_conditions],2,as.numeric)
  rownames(nres) <- rownames(result)
  if ('id' %in% colnames(result))
    ids = result$id
  else
    ids = rownames(result)
  return(data.frame(id=ids,nres,check.names=FALSE))
}
