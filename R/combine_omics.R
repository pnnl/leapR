#' combine_omics
#'
#' Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.
#'
#' @import SummarizedExperiment
#' @param omics_list Is a list of \code{SummarizedExperiment} each with one assay
#' @param id_list List of identifiers to use, in the same order as the omics_list elements. If an element
#' is `NA`, then rownames are used.
#' @return \code{SummarizedExperiment} with an additional assay called `combined`
#' @details This combines matrices of different omics types together and adds prefix tags to the ids.
#'
#' @examples
#'         library(leapR)
#'
#'
#'         pdata <- download.file('https://api.figshare.com/v2/file/download/56536217',method='libcurl',destfile='protData.rda')
#'         load('protData.rda')
#'         p <- file.remove("protData.rda")
#'         
#'         tdata <- download.file("https://api.figshare.com/v2/file/download/56536214",method='libcurl',destfile='transData.rda')
#'         load('transData.rda')
#'         p <- file.remove("transData.rda")
#'         
#'         phdata<-download.file('https://api.figshare.com/v2/file/download/56536211',method='libcurl',destfile = 'phosData.rda')
#'         #phosphodata<-read.csv("phdata",check.names=FALSE,row.names=1)
#'         load('phosData.rda')
#'         p <- file.remove('phosData.rda')# read in the example protein data
#'         
#'
#'         # merge the three datasets by rows and add prefix tags for different omics types
#'         multi_omics = combine_omics(list(pset, tset, phset), list(NA,NA,'hgnc_id'))
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
  allfeat <- NA
  for (i in 1:length(omics_list)) {
    this = omics_list[[i]]
    id_column = id_list[[i]]
    
    ##first lets pull some metadata out
    feat =  names(SummarizedExperiment::assays(this))[1]
    tag = paste0(substr(feat,0,3),'_')
    
    #now let's get the common conditions
    this = this[,common_conditions]
      
    if (!is.na(id_column)) {
        # we need to add an id_column or use one that's here
          # add the idcolumn from the input phospho data
      ids <- SummarizedExperiment::rowData(this)[,id_column]
    }else {
#        this = cbind(rownames(this),Biobase::exprs(this))
      ids <- rownames(this)
    }
        
        #else {
        # add an idcolumn that is the rownames
        #this = cbind(rownames(this), this)
        
     #colnames(this)[1] = "id"
      idval <- data.frame(id = ids)
      
      # tag all the ids appropriately
      expr <- SummarizedExperiment::assay(this,feat)
      rownames(expr) = sapply(rownames(expr), function(n) paste0(tag, n, sep = "_"))
      
      #  print(dim(this))
      #    print(dim(result))
      if (all(is.na(result))) result = expr
      else result = rbind(result, expr)
      
      if (all(is.na(allfeat))) allfeat = idval
      else allfeat = rbind(allfeat,idval)
      
  }
  
  allfeat <- data.frame(allfeat)
  rownames(allfeat) <- rownames(result)

  return(SummarizedExperiment::SummarizedExperiment(assays = as(list(combined=result),'SimpleList'), rowData = allfeat))
}
# 
# old_combine_omics = function(proteomics=NA, transcriptomics=NA, methylation=NA, cnv=NA, phospho=NA, proteomics_tag="prot_", 
#                          transcriptomics_tag="txn_", methylation_tag="meth_", cnv_tag="cnv_", phospho_tag="phospho_", 
#                          id_column=NA) {
#   
#   # find the common subset of colnames
#   common_conditions = NA
#   for (this in list(proteomics, transcriptomics, methylation, cnv, phospho)) {
#     if (!all(is.na(this))) {
#       if (all(is.na(common_conditions))) common_conditions = colnames(this)
#       that = colnames(this)
#       common_conditions = common_conditions[which(common_conditions %in% that)]
#     }
#   }
#   if (length(common_conditions)==0) stop("No common conditions found")
#   
#   result = NA
#   for (i in 1:5) {
#     this = list(proteomics, transcriptomics, methylation, cnv, phospho)[[i]]
#     tag = c(proteomics_tag, transcriptomics_tag, methylation_tag, cnv_tag, phospho_tag)[i]
#     
#     if (!all(is.na(this))) {
#       this = this[,common_conditions]
#       
#       if (!is.na(id_column)) {
#         # we need to add an id_column or use one that's here
#         if(!all(is.na(phospho))) {}
#         if (tag == phospho_tag) {
#           # add the idcolumn from the input phospho data
#           this = cbind(phospho[,id_column], this)
#         }else {
#           this = cbind(rownames(this),this)
#         }
#       
#       #else {
#           # add an idcolumn that is the rownames
#           #this = cbind(rownames(this), this)
#       
#       colnames(this)[1] = "id"
#       }
#       # tag all the ids appropriately
#       rownames(this) = sapply(rownames(this), function (n) paste(tag, n, sep=""))
#       
#       # we will also add tags to the idcolumn if necessary
#       if (!is.na(id_column)) {
#         this[,1] = sapply(this[,1], function(n) paste(tag, n, sep = ""))
#       }
#  
#       #  print(dim(this))
#     #    print(dim(result))
#       if (all(is.na(result))) result = this
#       else result = rbind(result, this)
#     }
#   }
#   
#   ##SG: added in second check to mak esure all values are numeric
#   nres <- apply(result[,common_conditions],2,as.numeric)
#   rownames(nres) <- rownames(result)
#   if ('id' %in% colnames(result))
#     ids = result$id
#   else
#     ids = rownames(result)
#   return(data.frame(id = ids,nres,check.names=FALSE))
# }
