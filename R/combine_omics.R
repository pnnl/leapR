#' combine_omics
#'
#' Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
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
#'         multi_omics <- combine_omics(list(pset, tset, phset), list(NA,NA,'hgnc_id'))
#'
#'
#' @export


combine_omics <- function(omics_list, id_list = rep(NA,length(omics_list))){

  # find the common subset of colnames
  common_conditions <- NA
  for (this in omics_list) { ##cannot seem to get rid of for loops
      if (all(is.na(common_conditions))) common_conditions <- colnames(this)
      that <- colnames(this)
      common_conditions <- common_conditions[which(common_conditions %in% that)]

  }

  if (length(common_conditions) == 0) stop("No common conditions found")

  result <- NA
  allfeat <- NA
  for (i in seq_along(1:length(omics_list))) { ##cannot seem to get rid of for loops for this test
    this <- omics_list[[i]]
    id_column <- id_list[[i]]

    ##first lets pull some metadata out
    feat <-  names(SummarizedExperiment::assays(this))[1]
    tag <- paste0(substr(feat,0,3),'_')

    #now let's get the common conditions
    this <- this[,common_conditions, drop = FALSE]

    if (!is.na(id_column)) {
        # we need to add an id_column or use one that's here
          # add the idcolumn from the input phospho data
      ids <- SummarizedExperiment::rowData(this)[,id_column, drop = TRUE]
    }else {
      ids <- rownames(this)
    }

      idval <- data.frame(id = ids)

      # tag all the ids appropriately
      expr <- SummarizedExperiment::assay(this,feat)
      rownames(expr) <- vapply(rownames(expr), function(n) paste0(tag, n, sep = "_"), "")

      if (all(is.na(result))) result <- expr
      else result <- rbind(result, expr)

      if (all(is.na(allfeat))) allfeat <- idval
      else allfeat <- rbind(allfeat,idval)

  }

  allfeat <- data.frame(allfeat)
  rownames(allfeat) <- rownames(result)

  return(SummarizedExperiment::SummarizedExperiment(assays = as(list(combined = result),'SimpleList'), rowData = allfeat))
}
