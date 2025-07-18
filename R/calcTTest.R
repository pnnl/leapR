#' calcTTest
#'
#' calculates a t-test for two distributions of data on a per-gene basis
#' append results to ExpressionSet with two extra columsn: `pvalue` and `difference` for each feature

#' @import Biobase
#' @param eset ExpressionSet
#' @param group1 List of samples comprising group 1
#' @param group2 List of samples comprising group 2
#' @return An Expression set with two colunns added to the featureData slot: pvalue, and estimate
#' @export
#' 
#' @examples
#' 
#'         library(leapR)
#'         tdata <- download.file("https://api.figshare.com/v2/file/download/55781153",method='libcurl',destfile='transData.rda')
#'         load('transData.rda')
#'         p <- file.remove("transData.rda")
#'
#'         # read in the pathways
#'         data("ncipid")
#'
#'         # read in the patient groups
#'         data("shortlist")
#'         data("longlist")
#'         calcTTest(tset, shortlist, longlist) 


calcTTest <- function(eset, group1, group2){
  exprdata <- Biobase::exprs(eset)
  
  res <- t(sapply(rownames(exprdata), function(r) {
    if (length(which(!is.na(exprdata[r,shortlist]))) > 2 && length(which(!is.na(exprdata[r,longlist]))) > 1) {
      res <- t.test(exprdata[r,shortlist], exprdata[r,longlist])
      #               error=return(c(NA,NA))) 
      #  if (is(res,'try-error')) return(c(NA, NA));
      return(c(pvalue = res$p.value, difference  = res$estimate[[1]] - res$estimate[[2]]))
    }
    else{
      return(c(pvalue = NA,difference = NA))
    }
  }))

  f <- Biobase::fData(eset)
  f <- cbind(f, res)

  Biobase::fData(eset) <- f
  return(eset)
}