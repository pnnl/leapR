#' calcTTest
#'
#' calculates a t-test for two distributions of data on a per-gene basis
#' append results to ExpressionSet with two extra columns: `pvalue` and
#' `difference` for each feature

#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @param eset SummarizedExperiment
#' @param assay_name name of assay
#' @param group1 List of samples comprising group 1
#' @param group2 List of samples comprising group 2
#' @return An Expression set with two columns added to the featureData
#' slot: pvalue, and estimate
#' @export
#'
#' @examples
#'
#'         library(leapR)
#'         url <- "https://api.figshare.com/v2/file/download/56536214"
#'         tdata <- download.file(url,method='libcurl',destfile='transData.rda')
#'         load('transData.rda')
#'         p <- file.remove("transData.rda")
#'
#'         # read in the pathways
#'         data("ncipid")
#'
#'         # read in the patient groups
#'         data("shortlist")
#'         data("longlist")
#'         calcTTest(tset, 'transcriptomics', shortlist, longlist)


calcTTest <- function(eset, assay_name, group1, group2){
  stopifnot(is(eset,'SummarizedExperiment'),
            length(group1) > 0,
            length(group2) > 0)

  exprdata <- SummarizedExperiment::assay(eset, assay_name)

  res <- t(vapply(rownames(exprdata), function(r) {
    if (length(which(!is.na(exprdata[r,group1, drop = FALSE]))) > 2 &&
        length(which(!is.na(exprdata[r,group2, drop = FALSE]))) > 1) {
      res <- t.test(exprdata[r,group1, drop = FALSE],
                    exprdata[r,group2, drop = FALSE])

      return(c(pvalue = res$p.value,
               difference  = res$estimate[[1]] - res$estimate[[2]]))
    }
    else{
      return(c(pvalue = NA,difference = NA))
    }
  }, numeric(2)))

  f <- SummarizedExperiment::rowData(eset)
  f <- cbind(f, res)

  SummarizedExperiment::rowData(eset) <- f
  return(eset)
}
