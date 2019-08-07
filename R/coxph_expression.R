#' coxph_expression
#'
#' coxph_expression function description is...
#'
#' @param survivalmatrix is...
#' @param statusvar is...
#' @param timevar is...
#' @param startcolumn is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

coxph_expression <- function(survivalmatrix, statusvar, timevar, startcolumn) {
  # function to build a coxph model from each expression pattern
  # in a datamatrix.

  results = matrix(nrow=ncol(survivalmatrix[,startcolumn:ncol(survivalmatrix)]), ncol=2)
  rownames(results) = colnames(survivalmatrix[,startcolumn:ncol(survivalmatrix)])

  for (col in colnames(survivalmatrix)[startcolumn:ncol(survivalmatrix)]) {
    thisform = paste("Surv(",timevar,",",statusvar,") ~", col, sep="")
    model = try(coxph(as.formula(thisform), data=survivalmatrix), silent=F)
    if (class(model) != "try-error") {
      results[col,1] = summary(model)$waldtest[3]
      results[col,2] = summary(model)$logtest[3]
    }
  }
  return(results)
}
