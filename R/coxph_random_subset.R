#' coxph_random_subset
#'
#' coxph_random_subset function description...
#'
#' @param survivalmatrix is...
#' @param statusvar is...
#' @param timevar is...
#' @param startcolumn is...
#' @param size default is 5
#' @param ntimes default is 50
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

coxph_random_subset <- function(survivalmatrix, statusvar, timevar, startcolumn, size=5, ntimes=50) {
  # function to build a coxph model from each expression pattern
  # in a datamatrix.

  results = matrix(nrow=ntimes, ncol=2)
  ids = colnames(survivalmatrix)[startcolumn:ncol(survivalmatrix)]

  for (i in 1:ntimes) {
    randsub = sample(ids, size)
    thisform = paste("Surv(",timevar,",",statusvar,") ~", paste(randsub, collapse= "+"), sep="")
    model = try(coxph(as.formula(thisform), data=survivalmatrix))
    if (class(model) != "try-error") {
      results[i,1] = summary(model)$waldtest[3]
      results[i,2] = summary(model)$logtest[3]
    }
  }
  return(results)
}
