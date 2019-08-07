#' print.model
#'
#' print.model function description is...
#'
#' @param fit is...
#' @param level defaults to NULL
#' @param pvalues defaults to NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

print.model <- function(fit, level=NULL, pvalues=NULL) {
  if (is.null(level)) {
    level = which.min(pvalues)
  }
  model = fit$b.corrector[[level]]
  names(model) = line.to.gene.symbol(names(model), fit$gene.symbols)

  return(model)
}
