#'print.model
#'
#'
#'

print.model <- function(fit, level=NULL, pvalues=NULL) {
  if (is.null(level)) {
    level = which.min(pvalues)
  }
  model = fit$b.corrector[[level]]
  names(model) = line.to.gene.symbol(names(model), fit$gene.symbols)

  return(model)
}
