#'coxph_subset
#'
#'
#'

coxph_subset <- function(survivalmatrix, statusvar="status", timevar="time", columns, rows, verbose=F) {
  # function to build a coxph model from an input data matrix given a range of columns and rows
  # the orientation of the matrix is columns=proteins, rows=examples
  ids = columns
  results = NA

  thisform = paste("Surv(",timevar,",",statusvar,") ~", paste(ids, collapse= "+"), sep="")
  if (verbose) cat(thisform, "\n")

  model = try(coxph(as.formula(thisform), data=survivalmatrix[rows,]), silent=T)

  if (class(model) != "try-error") {
    results = list(waldp=summary(model)$waldtest[3], logp=summary(model)$logtest[3], model=model)
  }
  return(results)
}
