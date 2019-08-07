#' pathway_model_predictor_matrix
#'
#' pathway_model_predictor_matrix function description is...
#'
#' @param pathfit is...
#' @param validation_data defaults to NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

pathway_model_predictor_matrix = function(pathfit, validation_data=NULL) {
  # this builds a matrix of predictions of pathway models from the data
  np = length(pathfit$fits)
  result = matrix(ncol=nrow(pathfit$data[[1]][[1]]), nrow=np)
  colnames(result) = rownames(pathfit$data[[1]][[1]])
  rownames(result) = names(pathfit$fits)

  pvalues = rep(1, length(names(pathfit$fits)))
  names(pvalues) = names(pathfit$fits)

  for (i in 1:np) {
    thisfit = pathfit$fits[[i]]

    pathname = names(pathfit$bestmodels[i])
    bestx = which.min(pathfit$pvalues[i,])

    if (any(sapply(1:ncol(pathfit$pvalues), function (c) !is.na(min(pathfit$pvalues[i,c]))))) {
      #cat(pathname, " - ", bestx, "\n")
      if (!is.null(validation_data)) {
        sig = build_predictor(thisfit, validation_data[[1]], bestx,
                              genesymbols_from=thisfit$gene.symbols, genesymbols_to=validation_data[[2]])
        model_pvalue = plot_model_kaplan_meier(thisfit,
                                               validation_data[[1]], 0,
                                               prediction=sig, plotthis=F)
        pvalues[pathname] = model_pvalue
      }
      sig = build_predictor(thisfit, pathfit$data[[i]][[1]], bestx)
      result[pathname,] = sig
    }
  }
  if (!is.null(validation_data)) {
    return(list(matrix=result, pvalues=pvalues))
  }
  return(result)
}
