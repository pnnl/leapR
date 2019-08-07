#' cv_selected_variable_analysis
#'
#' cv_selected_varibale_analysis function description is...
#'
#' @param fitcv is...
#' @param level is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

cv_selected_variable_analysis = function(fitcv, level) {
  varsrle = rle(sort(unlist(sapply(1:length(fitcv), function (i) names(fitcv[[i]][[1]]$b.corrector[[level]])))))

  results = matrix(0, nrow=length(fitcv[[1]][[2]]), ncol=length(varsrle$values))
  colnames(results) = varsrle$values
  rownames(results) = names(fitcv[[1]][[2]])
  for (i in 1:length(fitcv)) {
    for (j in names(fitcv[[i]][[1]]$b.corrector[[level]])) {
      results[,j] = results[,j] + as.integer(fitcv[[i]][[2]])
    }
  }
  return(results)
}
