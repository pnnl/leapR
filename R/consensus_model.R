#' consensus_model
#'
#' consensus_model function description is...
#'
#' @param model_list is...
#' @param return.all default is TRUE
#' @param filter_count default is NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

consensus_model = function(model_list, return.all=T, filter_count=NULL) {
  # similar to the consensus_coefficients above but works with a
  #         list of models and outputs helpful things :)

  varsrle = rle(sort(unlist(sapply(1:length(model_list), function (i) names(model_list[[i]])))))

  results = rep(list(NA), length(varsrle$values))
  names(results) = varsrle$values

  for (i in 1:length(model_list)) {
    for (j in names(model_list[[i]])) {
      val = model_list[[i]][j]
      results[[j]] = c(results[[j]], val)
    }
  }

  meancoefficients = sapply(names(results), function (n) mean(results[[n]], na.rm=T))
  stddev = sapply(names(results), function (n) sd(results[[n]], na.rm=T))
  counts = sapply(names(results), function (n) sum(!is.na(results[[n]])))
  return(list(mean=meancoefficients, sd=stddev, counts=counts))
}
