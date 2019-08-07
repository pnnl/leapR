#' ks_enrichment_vector
#'
#' ks_enrichment_vector function description is...
#'
#' @param groupgeneids is...
#' @param backgroundgenes is...
#' @param map_column defaults to 0
#' @param sort_col defaults to 1
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'
#'

ks_enrichment_vector <- function(groupgeneids, backgroundgenes, map_column=0, sort_col=1) {
  # returns a vector that is length(backgroundgenes) long of the
  #     enrichment of the groupgenes at each step in the ordered
  #     list of background genes
  if (map_column == 0) backgroundgenes = rownames(backgroundgenes[order(backgroundgenes[,sort_col]),])
  else backgroundgenes = backgroundgenes[order(backgroundgenes[,sort_col]),map_column]

  overallen = sum(groupgeneids %in% backgroundgenes)/length(backgroundgenes)

  result = c()
  for (i in 1:length(backgroundgenes)) {
    x = sum(groupgeneids %in% backgroundgenes[1:i])

    thisen = x/i

    result = c(result, thisen/overallen)
  }
  return(result)
}
