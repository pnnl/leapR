#' overlap_gene_sets
#'
#' overlap_gene_sets function description is...
#'
#' @param geneset is...
#' @param set1 is...
#' @param set2 is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

overlap_gene_sets <- function(geneset, set1, set2) {
  a = get_pathway_information(geneset, set1)$geneset
  b = get_pathway_information(geneset, set2)$geneset

  overlap = a[a %in% b]
  uniquea = a[!a %in% b]
  uniqueb = b[!b %in% a]

  return(list(overlap=overlap, set1_unique=uniquea, set2_unique=uniqueb))

}
