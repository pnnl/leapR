#' pvalue_from_random_gene_set_overlap_db
#'
#' pvalue_from_random_gene_set_overlap_db function description is...
#'
#' @param inmatrix_1 is...
#' @param inmatrix_2 is...
#' @param names_1 is...
#' @param names_2 is...
#' @param rep defaults to 1000
#' @param longform defaults to FALSE
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

pvalue_from_random_gene_set_overlap_db <- function(inmatrix_1, inmatrix_2, names_1, names_2, rep=1000, longform=F) {
  # calculates a pvalue from random iterations for overlap between two sets
  x = sum(names_1 %in% names_2)

  pvect = random_gene_set_overlap_db(inmatrix_1, inmatrix_2, length(names_1), length(names_2), rep=rep)

  pvalue = sum(pvect>=x)/length(pvect)

  if (longform) return(list(x=x, pvalue=pvalue, pvect=pvect, n1=length(names_1), n2=length(names_2)))

  return(pvalue)
}
