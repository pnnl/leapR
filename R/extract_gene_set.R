#' extract_gene_set
#'
#' extract_gene_set function description is...
#'
#' @param geneset is...
#' @param setlist is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

extract_gene_set <- function(geneset, setlist) {
  x = which(geneset$names %in% setlist)
  names = geneset$names[x]
  notfound = length(setlist[which(!setlist %in% geneset$names)])
  cat(notfound, "sets not found in geneset\n")

  desc = geneset$desc[x]
  size = geneset$sizes[x]
  mat = geneset$matrix[x,]

  return(list(names=names, desc=desc, sizes=size, matrix=mat))
}
