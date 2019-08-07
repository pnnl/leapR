#' make_gene_set_matrix
#'
#' make_gene_set_matrix function description is...
#'
#' @param database is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

make_gene_set_matrix <- function(database) {
  # takes a gene set database structure and returns
  #      a matrix that relates genes (rows) to their
  #      pathways (columns). This is useful because
  #      there can be redundance across different pathways
  #      and this structure can be used to figure out
  #      commonalities.
  ids = unique(as.vector(database$matrix))
  ids = ids[which(!ids == "null")]

  result = matrix(0, nrow=length(ids), ncol=length(database$names), dimnames=list(ids, database$names))

  x = 1
  for (this in database$names) {
    result[,this] = as.numeric(!rownames(result) %in% database$matrix[x,])
    x = x + 1
  }

  return(result)
}
