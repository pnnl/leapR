#' add_gene_set
#'
#' add_gene_set function description...
#'
#' @param geneset is...
#' @param name default is an NA
#' @param genelist default is an NA
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export


add_gene_set <- function(geneset, name=NA, desc=NA, genelist=NA) {
  geneset$names = c(geneset$names, name)
  geneset$desc = c(geneset$desc, desc)
  geneset$size = c(geneset$size, length(genelist))
  if (ncol(geneset$matrix)<length(genelist)) {
    geneset$matrix = cbind(geneset$matrix, rep("null", nrow(geneset$matrix)*(length(genelist)-ncol(geneset$matrix))))
  }
  add = genelist
  if ((ncol(geneset$matrix)-length(genelist))>0) add = c(genelist, rep("null", (ncol(geneset$matrix)-length(genelist))))

  geneset$matrix = rbind(geneset$matrix, add)

  return (geneset)
}
