#' add_gene_set
#'
#' add_gene_set function Add a gene set to a GeneSet object.
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param name is the name of the gene set to be added
#' @param genelist is a list of the gene ids to be added
#' @param desc is a short description of the gene set 
#' 
#' @details This function adds a single new gene set to an existing GeneSet object. It returns
#'              the same (now modified) GeneSet object.
#'
#' @examples
#' dontrun{
#'   mygeneset = make_gene_set()
#'   mygeneset = add_gene_set(mygeneset, name="MyPathway", 
#'                            desc="A pathway describing an important thing",
#'                            genelist=c("geneA","geneB","geneX","geneY"))
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
