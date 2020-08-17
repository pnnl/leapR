#' enrichment_in_sets
#'
#' Enrichment in sets is a wrapper function for enrichment_in_groups where 'method' argument is set to "fishers" and "abundance_column" is set to NA
#'
#' @param genesets is a GeneSet object for pathway annotation
#' @param targets defaults to NULL
#' @param background Is a \emph{mxn} matrix of gene expression data, with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param minsize Is the minimum size of a gene set that will be considered, defaults to 5
#' @param mapping_column Is an optional character string, a column name of \code{abundance}, for specifying gene mapping (e.g. for phosphoproteomics data). Defaults to NULL
#' @param randomize is a logical, defaults to FALSE
#' 
#' @export

enrichment_in_sets <- function(geneset, targets=NULL, background=NULL, minsize=5,
                                 mapping_column=NULL, randomize=F){
  
  result = enrichment_in_groups(geneset=geneset, targets=targets, background=background, method="fishers", minsize=5, mapping_column=mapping_column, abundance_column=NA, randomize=randomize)
  
  return(result)
  
}