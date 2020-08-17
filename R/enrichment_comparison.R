#' enrichment_comparison
#'
#' Enrichment comparison is a wrapper function for enrichment_in_abundance where 'sample_comparison' argument is a required input.
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param Is a \emph{mxn} matrix of abundance data (protein, gene, etc.), with \emph{m} gene names (rows) and \emph{n} sample/condition (columns).
#' @param mapping_column Is an optional character string, a column name of \code{abundance}, for specifying gene mapping (e.g. for phosphoproteomics data). Defaults to NULL
#' @param abundance_column Is a character vector composed of column names from \code{abundance}, that identifies which conditions should be tested.
#' @param fdr A numerical value which specifies how many times to randomly sample genes, defaults to 0.
#' @param matchset defaults to NULL
#' @param sample_comparison Is a character vector composed of column names from \code{abundance}, that identifies which conditions should be compared.
#' @param min_p_threshold Is a numeric value, a lower p-value threshold for returning results, defaults to NULL.
#' @param a prefix tag (optional) for specifying source of data default is NA
#' @param sample_n If set specifies the number of randomizations to perform to calculate an FDR [warning: time intensive], default is NA
#' 
#' @export


enrichment_comparison <- function(geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                  fdr=0, matchset=NULL, longform=F, sample_comparison,
                                  min_p_threshold=NULL, tag=NA, sample_n=NULL){
  
  result = enrichment_in_abundance(geneset=geneset, abundance=abundance, mapping_column=mapping_column, abundance_column=abundance_column, fdr=fdr, matchset=matchset, longform=longform, sample_comparison=sample_comparison,min_p_threshold=min_p_threshold, tag=tag, sample_n=sample_n)
  
  return(result)
  
}