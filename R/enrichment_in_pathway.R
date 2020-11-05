#' enrichment_in_pathway
#'
#' Enrichment in pathway is a wrapper function for enrichment_in_abundance where 'sample_comparison' argument is set to NULL
#'
#' @noRd

enrichment_in_pathway <- function(geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                    fdr=0, matchset=NULL, longform=F, min_p_threshold=NULL, tag=NA, sample_n=NULL){
  
  result = enrichment_in_abundance(geneset=geneset, abundance=abundance, mapping_column=mapping_column, abundance_column=abundance_column, fdr=fdr, matchset=matchset, longform=longform, sample_comparison=NULL,min_p_threshold=min_p_threshold, tag=tag, sample_n=sample_n)
  
  return(result)
  
}