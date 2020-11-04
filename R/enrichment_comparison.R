#' enrichment_comparison
#'
#' # access through leapr wrapper


enrichment_comparison <- function(geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                  fdr=0, matchset=NULL, longform=F, sample_comparison,
                                  min_p_threshold=NULL, tag=NA, sample_n=NULL){
  
  result = enrichment_in_abundance(geneset=geneset, abundance=abundance, mapping_column=mapping_column, abundance_column=abundance_column, fdr=fdr, matchset=matchset, longform=longform, sample_comparison=sample_comparison,min_p_threshold=min_p_threshold, tag=tag, sample_n=sample_n)
  
  return(result)
  
}