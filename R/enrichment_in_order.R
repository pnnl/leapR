#' enrichment_in_order
#'
#' Enrichment in order is a wrapper function for enrichment_in_groups where 'method' argument is set to "ks" and "targets" is set to NA
#'
#' 

enrichment_in_order <- function(geneset, background=NULL, minsize=5,
                               mapping_column=NULL, abundance_column=NULL, randomize=F){
  
  result = enrichment_in_groups(geneset=geneset, targets=NA, background=background, method="ks", minsize=5, mapping_column=mapping_column, abundance_column=abundance_column, randomize=randomize)
  
  return(result)
  
}