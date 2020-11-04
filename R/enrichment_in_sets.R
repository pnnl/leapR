#' enrichment_in_sets
#'
#' Enrichment in sets is a wrapper function for enrichment_in_groups where 'method' argument is set to "fishers" and "abundance_column" is set to NA
#'# access through leapr wrapper
#'
enrichment_in_sets <- function(geneset, targets=NULL, background=NULL, minsize=5,
                                 mapping_column=NULL, randomize=F){
  
  result = enrichment_in_groups(geneset=geneset, targets=targets, background=background, method="fishers", minsize=5, mapping_column=mapping_column, abundance_column=NA, randomize=randomize)
  
  return(result)
  
}