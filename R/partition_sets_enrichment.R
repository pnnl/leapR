#'partition_sets_enrichment
#'
#'
#'
#'

partition_sets_enrichment <- function(partition_sets, geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                      sample_comparison=NULL, background_comparison=NULL) {
  results = list()

  for (i in 1:length(partition_sets$names)) {
    thisname = partition_sets$names[i]
    thissize = partition_sets$size[i]
    description = partition_sets$desc[i]
    grouplist = partition_sets$matrix[i,1:thissize]
    cat(thisname, " : ", description, "\n")

    result = enrichment_in_abundance(geneset, abundance[grouplist[which(grouplist %in% rownames(abundance))],],
                                     mapping_column=mapping_column, abundance_column=abundance_column,
                                     sample_comparison=sample_comparison, background_comparison=background_comparison)

    results = c(results, list(result))
  }
  names(results) = partition_sets$names
  return(results)
}
