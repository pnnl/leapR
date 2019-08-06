#'partition_sets_ksenrich
#'
#'
#'

partition_sets_ksenrich <- function(partition_sets, geneset, abundance, mapping_column=NULL) {
  results = list()

  for (i in 1:length(partition_sets$names)) {
    thisname = partition_sets$names[i]
    thissize = partition_sets$size[i]
    description = partition_sets$desc[i]
    grouplist = partition_sets$matrix[i,1:thissize]
    cat(thisname, " : ", description, "\n")

    result = sapply(colnames(abundance), function (c)
      enrichment_in_groups(geneset, background=abundance[grouplist[which(grouplist %in% rownames(abundance))],c],
                           method="ks", mapping_column=mapping_column)[,"Signed_AdjP"])


    results = c(results, list(result))
  }
  names(results) = partition_sets$names
  return(results)
}
