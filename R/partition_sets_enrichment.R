#' partition_sets_enrichment
#'
#' partiton_sets_enrichment function description is...
#'
#' @param partition_sets is...
#' @param geneset is...
#' @param abundance is...
#' @param mapping_column defaults to NULL
#' @param abundance_column defaults to NULL
#' @param sample_comparison defaults to NULL
#' @param background_comparison defaults to NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
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
