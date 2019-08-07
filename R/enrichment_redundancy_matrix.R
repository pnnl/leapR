#' enrichment_redundancy_matrix
#'
#' enrichment_redundancy_matrix function description is...
#'
#' @param geneset is...
#' @param dataframe defaults to NA
#' @param enrichment_results defaults to NA
#' @param significance_threshold defaults to NA
#' @param pathway_list defaults to NA
#' @param method default is 'jaccard'
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

enrichment_redundancy_matrix <- function(geneset, dataframe=NA, enrichment_results=NA, significance_threshold=NA,
                                         pathway_list=NA, method="jaccard") {
  # takes an enrichment table and the data frame that created it- with the rownames corresponding to gene names
  #    and calculates a pairwise overlap matrix- that is, for each pair of pathways what is the overlap
  #    between the genes represented in those pathways?

  if (!is.na(significance_threshold)) enrichment_results = enrichment_results[which(enrichment_results[,7]<significance_threshold),]
  if (!is.na(pathway_list)) {
    enrichment_results = matrix(nrow=length(pathway_list), ncol=1)
    rownames(enrichment_results) = pathway_list
  }

  result_matrix = matrix(nrow=nrow(enrichment_results), ncol=nrow(enrichment_results))
  rownames(result_matrix) = rownames(enrichment_results)
  colnames(result_matrix) = rownames(enrichment_results)

  for (i in rownames(enrichment_results)) {
    gene_list_i = get_pathway_information(geneset, i)$geneset
    if (!is.na(dataframe)) gene_list_i = gene_list_i[which(gene_list_i %in% rownames(dataframe))]

    for (j in rownames(enrichment_results)) {
      if (j==i) overlap = 1.0
      else {
        gene_list_j = get_pathway_information(geneset, j)$geneset
        if (!is.na(dataframe)) gene_list_j = gene_list_j[which(gene_list_j %in% rownames(dataframe))]

        # for first pass calculate jaccard index
        total = length(unique(c(gene_list_i, gene_list_j)))
        over = sum(gene_list_i %in% gene_list_j)
        if (method == "jaccard") overlap = over/total
        else overlap = over/length(gene_list_j)
      }
      result_matrix[i,j] = overlap
    }
  }
  return(result_matrix)
}
