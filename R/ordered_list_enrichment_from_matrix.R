#' ordered_list_enrichment_from_matrix
#'
#' ordered_list_enrichment_from_matrix function description is...
#'
#' @param geneset is...
#' @param datamatrix is...
#' @param starting_col defaults to NULL
#' @param ending_col defaults to NULL
#' @param idcolumn defaults to NULL
#' @param threshold defaults to 0.02
#' @param background_list defaults to NULL
#' @param enrichment_threshold defaults to 0.05
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

ordered_list_enrichment_from_matrix <- function(geneset, datamatrix, starting_col=NULL, ending_col=NULL, idcolumn=NULL,
                                                threshold=0.02, background_list=NULL, enrichment_threshold=0.05) {

  # currently the thresholding is setup to handle fdrs/pvalues where lower is better
  #            but this should be extended
  if (is.null(starting_col)) starting_col = 1
  if (is.null(ending_col)) ending_col = ncol(datamatrix)
  if (is.null(background_list)) {
    if (is.null(idcolumn)) background_list = rownames(datamatrix)
    else background_list = datamatrix[,idcolumn]
  }

  thesenames = rownames(datamatrix)
  if (!is.null(idcolumn)) thesenames = datamatrix[,idcolumn]

  results = list()
  idlists = list()
  plists = list()
  for (c in colnames(datamatrix)[starting_col:ending_col]) {
    exlist = thesenames[which(datamatrix[,c]<threshold)]
    plist = rownames(datamatrix)[which(datamatrix[,c]<threshold)]
    enriched = enrichment_in_groups(geneset, exlist, background_list)
    results = c(results, list(enriched[which(enriched[,"Adjusted_pvalue"]<enrichment_threshold),]))
    idlists = c(idlists, list(exlist))
    plists = c(plists, list(plist))
  }
  names(results) = colnames(datamatrix)[starting_col:ending_col]
  names(idlists) = colnames(datamatrix)[starting_col:ending_col]
  names(plists) = colnames(datamatrix)[starting_col:ending_col]
  return(list(enrichment=results, idlists=idlists, names=plists))
}
