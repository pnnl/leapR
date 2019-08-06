#'gene_set_overlap
#'
#'
#'

gene_set_overlap <- function(inmatrix, colsets=NA, filter=1.5) {
  # this will simply examine the overlap in successively more sets
  #      but will ignore which specific sets are overlapping. That is
  # 	   it will identify genes that are DE in ANY N conditions.
  #
  #  inmatrix : matrix of values to examine
  #  colsets  : a list of vectors describing groups to look for DE in
  #  filter   : what DE do we consider to be DE?
  #
  # To create a colset list, e.g.:
  #     conds = list(colnames(ratios)[45:59], colnames(ratios)[61:66], colnames(ratios)[77:86])
  outmatrix = matrix(nrow=nrow(inmatrix), ncol=length(colsets))
  rownames(outmatrix) = rownames(inmatrix)
  ids = c()
  for (gene in rownames(inmatrix)) {
    for (i in 1:length(colsets)) {
      if (sum(abs(inmatrix[gene,colsets[[i]]])>filter)) {
        outmatrix[gene,i] = TRUE
      } else {
        outmatrix[gene,i] = FALSE
      }
    }
    ids = c(ids, sum(outmatrix[gene,])>0)
  }
  # filter the matrix so it only has shared DE genes
  outmatrix = outmatrix[ids,]
  return(outmatrix)
}
