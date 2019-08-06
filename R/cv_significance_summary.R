#'cv_significance_summary
#'
#'
#'
#'

cv_significance_summary <- function(mintable, results_struct, pthresh=0.05, return.percentage=T) {
  # mintable is a matrix of pathways (rows) by cv runs (cols) with the min pvalue
  #         from the cv in each cell
  # results is the results structure from the LEAP run which includes which tumors
  #         were included in the training and testing sets
  # The output is a table of pathways by tumors that provides the counts of the times
  #         where the tumor was included in a cv testing set that gave the pathway
  #         as significant by the pthresh
  results = matrix(nrow=nrow(mintable), ncol=length(results_struct$fits[[1]][[2]]$training.set))
  rownames(results) = rownames(mintable)
  colnames(results) = rownames(results_struct$data[[1]][[1]])

  counts = results

  rmatrix = sapply(1:ncol(mintable), function (i) (!results_struct$fits[[i]][[2]]$training.set))

  for (path in rownames(mintable)) {
    sigvec = matrix(mintable[path,]<pthresh, ncol=ncol(mintable), nrow=ncol(results))
    if (return.percentage == T) results[path,] = rowSums(sigvec * rmatrix)/rowSums(rmatrix)
    else results[path,] = rowSums(sigvec * rmatrix)
    counts[path,] = rowSums(rmatrix)
  }
  return(list(results, counts))
}
