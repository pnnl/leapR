#'random_gene_set_overlap
#'
#'
#'

random_gene_set_overlap <- function(inmatrix, colsets=NA, filter=1.5, N=NA, rep=100) {
  # this just returns the number of genes DE in N or more conditions
  #   if N==NA, N defaults to all conditions

  # this returns a list of the number of genes that are DE in each condition
  de_list = c()
  for (cond in colsets) {
    x = sapply(rownames(inmatrix), function (r) sum(abs(inmatrix[r,cond[1]]))>filter)
    de_list = c(de_list, sum(x))
    #cat("Condition", str(cond), "\n")
  }

  #cat("here\n")
  #return(de_list)
  results = c()

  for (i in 1:rep) {
    # now randomly select genes
    de_names = sapply(de_list, function (dn) sample(rownames(inmatrix), dn))

    # count the number of occurrences of each identifier
    de_occ = rle(sort(unlist(de_names)))

    if (is.na(N)) {
      N = length(de_list)
    }
    result = sum(de_occ$lengths >= N)
    results = c(results, result)
  }
  return(results)
}
