#'random_gene_set_overlap_db
#'
#'
#'
#'

random_gene_set_overlap_db <- function(inmatrix_1, inmatrix_2, n1, n2, rep=100) {
  # This is a variant of the above to test 2 lists of different lengths
  #      to see how much of an overlap there is- and replicate randomly a bunch of times

  if (is.matrix(inmatrix_1)) {
    names1 = rownames(inmatrix_1)
  } else names1 = names(inmatrix_1)

  if (is.matrix(inmatrix_2)) {
    names2 = rownames(inmatrix_2)
  } else names2 = names(inmatrix_2)


  results = c()

  for (i in 1:rep) {
    sn1 = sample(names1, n1)
    sn2 = sample(names2, n2)
    x = sum(sn1 %in% sn2)

    results = c(results, x)
  }
  return(results)
}
