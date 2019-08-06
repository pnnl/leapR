#'random_gene_set_overlap_x3
#'
#'
#'
#'

random_gene_set_overlap_x3 <- function(list_1, list_2, list_3, n1, n2, n3, rep=100) {
  # This is a variant of the above to test 3 lists of different lengths
  #      to see how much of an overlap there is- and replicate randomly a bunch of times

  results = c()

  for (i in 1:rep) {
    sn1 = sample(list_1, n1)
    sn2 = sample(list_2, n2)
    sn3 = sample(list_3, n3)
    x1 = sn1[which(sn1 %in% sn2)]
    x = sum(x1 %in% sn3)

    results = c(results, x)
  }
  return(results)
}
