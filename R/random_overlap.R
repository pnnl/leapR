#'random_overlap
#'
#'

random_overlap <- function(items, counts=c(10,50,100), N=3,rep=100) {
  # this just returns the number of items that match when randomly drawn
  #     with the numbers indicated in counts
  results = c()

  for (i in 1:rep) {
    # now randomly select genes
    bats = c()

    for (this in counts) {
      bits = sample(items, this)
      bats = c(bats, bits)
    }

    # count the number of occurrences of each item
    de_occ = rle(sort(bats))

    result = sum(de_occ$lengths >= N)
    results = c(results, result)
  }
  return(results)
}
