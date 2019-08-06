#'progressive_enrichment
#'
#'
#'

progressive_enrichment <- function(geneset, rankedlist, decreasing=T, ntimes=10, nstep=10, startat=10, endat=0.5) {
  # Progressive enrichment performs a simple Fisher's test on successively more and more
  #      of a ranked list. At each step the list is randomized ntimes and
  #      the percentage of times that a larger enrichment score is obtained is
  #      summed to provide an overall enrichment score

  rankedlist = sort(rankedlist, decreasing=decreasing)
  endatx = as.integer(length(rankedlist)*endat)
  presults_p = matrix(nrow=length(geneset$names), ncol=length(seq(startat, endatx, nstep)), dimnames=list(geneset$names, seq(startat, endatx, nstep)))
  presults_foldx = presults_p

  x = 1
  for (i in seq(startat, endatx, nstep)) {
    cat("Status", i, "\n")
    result = enrichment_in_groups(geneset, names(rankedlist[1:i]), names(rankedlist))

    randmat = sapply(1:ntimes, function (j) {randlist=sample(names(rankedlist), i); enrichment_in_groups(geneset, randlist, names(rankedlist))[,5]>result[,5]})
    presults_p[,x] = rowSums(randmat)/ntimes
    presults_foldx[,x] = result[,5]
    x = x + 1
  }
  return(list(pvalues=presults_p, foldx=presults_foldx))
}
