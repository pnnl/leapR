#'cluster_enrichment
#'
#'
#'

cluster_enrichment <- function(database, clusters, background=NA, sigfilter=0.05) {
  x = length(clusters)-1
  if (is.na(background))  background = clusters[[x+1]]
  this = sapply(1:x, function (i) list(enrichment_in_groups(database, clusters[[i]], background)))
  outlist = list()
  for (i in 1:x) {
    these = this[[i]][which(this[[i]][,7]<sigfilter),]
    if (length(which(this[[i]][,7]<sigfilter)) == 1) {
      these = as.matrix(t(these))
      rownames(these)[[1]] = rownames(this[[i]])[which(this[[i]][,7]<sigfilter)]
    }
    outlist = c(outlist, list(these))
  }
  return(outlist)
}
