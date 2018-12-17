# ClusterEnrichment

cluster_enrichment <- function(hclust, annotation, cutrange=2:10) {
	results = list()
	for (i in cutrange) {
		cc = cutree(hclust, k=i)
		res = c()
		rest = list()
		for (j in 1:i) {
			f = enrichment_by_fishers(names(which(cc==j)), names(cc), annotation)
			rest = c(rest, list(f))
			res = c(res, f$mat[1,1])
		}
		x = rest[which(res==max(res))]
		results = c(results, list(x))
	}
	return(results)
}