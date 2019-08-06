#'geneset_pca
#'
#'
#'



geneset_pca <- function(geneset, datatable, minsize=5, ncomp=1) {
  out = list()
  for (c in 1:ncomp) {
    results = matrix(nrow=length(geneset$names), ncol=ncol(datatable))
    rownames(results) = geneset$names
    colnames(results) = colnames(datatable)

    for (i in 1:length(geneset$names)) {
      thisname = geneset$names[i]
      thissize = geneset$size[i]
      thisdesc = geneset$desc[i]
      grouplist = geneset$matrix[i,1:thissize]

      in_group = datatable[grouplist[which(grouplist %in% rownames(datatable))],]

      if (nrow(in_group) < minsize) next

      prc = prcomp(in_group)
      #return(list(thisname, in_group, prc))

      results[thisname,] = prc$rotation[,c]
    }
    out = c(out, list(results))
  }
  return(out)
}
