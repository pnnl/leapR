#'get_pathway_information
#'
#'

get_pathway_information <- function(geneset, path, remove.tags=F) {
  i = which(geneset$names == path)
  thissize = geneset$size[i]
  thisdesc = geneset$desc[i]

  # default is to look through all members for effects on the pathway
  grouplist = geneset$matrix[i,1:thissize]
  grouplist = unique(sort(grouplist[which(grouplist != "")]))

  if (remove.tags) {
    grouplist = unique(sapply(grouplist, function (n) strsplit(n, "_")[[1]][2]))
  }
  thissize = length(grouplist)

  return(list(name=path, size=thissize, description=thisdesc, geneset=grouplist))
}
