#' get_pathway_information
#'
#' get_pathway_information extracts information about a pathway from a GeneSet object
#'
#' @param geneset is a GeneSet object for pathway annotation
#' @param name is the name of the gene set to be returned
#' 
#' @examples
#' dontrun{
#'      library(leapr)
#'      
#'      # load example gene set
#'      data("ncipid")
#'      
#'      tnfpathway = get_pathway_information(ncipid, "tnfpathway")
#'
#' }
#'
#' @export
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
