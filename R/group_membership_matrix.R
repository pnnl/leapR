#'group_membership_matrix
#'
#'
#'

group_membership_matrix <- function(geneset, pathwayname, sources, mode="original") {
  # for the named pathway creates a matrix that indicates membership
  #     from each source where source is a list of gene/protein lists
  members = geneset$matrix[which(geneset$names==pathwayname),1:geneset$sizes[which(geneset$names==pathwayname)]]

  results = matrix(nrow=length(members), ncol=length(sources))
  rownames(results) = members
  colnames(results) = names(sources)

  for (i in 1:length(sources)) {
    if (mode == "original") results[,i] = members %in% sources[i][[1]]
    else {
      results[members[which(members %in% names(sources[i][[1]]))],i] = sources[i][[1]][members[which(members %in% names(sources[i][[1]]))]]
    }
  }
  return(results)
}
