#'average_redundant_rows
#'
#'

average_redundant_rows <- function(data) {
  thesenames = unique(rownames(data))
  uniquedata = matrix(0, ncol=ncol(data), nrow=length(thesenames))
  for(i in 1:length(thesenames)) {
    lines = which(rownames(data) == thesenames[i])
    uniquedata[i,] = sapply(colnames(data), function (i) mean(data[lines,i], na.rm=T))
  }
  rownames(uniquedata) = thesenames
  colnames(uniquedata) = colnames(data)
  return(uniquedata)
}
