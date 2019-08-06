#'bottleneck_overlap
#'
#'
#'

bottleneck_overlap <- function(set_1, set_2, percentage=0.2) {
  # takes two matrices/data frames that are ordered and returns a
  # matrix of overlapping rows by rowname
  ret_set = set_1[which(rownames(set_1)[1:(nrow(set_1)*percentage)] %in% rownames(set_2)[1:(nrow(set_2)*percentage)]),]
  return(ret_set)
}
