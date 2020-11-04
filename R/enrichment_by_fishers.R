#' enrichment_by_fishers
#'
#' # helper function for enrichment functions
#'
#'

enrichment_by_fishers <- function(group, background, annotation) {
  # calculate the fishers exact on a group of things versus a background
  #     for those things with a particular annotation (on another list)
  non_group = background[!background %in% group]

  group_annot = sum(group %in% annotation)
  group_annot_names = paste(intersect(group, annotation), collapse = ", ")
  non_group_annot = sum(non_group %in% annotation)

  group_nonannot = sum(!group %in% annotation)
  non_group_nonannot = sum(!non_group %in% annotation)

  test = matrix(c(group_annot, non_group_annot, group_nonannot, non_group_nonannot), nr=2,
                dimnames=list(c("Group", "NonGroup"), c("Annotated", "NonAnnotated")))

  per = c(test[1,1]/(test[1,1]+test[1,2]), test[2,1]/(test[2,1]+test[2,2]))

  ft = fisher.test(test)
  test = cbind(test, per)
  #print(test)
  #cat(sum(test), length(background), "\n")

  fold = per[1]/per[2]

  return(list(fisher=ft, mat=test, foldx=fold, in_path_names=group_annot_names))
}
