#' enrichment_by_fishers
#'
#' Calculate statistical enrichment by membership (Fisher's exact test) between a group and background
#' based on an annotation.
#'
#' @param group is a list of identifiers (e.g. differentially expressed gene ids) to test for enrichment
#' @param background is the total list of identifiers to consider as a background
#' @param annotation is a list of identifiers to compare with group (e.g. a pathway)
#'
#' @details This function uses the Fisher's exact test to test a contingency table formed by the group, the annotation,
#' @details and the background - a lists of identifiers. Also known as a hypergeometric test, this essentially tests
#' @details how surprising it is that one group (group) overlaps with another group (annotation) in the context of the
#' @details background of all identifiers considered (background; e.g. all the proteins identified in the assay), to the
#' @details extent observed in a dataset.
#'
#' @examples
#' dontrun{
#' 
#'  background = c("A","B","C","D","E","F","G","H","I","J","K","L")
#'  group = c("A","B","C","D")
#'  annotation = c("A", "B","C","D","E")
#'  
#'  result = enrichment_by_fishers(group, background, annotation)
#'  result$fisher$p.value
#'
#'
#' }
#'
#' @export
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
