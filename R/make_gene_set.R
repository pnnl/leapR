#' make_gene_set
#'
#' make_gene_set function description is...
#'
#' @param name is...
#' @param desc is...
#' @param genelist is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

make_gene_set <- function(name, desc, genelist) {
  return(list(names=c(name,"dum"),
              desc=c(desc,"dum"),
              size=c(length(genelist),1),
              matrix=t(data.frame(genelist, rep("null", length(genelist)), ncol=length(genelist)))))
}
