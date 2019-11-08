#' pairwise_overlap_enrichment
#'
#' pairwise_overlap_enrichment function description is...
#'
#' @param geneset is...
#' @param datamatrix is...
#' @param condition1 is...
#' @param condition2 is...
#' @param mapping_column defaults to NULL
#' @param subsample_components defaults to NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

pairwise_overlap_enrichment <- function(geneset, datamatrix, condition1, condition2, mapping_column=NULL,
                                        subsample_components=NULL) {
  # iterate through a geneset and calculate the enrichment (abundance in condition1 v condition2)
  #      in the overlapping and unique gene sets for every combination

  results = list()

  overlap = matrix(nrow=length(geneset$names), ncol=length(geneset$names))
  rownames(overlap) = geneset$names
  colnames(overlap) = geneset$names

  # this wouldn't work in Python - I think it does in R
  # - that is the resulting matrices are new
  uniquea = unique_count = overlap_count = overlap

  for (path1 in geneset$names) {
    #cat(path1, "\n")
    for (path2 in geneset$names) {
      if (path1 != path2) {
        this = overlap_gene_sets(geneset, path1, path2)

        if (is.null(mapping_column)) {
          theserows_o = rownames(datamatrix)[which(rownames(datamatrix) %in% this$overlap)]
          theserows_a = rownames(datamatrix)[which(rownames(datamatrix) %in% this$set1_unique)]
        }
        else {
          theserows_o = rownames(datamatrix)[which(datamatrix[,mapping_column] %in% this$overlap)]
          theserows_a = rownames(datamatrix)[which(datamatrix[,mapping_column] %in% this$set1_unique)]
        }

        if (!is.null(subsample_components)) {
          if (length(theserows_o)>subsample_components) theserows_o = sample(theserows_o, subsample_components)
          else theserows_o = c()
          if (length(theserows_a)>subsample_components) theserows_a = sample(theserows_a, subsample_components)
          else theserows_a = c()
        }
        #return(theserows_b)
        o = a = b = NA
        if (length(theserows_o)>5) o = try(t.test(datamatrix[theserows_o, condition1], datamatrix[theserows_o, condition2])$p.value, silent=T)
        if (length(theserows_a)>5) a = try(t.test(datamatrix[theserows_a, condition1], datamatrix[theserows_a, condition2])$p.value, silent=T)
        if (!class(o)=="try-error") overlap[path1,path2] = o
        if (!class(a)=="try-error") uniquea[path1,path2] = a

        unique_count[path1,path2] = length(theserows_a)
        overlap_count[path1,path2] = length(theserows_o)

        #cat(path1, path2, o, a, b, "\n", sep="\t")
      }
    }
  }
  return(list(overlap_pmat=overlap, unique_pmat=uniquea,
              overlap_counts=overlap_count, unique_counts=unique_count))
}
