#' read_gene_sets
#'
#' read_gene_sets is a function to import external pathway database files in .gmt format
#' @importFrom readr  read_lines
#' @param gsfile is a gene set file, for example a .gmt file (gene matrix transposed file format)
#' @param gene.labels defaults to NA
#' @param gs.size.threshold.min defaults to 5
#' @param gs.size.threshold.max defaults to 15000
#' @return gene set object
#' @export
#' @return geneset list
#' @examples
#' gfile <- system.file('extdata','h.all.v2024.1.Hs.symbols.gmt',package='leapR')
#' glist <- read_gene_sets(gfile)
#'

read_gene_sets <- function(gsfile, gene.labels=NA, gs.size.threshold.min=5, gs.size.threshold.max=15000) {
  # Read input gene set database

  temp <- readr::read_lines(gsfile)

  max.Ng <- length(temp)
  temp.size.G <- vapply(seq_along(1:max.Ng), function(i){ #vector(length = max.Ng, mode = "numeric")
    length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }, double(1))

  max.size.G <- max(temp.size.G)
  gs_stats <- vapply(seq_along(1:max.Ng),function(i){
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[2]
    gene.set.tags <- vapply(seq_along(1:gene.set.size), function(j) gs.line[j + 2], "")
      # we just want to load it all (unless we want to filter here?)
    existing.set <- rep(TRUE, length(gene.set.tags))

    set.size <- length(existing.set[existing.set == TRUE])
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next

    return(list(gvals = c(gene.set.tags[existing.set], rep(NA, max.size.G - temp.size.G[i])),
                gsnames = gene.set.name,
                gsdesc = gene.set.desc,
                gsize = set.size))
  }, list(gvals = list(),gsnames = "",gsdesc = "",gsize = numeric(1)))

  gs.names <- apply(gs_stats, 2,function(g) return(g$gsnames))
  gs.desc <- apply(gs_stats,2,function(g) return(g$gsdesc))
  sizes <- apply(gs_stats, 2, function(g) return(g$gsize))
  gs <- apply(gs_stats,2, function(g) return(g$gvals))
  result <- list(names = gs.names, desc = gs.desc, sizes = sizes, matrix = gs)

  class(result) <- c("geneset_data", "list")
  return(result)
}
