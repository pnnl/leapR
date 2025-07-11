#' read_gene_sets
#'
#' read_gene_sets is a function to import external pathway database files in .gmt format
#' @import readr
#' @param gsfile is a gene set file, for example a .gmt file (gene matrix transposed file format)
#' @param gene.labels defaults to NA
#' @param gs.size.threshold.min defaults to 5
#' @param gs.size.threshold.max defaults to 15000
#' @return gene set object
#' @export
#' @return geneset list
#' @examples
#' gfile = system.file('extdata','h.all.v2024.1.Hs.symbols.gmt',package='leapR')
#' glist = read_gene_sets(gfile)
#'

read_gene_sets <- function(gsfile, gene.labels=NA, gs.size.threshold.min=5, gs.size.threshold.max=15000) {
  # Read input gene set database
  
  temp <- readr::read_lines(gsfile)

  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
  }

  max.size.G <- max(temp.size.G)
  gs <- matrix(rep(NA, max.Ng*max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
    gene.set.name <- gs.line[1]
    gene.set.desc <- gs.line[2]
    #cat(gene.set.name, "\n")
    #cat(gs.line)
    gene.set.tags <- vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j] <- gs.line[j + 2]
    }
    if (is.na(gene.labels)) {
      # we just want to load it all (unless we want to filter here?)
      existing.set = rep(TRUE, length(gene.set.tags))
    } else {
      existing.set <- is.element(gene.set.tags, gene.labels)
    }

    set.size <- length(existing.set[existing.set == TRUE])
    if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
    temp.size.G[gs.count] <- set.size
    gs[gs.count,] <- c(gene.set.tags[existing.set], rep(NA, max.size.G - temp.size.G[gs.count]))
    temp.names[gs.count] <- gene.set.name
    temp.desc[gs.count] <- gene.set.desc
    gs.count <- gs.count + 1
  }
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]

  result = list(names = gs.names, desc = gs.desc, sizes = size.G, matrix = gs)

  class(result) = c("geneset_data", "list")
  return(result)
}
