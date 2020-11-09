#' Protein expression data
#' 
#' A dataset containing protein expression data containing 1999 proteins and 174 patient samples. 
#' Values represent the logratio of the protein expression compared to a control
#' @format A data frame with 1999 rows and 174 columns
#' @source H. Zhang et al., Integrated Proteogenomic Characterization of Human High-Grade Serous Ovarian Cancer. Cell 166, 755-765 (2016)
"protdata"

#' Transcriptomic data
#' 
#' A dataset containing gene transcriptomic data containing 18632 genes and 174 patient samples. 
#' Values represent normalized (Z score) gene expression relative to the mean expression for that gene
#' @format A data frame with 18632 rows and 174 columns
#' @source H. Zhang et al., Integrated Proteogenomic Characterization of Human High-Grade Serous Ovarian Cancer. Cell 166, 755-765 (2016)
"transdata"

#' Phosphoproteomics data
#' 
#' A dataset containing phosphoproteomic data containing 20732 phosphosites and 69 patient samples. 
#' Values represent normalized (Z score) phosphorylation adjusted to the protein abundance
#' @format A data frame with 20732 rows and 70 columns - first column is protein id
#' @source H. Zhang et al., Integrated Proteogenomic Characterization of Human High-Grade Serous Ovarian Cancer. Cell 166, 755-765 (2016)
"transdata"

#'NCI Gene lists
#'
#'A list of pathways and the genes that comprise these pathways
#'@format A list with 4 items
#'\describe{
#' \item{names} The names of the signling pathways
#' \item{desc}
#' \item{sizes}
#' \item{matrix}
#' }
#' @source NCIPID
#'
"ncipid"

#' MSigDB Gene Lists
#' 
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the signaling pathways
#' \item{desc}
#' \item{sizes}
#' \item{matrix}
#' }
#' @source MSIGDB
#' 
"msigdb"

#' A list of pathways and genes that comprise these pathways from msigdb
#' @format a list with 4 items
#' 
#' Short list of patient samples
#' 
"shortlist"
#' Long list of patient samples
#' 
"longlist"

