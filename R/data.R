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
#' \item{desc} Short description of the pathways
#' \item{sizes} Number of genes in the signaling pathways
#' \item{matrix} Matrix containing the genes in the pathways
#' }
#' @source NCIPID
#'
"ncipid"

#' KEGG, Reactome, BioCarta Pathways
#'
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the pathways
#' \item{desc} Short description of the pathways
#' \item{sizes} Number of genes in the signaling pathways
#' \item{matrix} Matrix containing the genes in the pathways
#' }
#' @source MSIGDB
#'
"krbpaths"

#' Multi-Omic KEGG, Reactome, BioCarta Pathways
#'
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the pathways
#' \item{desc} Short description of the pathways
#' \item{sizes} Number of genes in the signaling pathways
#' \item{matrix} Matrix containing the genes in the pathways with multi-omic prefixes
#' }
#' @source MSIGDB
#'
"mo_krbpaths"

#' Kinase substrate lists
#'
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the kinases
#' \item{desc} Short description of the kinase
#' \item{sizes} Length of the substrate list
#' \item{matrix} Substrate list for the kinase
#' }
#' @source PhosphositePlus
#'
"kinasesubstrates"


#' A list of pathways and genes that comprise these pathways from msigdb
#' @format a list with 4 items
#'
#' Short list of patient samples
#'
"shortlist"
#' Long list of patient samples
#'
"longlist"

#' KEGG, Reactome, BioCarta Pathways For Mouse (Gene Symbols)
#'
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the pathways
#' \item{desc} Short description of the pathways
#' \item{sizes} Number of genes in the signaling pathways
#' \item{matrix} Matrix containing the genes in the pathways with gene symbols
#' }
#' @source MSIGDB
#'
"m2_mouse_symbols"

#' KEGG, Reactome, BioCarta Pathways For Mouse (Entrezs)
#'
#' @format A list with 4 items
#'\describe{
#' \item{names} The names of the pathways
#' \item{desc} Short description of the pathways
#' \item{sizes} Number of genes in the signaling pathways
#' \item{matrix} Matrix containing the genes in the pathways with entrez ids
#' }
#' @source MSIGDB
#'
"m2_mouse_entrez"
