

#'NCI Gene lists
#'
#'A list of pathways and the genes that comprise these pathways
#'@format A list with 4 items
#'\describe{
#' \item{names}{The names of the signaling pathways}
#' \item{desc}{Short description of the pathways}
#' \item{sizes}{Number of genes in the signaling pathways}
#' \item{matrix}{Matrix containing the genes in the pathways}
#' }
#' @source NCIPID
#'
"ncipid"

#' KEGG, Reactome, BioCarta Pathways
#' 
#' @format A list with 4 items
#'\describe{
#' \item{names}{The names of the pathways}
#' \item{desc}{Short description of the pathways}
#' \item{sizes}{Number of genes in the signaling pathways}
#' \item{matrix}{Matrix containing the genes in the pathways} 
#' }
#' @source MSIGDB
#' 
"krbpaths"

#' MSIGDB
#' Is this necessary?
"msigdb"

#' Kinase substrate lists
#' 
#' @format A list with 4 items
#'\describe{
#' \item{names}{The names of the kinases}
#' \item{desc}{Short description of the kinase}
#' \item{sizes}{Length of the substrate list}
#' \item{matrix}{Substrate list for the kinase}
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

