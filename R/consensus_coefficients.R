#' consensus_coefficients
#'
#' consensus_coefficients function description is...
#'
#' @param fitcv is...
#' @param level is...
#' @param return.all default is TRUE
#' @param filter_count default is NULL
#' @param meanweight default is NULL
#' @param subsample default is NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

consensus_coefficients = function(fitcv, level, return.all=T, filter_count=NULL, meanweight=NULL, subsample=NULL) {
  # get a list of the lengths of each model- then we can screen for which have the
  # appropriate number of levels (some are shorter than max)
  if (!is.null(subsample)) {
    usethese = sample(1:length(fitcv), subsample)
  }
  else {
    usethese = 1:length(fitcv)
  }

  fitlengths = sapply(usethese, function (i) length(fitcv[[i]][[1]]$b.corrector))

  passed = usethese[which(fitlengths>=level)]
  if (sum(fitlengths>=level)==0) return(NA)

  #varsrle = rle(sort(unlist(sapply(1:length(fitcv), function (i) names(fitcv[[i]][[1]]$b.corrector[[level]])))))
  varsrle = try(rle(sort(unlist(sapply(usethese[which(fitlengths>=level)],
                                       function (i) names(fitcv[[i]][[1]]$b.corrector[[level]]))))))

  if (class(varsrle) == "try-error") browser()

  #browser()
  results = matrix(nrow=length(fitlengths), ncol=length(varsrle$values))
  colnames(results) = varsrle$values

  #browser()
  for (i in 1:length(usethese)) {
    x = usethese[i]

    if (level > length(fitcv[[x]][[1]]$b.corrector)) next()

    for (j in 1:length(fitcv[[x]][[1]]$b.corrector[[level]])) {
      val = fitcv[[x]][[1]]$b.corrector[[level]][j]
      results[i,names(val)] = val
    }
  }
  meansig = colMeans(results, na.rm=T)
  names(meansig) = colnames(results)

  if (!is.null(meanweight)) {
    meansig = meansig*varsrle$lengths
  }

  if (!is.null(filter_count)) {
    meansig = meansig[which(varsrle$lengths>filter_count)]
  }

  #browser()

  if (return.all) {
    results_annotated = results
    colnames(results_annotated) = line.to.gene.symbol(colnames(results), fitcv[[level]][[1]]$gene.symbols)
    meansig_annotated = meansig
    names(meansig_annotated) = line.to.gene.symbol(names(meansig_annotated), fitcv[[level]][[1]]$gene.symbols)
    return(list(mat=results, meansig=meansig, counts=varsrle$lengths,
                meansig_annotated=meansig_annotated, mat_annotated=results_annotated))
  }

  return(meansig)
}
