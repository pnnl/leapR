#' build_predictor_save
#'
#' build_predictor_save function description is...
#'
#' @param fit is...
#' @param data is...
#' @param level default is NULL
#' @param signature default is NULL
#' @param genesymbols_from default is TRUE
#' @param genesymbols_to default is NULL
#' @param browser default is FALSE
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

build_predictor_save = function(fit, data, level=NULL, signature=NULL, genesymbols_from=TRUE, genesymbols_to=NULL, browser=F) {
  if (is.null(signature)) {
    coeffs = fit$b.corrector[[level]]
    #coeffs = fit$b.predictor[[level]]
    ids = names(coeffs)
  }
  else {
    coeffs = signature
    ids = names(signature)
  }
  inters = ids[grep("\\*", ids)]
  singles = ids[grep("\\*", ids, invert=TRUE)]

  # first process single variables
  idx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][2])
  udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])

  if (!is.null(genesymbols_from)) {
    # these are _levels_ for some odd reason
    #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
    #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
    idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])

  }
  browser()
  udx[which(udx == "dn")] = 1
  udx[which(udx == "up")] = -1

  sigd = data$x[,as.integer(idx)]
  sig = sapply(rownames(sigd), function (r) sum(sigd[r,]*coeffs[singles]*as.integer(udx), na.rm=T))

  # now process interactions
  for (intid in inters) {
    coeff = coeffs[intid]
    ids = strsplit(intid, "\\*")[[1]]

    idx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][2])
    udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])

    if (!is.null(genesymbols_from)) {
      #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
      #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
      idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
      idx = which(colnames(data[,1]) %in% idx)
    }
    udx[which(udx == "dn")] = 1
    udx[which(udx == "up")] = -1

    sigd = data$x[,as.integer(idx)]
    pdt = rep(1, nrow(sigd))

    # calculate column products
    for (i in 1:ncol(sigd)) {
      x = sigd[,i]
      x[which(is.na(x))] = 0
      pdt = pdt * x * as.integer(udx[i])
    }

    # add to model
    sig = sig + (pdt*coeff)
  }

  if (browser) browser()
  return(sig)
}
