#'build_predictor
#'
#'
#'

build_predictor = function(fit, data, level=NULL, signature=NULL, genesymbols_from=TRUE, genesymbols_to=NULL, browser=F) {
  if (is.null(signature)) {
    coeffs = fit$b.corrector[[level]]
    #coeffs = fit$b.predictor[[level]]
    ids = names(coeffs)
  }
  else {
    coeffs = signature
    ids = names(signature)
  }
  if (is.null(coeffs)) return(NA)

  inters = ids[grep("\\*", ids)]
  inters_coeffs = coeffs[inters]
  singles = ids[grep("\\*", ids, invert=TRUE)]
  singles_coeffs = coeffs[singles]

  # first process single variables
  idx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][2])
  udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])

  if (!is.null(genesymbols_from)) {
    # these are _levels_ for some odd reason
    #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
    #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
    idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])

    names(singles_coeffs) = idx
    names(udx) = idx
    idx = idx[which(idx %in% colnames(data[,1]))]
    udx = udx[which(idx %in% colnames(data[,1]))]
  }
  #browser()
  udx[which(udx == "dn")] = -1
  udx[which(udx == "up")] = 1

  sigd = data$x[,idx]
  sig = sapply(rownames(sigd), function (r) sum(sigd[r,]*singles_coeffs[idx]*as.integer(udx[idx]), na.rm=T))

  # now process interactions
  for (intid in inters) {
    inters_coeff = inters_coeffs[intid]
    ids = strsplit(intid, "\\*")[[1]]

    idx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][2])
    idxp = idx
    udx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][1])

    if (!is.null(genesymbols_from)) {
      #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
      #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
      idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
      names(udx) = idx
      idx = idx[which(idx %in% colnames(data[,1]))]
      udx = udx[which(idx %in% colnames(data[,1]))]
    }
    udx[which(udx == "dn")] = -1
    udx[which(udx == "up")] = 1

    if (length(idx)<2) next
    sigd = data$x[,idx]
    pdt = rep(1, nrow(sigd))

    # calculate column products
    for (i in 1:ncol(sigd)) {
      x = sigd[,i]
      x[which(is.na(x))] = 0
      pdt = pdt * x * as.integer(udx[i])
    }
    # add to model
    sig = sig + (pdt*inters_coeff)
  }
  if (browser) browser()
  return(sig)
}
