#' rpf_coxph_formatdata
#'
#' rpf_coxph_formatdata function description is...
#'
#' @param data is...
#' @param survdata is...
#' @param timevar is...
#' @param statusvar is...
#' @param extra.features defaults to NULL
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

rpf_coxph_formatdata <- function(data, survdata, timevar, statusvar, extra.features=NULL) {
  # format for suvival data analysis using the regularisation path following algorithm
  namemat = as.data.frame(cbind(1:nrow(data), rownames(data)))
  colnames(namemat) = c("ID", "Gene.Symbol")

  matchcols = colnames(data)[which(colnames(data) %in% rownames(survdata))]

  efeatures = NULL
  if (!is.null(extra.features)) {
    efeatures = survdata[matchcols,extra.features]
  }

  stimes = as.numeric(survdata[matchcols,timevar])

  status = as.numeric(survdata[matchcols,statusvar])

  # some of the status or time values are missing and thus the data can't be used- screen these out
  x = which(!is.na(status)&!is.na(stimes))

  matchcols = matchcols[x]

  stimes = stimes[x]
  status = status[x]
  efeatures = efeatures[x,]

  mat = I(t(data[,matchcols]))

  if (is.null(extra.features)) {
    data.struct = list(data.frame(x=mat,
                                  time=stimes,
                                  status=status),
                       namemat)
  }
  else {
    data.struct = list(list(x=mat,
                            time=stimes,
                            status=status,
                            extra.features=efeatures),
                       namemat)
  }

  return(data.struct)
}
