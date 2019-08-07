#' plot_expression
#'
#' plot_expression function description is...
#'
#' @param ratios is...
#' @param rowname is...
#' @param pfactor defaults to 1
#' @param xvector defaults to NA
#' @param xfactor defaults to 0
#' @param title defaults to NA
#' @param condrange is...
#' @param trendline defaults to FALSE
#' @param trendcolor defaults to 'red'
#' @param color defaults to 'grey'
#' @param legend defaults to FALSE
#' @param yr defaults to NA
#' @param mode defaults to "new"
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

plot_expression <- function(ratios, rowname, pfactor=1, xvector=NA, xfactor=0, title=NA, condrange=c(), trendline=FALSE, trendcolor="red", color="grey", legend=FALSE, yr=NA, mode="new") {
  # expand the y margin lower
  par(mar=c(15,4.1,4.2,2.1))

  #if (!class(title) == "character") {
  #	if (class(rowname) == "character") {
  #		title = rowname
  #	} else {
  #		title = rownames(ratios)[rowname]
  #	}
  #}
  if (!is.vector(rowname)) {
    rowname = c(rowname)
    if (length(condrange) == 0) {
      condrange = c(1:length(ratios[rowname,]))
    }
  } else {
    if (length(condrange) == 0) {
      condrange = c(1:length(ratios[rowname[1],]))
    }
  }

  name= rowname[1]
  if (mode == "new") {
    if (is.na(yr)) {
      if (is.na(xvector)) {
        #cat(name, str(ratios[name, condrange]),"\n")
        # plot the observed values in red with no axes
        plot(ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1)
      } else {
        plot(xvector, ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1)
      }

    } else {
      if (is.na(xvector)) {
        plot(ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1, ylim=yr)
      } else {
        plot(xvector, ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1, ylim=yr)
      }
    }
  }
  for (name in rowname[2:(length(rowname))]) {
    if (is.na(xvector)) {
      lines(ratios[name,condrange], type="l", col=color)
    } else {
      lines(xvector, ratios[name,condrange], type="l", col=color)
    }
  }
  if (trendline == TRUE) {
    trend = sapply(condrange, function (c) mean(as.numeric(ratios[rowname,c])))
    if (is.na(xvector)) {
      lines(trend, type="l", col=trendcolor)
    } else {
      lines(xvector, trend, type="l", col=trendcolor)
    }
  }
  # add categories on the x axis
  if (is.na(xvector)) {
    axis(1, labels=names(ratios[name,condrange]), at=1:length(ratios[name,condrange]), las=2)
  } else {
    axis(1, labels=names(ratios[name,condrange]), at=xvector, las=2)
  }
  axis(2)
}
