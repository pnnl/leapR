% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcTTest.R
\name{calcTTest}
\alias{calcTTest}
\title{calcTTest}
\usage{
calcTTest(eset, group1, group2)
}
\arguments{
\item{eset}{ExpressionSet}

\item{group1}{List of samples comprising group 1}

\item{group2}{List of samples comprising group 2}
}
\value{
An Expression set with two colunns added to the featureData slot: pvalue, and estimate
}
\description{
calculates a t-test for two distributions of data on a per-gene basis
append results to ExpressionSet with two extra columsn: `pvalue` and `difference` for each feature
}
\examples{

        library(leapR)
        tdata <- download.file("https://api.figshare.com/v2/file/download/55781153",method='libcurl',destfile='transData.rda')
        load('transData.rda')
        p <- file.remove("transData.rda")

        # read in the pathways
        data("ncipid")

        # read in the patient groups
        data("shortlist")
        data("longlist")
        calcTTest(tset, shortlist, longlist) 
}
