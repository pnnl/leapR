% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leapR.R
\name{leapR}
\alias{leapR}
\title{leapR}
\usage{
leapR(geneset, enrichment_method, ...)
}
\arguments{
\item{geneset}{is a list of four vectors, gene names, gene descriptions, gene sizes and a matrix of genes. It represents .gmt format pathway files.}

\item{enrichment_method}{is a character string specifying the method of enrichment to be performed, one of: "enrichment_comparison", "enrichment_in_order", "enrichment_in_sets", "enrichment_in_pathway", "correlation_enrichment".}

\item{...}{further arguments}
}
\value{
data frame with results
}
\description{
leapR is a wrapper function that consolidates multiple enrichment methods.
}
\details{
Further arguments and enrichment method optional argument information: \cr
\tabular{ll}{
eset \tab Is an \code{ExpressionSet} of expression data, with features as rows and \emph{n} sample/conditions as columns. 
The \code{Annotation} field ideally describes the data type (i.e. proteomics, phosphoproteomics), the \code{featureData} field describes any mapping 
identifiers and the \code{phenoData} field describes any phenotyptic data. We recommend that the `Annotation` slot contain the omics data type for when using with `combine_omics`
This is an required for all active enrichment methods with the exception of 'enrichment_in_relationships'. \cr
\cr
id_column \tab Is a character string, present in the \code{featureData} slot, that is used to specify a column for identifiers to map to enrichment libraries. 
If missing, the rownames of the ExpressionSet will be used. 
\cr
primary_columns \tab Is a character vector composed of column names from \code{eset} (either in the `exprs` or in the `featureData`), 
that specifies a set of primary columns to calculate enrichment on. 
The meaning of this varies according to the enrichment method used - see the descriptions for each method below. 
This is an optional argument used with 'enrichment_in_order', 'enrichment_in_sets', and 'enrichment_comparison' methods. \cr
\cr
secondary_columns \tab Is a character vector of column names for comparison, pulled from the `exprs` of the ExpressionSet. This is an optional argument used with 'enrichment_comparison' methods. \cr
\cr
threshold \tab Is a numeric value, an optional argument used with 'enrichment_in sets' method which filters out abundance values or p-values (depending on what `primary_columns` is used) 
either above or below it. \cr
\cr
greaterthan \tab Is a logical value that defaults to TRUE, it's used with 'enrichment_in_sets' method. 
When set to TRUE, genes with `primary_columns` value above the \code{threshold} argument are kept. 
When set to FALSE genes with `primary_columns` value below the \code{threshold} argument are kept. 
This is an optional argument used with 'enrichment_in_sets' method. \cr
\cr
minsize \tab Is a numeric value, an optional argument used with 'enrichment_in_sets' and 'enrichment_in_order". \cr
\cr
idmap \tab Is...??. This is an optional argument used with 'enrichment_in_relationships' method. \cr
\cr
fdr \tab A numerical value which specifies how many times to randomly sample genes to calculate an empirical false discovery rate, is an optional argument used with 'enrichment_comparison' method. \cr
\cr
min_p_threshold \tab Is a numeric value, a lower p-value threshold and is an optional argument used with 'enrichment_comparison' method. \cr
\cr
sample_n \tab Is a way to subsample the number of components considered for each calculation randomly. This is an optional argument used with 'enrichment_comparison' method. \cr
\cr
}

\strong{Enrichment Methods:}
\cr
\cr
enrichment_comparison
\cr
Compares the distribution of abundances between two sets of conditions for each pathway using a t test. For each pathway in \code{geneset}
uses a t test to compare the distribution of abundance values/numbers in \code{eset} \code{primary_columns} with those in
\code{eset} \code{secondary_columns}. Lower p-values for pathways indicate that the expression of the pathway is
significantly different between the set of conditions in primary_columns and the set of conditions in secondary_columns.
Optionally, users can specify \code{fdr} which will calculate an empirical p-value by randomizing abdunances
\code{fdr} number of times. If the \code{min_p_threshold} is specified the method will only return pathways with an
adjusted p-value lower than the specified threshold. If \code{sample_n} is specified the method will subsample the 
pathway members to the specified number of components.
\cr \cr
enrichment_in_order
\cr
Calculates enrichment of pathways based on a ranked list using the Kologmorov-Smirnov test. 
For each pathway in \code{geneset} uses a Kolgmorov-Smirnov test for rank order to test if the distribution
of ranked abundance values in the \code{eset} \code{primary_columns} is significant relative to a random
distribution. Note that currently \code{primary_columns} only accepts a single column for this method.
\cr \cr
enrichment_in_sets
\cr
Calculates enrichment in pathway membership in a list (e.g. highly differential proteins) relative to background using Fisher's exact test.
For each pathway in \code{geneset} uses a Fisher's exact test over- or under- representation of a list
of components specified. If \code{targets} are specified this must be a vector of identifiers to serve
as the target list for comparison. If \code{eset} and \code{primary_columns} are specified then \code{threshold} specifies
a threshold value for determining the target list of components to test. Specifying \code{greaterthan} to be False
will result in components with values lower than the specified \code{threshold}. If \code{eset} is
a data frame or matrix, the background used for calculation will be taken as the rownames of \code{eset}
\cr \cr
enrichment_in_pathway
\cr
Compares the distribution of abundances in a pathway with the background distribution of abundances using a t test
For each pathway in \code{geneset} calculates the signficance of the difference between the abundances
from pathway members versus abundance of non-pathway members in the set of conditions specified by \code{primary_columns}.
Optionally, users can specify \code{fdr} which will calculate an empirical p-value by randomizing abdunances
\code{fdr} number of times. If the \code{min_p_threshold} is specified the method will only return pathways with an
adjusted p-value lower than the specified threshold. If \code{sample_n} is specified the method will subsample the 
pathway members to the specified number of components.
\cr \cr
correlation_enrichment
Calculates the enrichment of a pathway based on correlation between pathway members across conditions versus correlation between members not in the pathway.
For each pathway in \code{geneset} calculates the pairwise correlation between all pathway members and non-pathway members
across the specified \code{primary_columns} conditions in \code{eset}. Note that for large matrices this can take a long
time. A p-value is calculated based on comparing the correlation within the members of a pathway with the correlation
values between members of the pathway and non-members of the pathway.
\cr
}
\examples{
        library(leapR)

        # read in the example abundance data
        # read in the example transcriptomic data
        tdata <- download.file("https://api.figshare.com/v2/file/download/55781153",method='libcurl',destfile='transData.rda')
        load('transData.rda')
        p <- file.remove("transData.rda")

        # read in the pathways
        data("ncipid")

        # read in the patient groups
        data("shortlist")
        data("longlist")
        
        # use enrichment_comparison to calculate enrichment in one set of conditions (shortlist) and another
        # (longlist)
        short_v_long = leapR(geneset=ncipid, enrichment_method='enrichment_comparison', 
              eset=tset, primary_columns=shortlist, secondary_columns=longlist)
        
        # use enrichment_in_sets to calculate the most enriched pathways from the highest abundance proteins
        #     from one condition
        onept_sets = leapR(geneset=ncipid, enrichment_method='enrichment_in_sets',
               eset=tset, primary_columns="TCGA-13-1484", threshold=1.5)
               
         # use enrichment_in_order to calculate the most enriched pathways from the same condition
         #     Note: that this uses the entire set of abundance values and their order - whereas
         #     the previous example uses a hard threshold to get a short list of most abundant proteins
         #     and calculates enrichment based on set overlap. The results are likely to be similar - but
         #     with some notable differences.
         onept_order = leapR(geneset=ncipid, enrichment_method='enrichment_in_order',
               eset=tset, primary_columns="TCGA-13-1484")
               
         # use enrichment_in_pathway to calculate the most enriched pathways in a set of conditions
         #     based on abundance in the pathway members versus abundance in non-pathway members
         short_pathways = leapR(geneset=ncipid, enrichment_method='enrichment_in_pathway',
               eset=tset, primary_columns=shortlist)
               
         # use correlation_enrichment to calculate the most enriched pathways in correlation across
         #     the shortlist conditions
         short_correlation_pathways = leapR(geneset=ncipid, enrichment_method='correlation_enrichment',
                eset=tset, primary_columns=shortlist)


}
