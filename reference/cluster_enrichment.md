# cluster_enrichment

Cluster enrichment Run enrichment (Fisher's exact) on clusters (lists of
identifier groups)

## Usage

``` r
cluster_enrichment(eset, assay_name, geneset, clusters, sigfilter = 0.05)
```

## Arguments

- eset:

  is an SummarizedExperiment containing data that is clustered

- assay_name:

  is the name of the assay

- geneset:

  is a GeneSet object for pathway annotation

- clusters:

  is a list of clusters (gene lists) to calculate enrichment on,
  generally the result of the \`cutree\` function

- sigfilter:

  minimum significance threshold default is .05

## Value

data frame with enrichment results

## Details

This function will calculate enrichment (Fisher's exact test for
membership overlap) on

a series of lists of genes, such as from a set of clusters. The results
are returned as

a list of results matrices in the order of the input clusters.

## Examples

``` r
        library(leapR)

        # read in the example transcriptomic data
        url <- "https://api.figshare.com/v2/file/download/56536214"
        tdata <- download.file(url,method='libcurl',destfile='transData.rda')
        load('transData.rda')
        p <- file.remove("transData.rda")

        # read in the pathways
        data("ncipid")

        # for the example we will limit the number of transcripts considered
        #- arbitrarily in this case
        transdata <- SummarizedExperiment::assay(tset,'transcriptomics')
        transdata[which(is.na(transdata),arr.ind=TRUE)]<-0.0
        # perform heirarchical clustering on the  data
        transdata.hc <- hclust(dist(transdata), method="ward.D2")

        transdata.hc.clusters <- cutree(transdata.hc, k=5)
        clust.list <- lapply(seq_len(5), function(x) {
           return(names(which(transdata.hc.clusters==x)))})
        #calculates enrichment for each of the clusters individually a
        #and returns a list of enrichment results
        transdata.hc.enrichment <- leapR::cluster_enrichment(eset=tset,
                assay_name='transcriptomics',
                geneset=ncipid,
                clusters=clust.list)


```
