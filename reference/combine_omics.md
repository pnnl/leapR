# combine_omics Combine two or more omics matrices into one multi-omics matrix with 'tagged' ids.

combine_omics Combine two or more omics matrices into one multi-omics
matrix with 'tagged' ids.

## Usage

``` r
combine_omics(omics_list, id_list = rep(NA, length(omics_list)))
```

## Arguments

- omics_list:

  Is a list of `SummarizedExperiment` each with one assay

- id_list:

  List of identifiers to use, in the same order as the omics_list
  elements. If an element is \`NA\`, then rownames are used.

## Value

`SummarizedExperiment` with an additional assay called \`combined\`

## Details

This combines matrices of different omics types together and adds prefix
tags to the ids.

## Examples

``` r
        library(leapR)
        url <- 'https://api.figshare.com/v2/file/download/56536217'

        pdata <- download.file(url,method='libcurl',destfile='protData.rda')
        load('protData.rda')
        p <- file.remove("protData.rda")

        url <- "https://api.figshare.com/v2/file/download/56536214"
        tdata <- download.file(url,method='libcurl',destfile='transData.rda')
        load('transData.rda')
        p <- file.remove("transData.rda")

        url <- 'https://api.figshare.com/v2/file/download/56536211'
        phdata<-download.file(url,method='libcurl',destfile = 'phosData.rda')
        #phosphodata<-read.csv("phdata",check.names=FALSE,row.names=1)
        load('phosData.rda')
        p <- file.remove('phosData.rda')# read in the example protein data


        # merge the three datasets by rows and add prefix tags for
        # different omics types
        multi_omics <- combine_omics(list(pset, tset, phset),
                    list(NA,NA,'hgnc_id'))

```
