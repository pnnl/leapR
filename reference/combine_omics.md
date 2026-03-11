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
        library(BiocFileCache)
        path <- tools::R_user_dir("leapR", which = "cache")
        bfc <- BiocFileCache(path, ask = FALSE)
        
        url <- "https://api.figshare.com/v2/file/download/56536217"
        pc <- bfcadd(bfc, "pdat", fpath = url)
#> 
#> Error while performing HEAD request.
#>    Proceeding without cache information.
        load(pc)
        
        url <- "https://api.figshare.com/v2/file/download/56536214"
        tc <- bfcadd(bfc, "tdat", fpath = url)
#> 
#> Error while performing HEAD request.
#>    Proceeding without cache information.
        load(tc)
        
        url <- "https://api.figshare.com/v2/file/download/56536211"
        phc <- bfcadd(bfc, "phdat", fpath = url)
#> 
#> Error while performing HEAD request.
#>    Proceeding without cache information.
        load(phc)

        # merge the three datasets by rows and add prefix tags for
        # different omics types
        multi_omics <- combine_omics(list(pset, tset, phset),
                    list(NA,NA,'hgnc_id'))

```
