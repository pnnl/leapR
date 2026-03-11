# calcTTest

calculates a t-test for two distributions of data on a per-gene basis
append results to ExpressionSet with two extra columns: \`pvalue\` and
\`difference\` for each feature

## Usage

``` r
calcTTest(eset, assay_name, group1, group2)
```

## Arguments

- eset:

  SummarizedExperiment

- assay_name:

  name of assay

- group1:

  List of samples comprising group 1

- group2:

  List of samples comprising group 2

## Value

An Expression set with two columns added to the featureData slot:
pvalue, and estimate

## Examples

``` r
        library(leapR)
        library(BiocFileCache)
#> Loading required package: dbplyr
        
        path <- tools::R_user_dir("leapR", which = "cache")
        bfc <- BiocFileCache(path, ask = FALSE)
        
        url <- "https://api.figshare.com/v2/file/download/56536214"
        tc <- bfcadd(bfc, "tdat", fpath = url)
#> 
#> Error while performing HEAD request.
#>    Proceeding without cache information.
        load(tc)
        
        # read in the pathways
        data("ncipid")

        # read in the patient groups
        data("shortlist")
        data("longlist")
        calcTTest(tset, 'transcriptomics', shortlist, longlist)
#> class: SummarizedExperiment 
#> dim: 1999 174 
#> metadata(0):
#> assays(1): transcriptomics
#> rownames(1999): NOC2L ISG15 ... ARL6 MINA
#> rowData names(2): pvalue difference
#> colnames(174): TCGA-13-1484 TCGA-13-1495 ... TCGA-61-1995 TCGA-61-2008
#> colData names(0):
```
