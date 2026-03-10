# correlation_comparison_enrichment

\# internal function to calculate enrichment in differences in
correlation \# between two groups \# access through the leapr wrapper

## Usage

``` r
correlation_comparison_enrichment(
  geneset,
  eset,
  assay_name,
  set1,
  set2,
  mapping_column = NA
)
```

## Arguments

- geneset:

  pathway to use for enrichment

- eset:

  SummarizedExperiment with abundance matrix

- assay_name:

  name of assay

- set1:

  first set to use

- set2:

  second set to use

- mapping_column:

  Column to use for id mapping within rowData

## Value

data frame with enrichment results
