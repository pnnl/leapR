# correlation_enrichment

\# calculate enrichment in correlation between pathway members \# access
through leapr wrapper

## Usage

``` r
correlation_enrichment(geneset, eset, assay_name, mapping_column = NA)
```

## Arguments

- geneset:

  Geneset list

- eset:

  a SummarizedExperiment object

- assay_name:

  name of assay

- mapping_column:

  Column to use to map identifiers, if not rownames

## Value

list of enrichment statistic table and correlation matrix
