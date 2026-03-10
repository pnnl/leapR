# enrichment_in_abundance

Enrichment in abundance calculates enrichment in pathways by the
difference in abundance of the pathway members.

## Usage

``` r
enrichment_in_abundance(
  geneset,
  eset,
  assay_name,
  mapping_column = NULL,
  abundance_column = NULL,
  fdr = 0,
  matchset = NULL,
  sample_comparison = NULL,
  min_p_threshold = NULL,
  sample_n = NULL,
  silence_try_errors = TRUE
)
```

## Arguments

- geneset:

  Gene set to calculate enrichment

- eset:

  Molecular abundance data in \`SummarizedExperiment\` format

- assay_name:

  Name of assay to compare

- mapping_column:

  Column to use to map identifiers

- abundance_column:

  Columns to use to quantify abundance

- fdr:

  number of times to sample for FDR value

- matchset:

  Name of a set to use for enrichment

- sample_comparison:

  list of samples to use as comparison. if missing background (eset) is
  used

- min_p_threshold:

  Only include p-values lower than this

- sample_n:

  size of sample to use

- silence_try_errors:

  set to true to silence try errors

## Value

data frame of enrichment result
