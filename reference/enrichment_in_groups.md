# enrichment_in_groups

Calculate the enrichment in pathways using Fisher's exact or
Kolmogorov-Smirnov test, using either the abundance column to identify
feature or the targets list. access through leapr wrapper

## Usage

``` r
enrichment_in_groups(
  geneset,
  targets = c(),
  background = NULL,
  assay_name = NULL,
  method = "fishers",
  minsize = 5,
  mapping_column = NULL,
  log_transformed = FALSE,
  abundance_column = NULL,
  randomize = FALSE,
  silence_try_errors = TRUE
)
```

## Arguments

- geneset:

  geneset to use for enrichment

- targets:

  targets to use for enrichment

- background:

  \`SummarizedExperiment\` describing background to use

- assay_name:

  is the name of the assay to use from the background

- method:

  method to use for statistical test, options are 'fishers', 'ks',
  'ztest', or 'chisq'. Remember that KS test assumes normality, so it
  would be good to log your data before calling. NOTE: if you do not
  call \`suppressWarnings\` then the KS test will warn you about ties.

- minsize:

  minimum size of set

- mapping_column:

  column name of mapping identifiers

- log_transformed:

  Set to TRUE if data is already log-transformed

- abundance_column:

  columns mapping abundance, either in the \`assay\` matrix or
  \`rowData\`

- randomize:

  true/false whether to randomize

- silence_try_errors:

  true/false to silence errors

## Value

data frame with enrichment results
