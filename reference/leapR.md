# leapR

leapR is a wrapper function that consolidates multiple enrichment
methods.

## Usage

``` r
leapR(geneset, enrichment_method, eset, assay_name, ...)
```

## Arguments

- geneset:

  is a list of four vectors, gene names, gene descriptions, gene sizes
  and a matrix of genes. It represents .gmt format pathway files.

- enrichment_method:

  is a character string specifying the method of enrichment to be
  performed, one of: "enrichment_comparison", "enrichment_in_order",
  "enrichment_in_sets", "enrichment_in_pathway",
  "correlation_enrichment".

- eset:

  is a \`SummarizedExperiment\` object containing expression data, with
  features as rows and *n* sample/conditions as columns.

- assay_name:

  is the assay to be analyzed within the \`eset\`. Recommended to
  describe the data type (e.g. transcriptomics, proteomics) so that it
  can be integrated in \`combine_omics\`

- ...:

  further arguments

## Value

data frame with results

## Details

Further arguments and enrichment method optional argument information:  

|                                                                                                                                                                                                                                                                                                                                                           |                                                                                                                                                                                                                                                                                                                                                                                                                           |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| id_column                                                                                                                                                                                                                                                                                                                                                 | Is a character string, present in the `rowData` slot, that is used to specify a column for identifiers to map to enrichment libraries. If missing, the rownames of the SummarizedExperiment assay will be used.                                                                                                                                                                                                           |
| primary_columns                                                                                                                                                                                                                                                                                                                                           | Is a character vector composed of column names from `eset` (either in the \`assay\` or in the \`rowData\`), that specifies a set of primary columns to calculate enrichment on. The meaning of this varies according to the enrichment method used - see the descriptions for each method below. This is an optional argument used with 'enrichment_in_order', 'enrichment_in_sets', and 'enrichment_comparison' methods. |
|                                                                                                                                                                                                                                                                                                                                                           | secondary_columns                                                                                                                                                                                                                                                                                                                                                                                                         |
| Is a character vector of column names for comparison, pulled from the \`assay\` of the SummarizedExperiment. This is an optional argument used with 'enrichment_comparison' methods.                                                                                                                                                                      |                                                                                                                                                                                                                                                                                                                                                                                                                           |
| threshold                                                                                                                                                                                                                                                                                                                                                 | Is a numeric value, an optional argument used with 'enrichment_in sets' method which filters out abundance values or p-values (depending on what \`primary_columns\` is used) either above or below it.                                                                                                                                                                                                                   |
|                                                                                                                                                                                                                                                                                                                                                           | greaterthan                                                                                                                                                                                                                                                                                                                                                                                                               |
| Is a logical value that defaults to TRUE, it's used with 'enrichment_in_sets' method. When set to TRUE, genes with \`primary_columns\` value above the `threshold` argument are kept. When set to FALSE genes with \`primary_columns\` value below the `threshold` argument are kept. This is an optional argument used with 'enrichment_in_sets' method. |                                                                                                                                                                                                                                                                                                                                                                                                                           |
| minsize                                                                                                                                                                                                                                                                                                                                                   | Is a numeric value, an optional argument used with 'enrichment_in_sets' and 'enrichment_in_order".                                                                                                                                                                                                                                                                                                                        |
|                                                                                                                                                                                                                                                                                                                                                           | fdr                                                                                                                                                                                                                                                                                                                                                                                                                       |
| A numerical value which specifies how many times to randomly sample genes to calculate an empirical false discovery rate, is an optional argument used with 'enrichment_comparison' method.                                                                                                                                                               |                                                                                                                                                                                                                                                                                                                                                                                                                           |
| min_p_threshold                                                                                                                                                                                                                                                                                                                                           | Is a numeric value, a lower p-value threshold and is an optional argument used with 'enrichment_comparison' method.                                                                                                                                                                                                                                                                                                       |
|                                                                                                                                                                                                                                                                                                                                                           | sample_n                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Is a way to subsample the number of components considered for each calculation randomly. This is an optional argument used with 'enrichment_comparison' method.                                                                                                                                                                                           |                                                                                                                                                                                                                                                                                                                                                                                                                           |

**Enrichment Methods:**  
  
enrichment_comparison  
Compares the distribution of abundances between two sets of conditions
for each pathway using a t test. For each pathway in `geneset` uses a t
test to compare the distribution of abundance values/numbers in `eset`
`primary_columns` with those in `eset` `secondary_columns`. Lower
p-values for pathways indicate that the expression of the pathway is
significantly different between the set of conditions in primary_columns
and the set of conditions in secondary_columns. Optionally, users can
specify `fdr` which will calculate an empirical p-value by randomizing
abundances `fdr` number of times. If the `min_p_threshold` is specified
the method will only return pathways with an adjusted p-value lower than
the specified threshold. If `sample_n` is specified the method will
subsample the pathway members to the specified number of components.  
  
enrichment_in_order  
Calculates enrichment of pathways based on a ranked list using the
Kolmogorov-Smirnov test. For each pathway in `geneset` uses a
Kolmogorov-Smirnov test for rank order to test if the distribution of
ranked abundance values in the `eset` `primary_columns` is significant
relative to a random distribution. Note that currently `primary_columns`
only accepts a single column for this method.  
  
enrichment_in_sets  
Calculates enrichment in pathway membership in a list (e.g. highly
differential proteins) relative to background using Fisher's exact test.
For each pathway in `geneset` uses a Fisher's exact test over- or under-
representation of a list of components specified. If `targets` are
specified this must be a vector of identifiers to serve as the target
list for comparison. If `eset` and `primary_columns` are specified then
`threshold` specifies a threshold value for determining the target list
of components to test. Specifying `greaterthan` to be False will result
in components with values lower than the specified `threshold`. If
`eset` is a data frame or matrix, the background used for calculation
will be taken as the rownames of `eset`  
  
enrichment_in_pathway  
Compares the distribution of abundances in a pathway with the background
distribution of abundances using a t test. For each pathway in `geneset`
calculates the significance of the difference between the abundances
from pathway members versus abundance of non-pathway members in the set
of conditions specified by `primary_columns`. Optionally, users can
specify `fdr` which will calculate an empirical p-value by randomizing
abundances `fdr` number of times. If the `min_p_threshold` is specified
the method will only return pathways with an adjusted p-value lower than
the specified threshold. If `sample_n` is specified the method will
subsample the pathway members to the specified number of components.  
  
correlation_enrichment  
Calculates the enrichment of a pathway based on correlation between
pathway members across conditions versus correlation between members not
in the pathway. For each pathway in `geneset` calculates the pairwise
correlation between all pathway members and non-pathway members across
the specified `primary_columns` conditions in `eset`. Note that for
large matrices this can take a long time. A p-value is calculated based
on comparing the correlation within the members of a pathway with the
correlation values between members of the pathway and non-members of the
pathway.  

## Examples

``` r
        library(leapR)

 # read in the example abundance data
 # read in the example transcriptomic data
 tdata <- download.file("https://api.figshare.com/v2/file/download/56536214",
      method='libcurl',destfile='transData.rda')
 load('transData.rda')
 p <- file.remove("transData.rda")

 # read in the pathways
 data("ncipid")

 # read in the patient groups
 data("shortlist")
 data("longlist")

 # use enrichment_comparison to calculate enrichment in one set of
 # conditions (shortlist) and another (longlist)
 short_v_long = leapR(geneset=ncipid, assay_name='transcriptomics',
              enrichment_method='enrichment_comparison',
              eset=tset, primary_columns=shortlist,
               secondary_columns=longlist)

 # use enrichment_in_sets to calculate the most enriched pathways
 # from the highest abundance proteins
 #     from one condition
 onept_sets = leapR(geneset=ncipid, assay_name='transcriptomics',
               enrichment_method='enrichment_in_sets',
               eset=tset, primary_columns="TCGA-13-1484", threshold=1.5)

 # use enrichment_in_order to calculate the most enriched pathways from the
 # same condition
 # Note: that this uses the entire set of abundance values and their order -
 # whereas the previous example uses a hard threshold to get a short list of
 # most abundant proteins and calculates enrichment based on set overlap.
 # The results are likely to be similar - but with some notable differences.
 onept_order = leapR(geneset=ncipid, assay_name='transcriptomics',
               enrichment_method='enrichment_in_order',
               eset=tset, primary_columns="TCGA-13-1484")
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties
#> Warning: p-value will be approximate in the presence of ties

 # use enrichment_in_pathway to calculate the most enriched pathways in a
 # set of conditions based on abundance in the pathway members versus
 # abundance in non-pathway members
 short_pathways = leapR(geneset=ncipid, assay_name='transcriptomics',
               enrichment_method='enrichment_in_pathway',
               eset=tset, primary_columns=shortlist)

 # use correlation_enrichment to calculate the most enriched pathways in
 # correlation across the shortlist conditions
 short_correlation_pathways = leapR(geneset=ncipid,
                assay_name='transcriptomics',
                enrichment_method='correlation_enrichment',
                eset=tset, primary_columns=shortlist)

```
