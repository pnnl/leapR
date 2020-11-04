#leapR

Layered Enrichment Analysis of Pathways in R (leapR) a tool that carries out statistical enrichment analysis on single- or multi-omics data.

## Installation
To install leapR, you can use the `devtools` package as follows:

``` R
install.packages("devtools")
devtools::install_github("biodataganache/leapR")
```

Once you have successfully installed the package you can load the vignette to read examples using the `vignette('leapR')` command.

## Basic Usage

The primary function of the `leapR` package is the `leapR` function itself. This function serves a wrapper to run different styles of enrichment functions on the data.

### Enrichment calls

Here is a list of enrichment arguments that can be called with the `leapR` command.

| Argument                                        | Description |   |
| ---                                             | ----        |   |
| `correlation_comparison_enrichment`             |             |   |
| `correlation_enrichment`                        |             |   |
| `difference_enrichment_in_relationships`        |             |   |
| `enrichment_in_abundance`                       |             |   |
| `enrichment_by_fishers`                         |             |   |
| `enrichment_by_ks`||
| `enrichment_in_relationships`||
|`enrichment_redundancy_matrix`||
|`pairwise_overlap_enrichment`||
|`permute_enrichment_in_groups`||


### Data examples
We included three examples of data.

### Gene pathway examples
We included two different gene pathways.
