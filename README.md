# leapR

Layered Enrichment Analysis of Pathways in R (leapR) a tool that carries out statistical enrichment analysis on single- or multi-omics data.

## Installation
To install leapR, you can use the `devtools` package as follows:

``` R
install.packages("devtools")
devtools::install_github("PNNL-CompBio/leapR",build_vignette=TRUE)
```

Once you have successfully installed the package you can load the vignette to read examples using the `vignette('leapR')` command.

## Basic Usage

The primary function of the `leapR` package is the `leapR` function itself. This function serves a wrapper to run different styles of enrichment functions on the data. The package contains other functions to support pathway information and multi-omics datasets.

### Enrichment calls

Here is a list of enrichment arguments that can be called with the `leapR` command.

| Argument                                        | Description |   
| ---                                             | ----        |   
| `enrichment_in_sets`                           | Calculates enrichment in pathway membership in a list (e.g. highly differential proteins) relative to background using Fisher's exact test. |   
| `enrichment_in_order`                         | Calculates enrichment of pathways based on a ranked list using the Kologmorov-Smirnov test |   
| `enrichment_comparison`                         | Compares the distribution of abundances between two sets of conditions for each pathway using a t test  |
| `enrichment_in_pathways`                        | Compares the distribution of abundances in a pathway with the background distribution of abundances using a t test |   
|`correlation_enrichment`              | Calculates the enrichment of a pathway based on correlation between pathway members across conditions versus correlation between members not in the pathway             |   
| `enrichment_in_relationships`| Calculates the enrichment of a pathway in specified interactions relative to non-pathway members |

### Data examples
We included examples of including proteomics data and transcriptomics data from 169 high-grade serous ovarian cancer (HGSOC) tumors previously studied and lists of the short- and long- surviving patients from that cohort.

### Gene pathway examples
We included two different gene pathways. An NCI pathway database (Pathway Information Database; PID) of signaling pathways and the MSIGDB set of gene collections from various sources.
