# get_pathway_information

get_pathway_information extracts information about a pathway from a
GeneSet object

## Usage

``` r
get_pathway_information(geneset, path, remove.tags = FALSE)
```

## Arguments

- geneset:

  is a GeneSet object for pathway annotation

- path:

  is the name of the gene set pathway to be return

- remove.tags:

  boolean indicating whether to remove tags

## Value

list of pathway information

## Examples

``` r
     library(leapR)

     # load example gene set
     data("ncipid")

     tnfpathway = get_pathway_information(ncipid, "tnfpathway")

```
