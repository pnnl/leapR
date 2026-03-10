# enrichment_in_relationships

enrichment_in_relationships function description is a general way to
determine if a pathway is enriched in relationships (interactions,
correlation) between its members \# access through leapr wrapper

## Usage

``` r
enrichment_in_relationships(
  geneset,
  relationships,
  idmap = NA,
  silence_try_errors = TRUE
)
```

## Arguments

- geneset:

  List of pathways in gmt format

- relationships:

  table of relationship information, e.g. correlation

- idmap:

  list of identifiers to use for mapping, the names of the items should
  agree with names of features in matrix

- silence_try_errors:

  boolean to silence errors

## Value

table of enrichment statistics
