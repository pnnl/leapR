# read_gene_sets

read_gene_sets is a function to import external pathway database files
in .gmt format

## Usage

``` r
read_gene_sets(
  gsfile,
  gene.labels = NA,
  gs.size.threshold.min = 5,
  gs.size.threshold.max = 15000
)
```

## Arguments

- gsfile:

  is a gene set file, for example a .gmt file (gene matrix transposed
  file format)

- gene.labels:

  defaults to NA

- gs.size.threshold.min:

  defaults to 5

- gs.size.threshold.max:

  defaults to 15000

## Value

gene set object

geneset list

## Examples

``` r
gfile <- system.file('extdata','h.all.v2024.1.Hs.symbols.gmt',
          package='leapR')
glist <- read_gene_sets(gfile)
```
