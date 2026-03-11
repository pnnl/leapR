# leapR Order Enrichment Tests

## Load libraries needed

``` r
# load the core libraries
library(leapR)
library(gplots)
library(rmarkdown)
# plotting helpers used in this vignette
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(BiocFileCache)
```

## Load our test proteomics dataset

``` r


#currently using the BiocFileCache though i'm not sure it helps
path <- tools::R_user_dir("leapR", which = "cache")
bfc <- BiocFileCache(path, ask = FALSE)

url <- "https://api.figshare.com/v2/file/download/56536217"
pd <- bfcadd(bfc, 'protdat', url)
#> Error while performing HEAD request.
#>    Proceeding without cache information.
load(pd)
#pdata <- download.file(url, method = "libcurl", destfile = "protData.rda")
#  as.matrix()
#load("protData.rda")

#p <- file.remove("protData.rda")

data(shortlist)
data(longlist)

## columns that we want to use for results

cols_to_display <- c("ingroup_n", "outgroup_n", "background_n", 
                     "pvalue", "BH_pvalue")
```

## Compare enrichment methods

We want to evaluate how various samples compare in enrichmnet

``` r

i = 8
data("ncipid")

cor.res <- do.call(rbind,lapply(1:length(shortlist), function (i) {
  protdata.enrichment.ks <- leapR::leapR(
    geneset = ncipid, "enrichment_in_order",
    eset = pset,
    minsize = 5,
    assay_name = "proteomics",
    primary_columns = shortlist[i]
  )
  
  colnames(protdata.enrichment.ks) <- paste('ks',
                                            colnames(protdata.enrichment.ks),
                                            sep = '.')
  
  
  protdata.enrichment.cs <- leapR::leapR(
    geneset = ncipid, "enrichment_in_order",
    method = 'chisq',
    eset = pset,  minsize = 5,
    assay_name = "proteomics",
    primary_columns = shortlist[i]
  )
  colnames(protdata.enrichment.cs) <- paste('chisq',
                                            colnames(protdata.enrichment.cs),
                                            sep = '.')
  
  protdata.enrichment.zt <- leapR::leapR(
    geneset = ncipid, "enrichment_in_order",
    eset = pset,  minsize = 5,
    method = 'ztest',
    assay_name = "proteomics",
    primary_columns = shortlist[i]
  )
  
  colnames(protdata.enrichment.zt) <- paste('ztest',
                                            colnames(protdata.enrichment.zt),
                                            sep = '.')
  
  paths <- rownames(protdata.enrichment.ks)
  allvals <- cbind(protdata.enrichment.cs[paths,], 
                   protdata.enrichment.ks[paths,], 
                   protdata.enrichment.zt[paths,])

  kz <- cor(allvals$ks.BH_pvalue, allvals$ztest.BH_pvalue, use = 'p')
  
  kc <- cor(allvals$ks.BH_pvalue, allvals$chisq.BH_pvalue, use = 'p')
  
  cz <- cor(allvals$chisq.BH_pvalue, allvals$ztest.BH_pvalue, use = 'p')
  
    kz <- cor(allvals$ks.pvalue, allvals$ztest.pvalue, use = 'p')
  
  kc <- cor(allvals$ks.pvalue, allvals$chisq.pvalue, use = 'p')
  
  cz <- cor(allvals$chisq.pvalue, allvals$ztest.pvalue, use = 'p')
  
  return(c(kz = kc, kc=kc, cz = cz))
}))
```

We can compare the corrected p-values

``` r
library(ggplot2)

cor.res |> as.data.frame() |>
  tidyr::pivot_longer(c(1:3),names_to='test_comparison',values_to='p cor') |>
  ggplot(aes(x=test_comparison, y= `p cor`)) + geom_jitter()
```

![](order-enrichment_files/figure-html/compare%20pvals-1.png)
