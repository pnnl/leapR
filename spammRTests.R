###quick spammr Tests
library(spammR)
library(leapR)

data(smallPancData)
data(pancMeta)
data(protMeta)
pooledPanc <- dplyr::bind_cols(smallPancData)
panc.spe <- convert_to_spe(pooledPanc,pancMeta,protMeta,feature_meta_colname='pancProts',samples_common_identifier='')
diffex.spe <- calc_spatial_diff_ex(panc.spe,category_col='IsletOrNot')
res<-as(diffex.spe,'ExpressionSet')
data(krbpaths)

ora.res <- enrich_ora(res, geneset=krbpaths,geneset_name='krbpaths',id_column='PrimaryGeneName',
                      feature_column="NonIslet_vs_Islet.adj.P.Val.limma",threshold=0.1,greaterthan=FALSE)

ora.res <- leapR(eset=res, geneset=krbpaths,enrichment_method='enrichment_in_sets',id_column='PrimaryGeneName',
                 primary_column="NonIslet_vs_Islet.adj.P.Val.limma",threshold=0.05,greaterthan=FALSE)