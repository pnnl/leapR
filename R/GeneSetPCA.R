
pca_cluster_distinctness <- function(geneset, datatable, groups, labels=NULL, includeall=F, mode="original") {
  if (includeall) {
    # currently doesn't work- plus the data structure for genesets is way clunky
    geneset = add_gene_set(geneset, name="ALLGenes", desc="All genes", genelist=rownames(datatable))
  }
  
  gspca = geneset_pca(geneset, datatable, ncomp=2)
  results = list()
  for (j in geneset$names) {
    results = c(results, list(matrix(nrow=length(groups), ncol=length(groups), dimnames=list(labels, labels))))
  }
  names(results) = geneset$names
  
  for (g in 1:length(groups)) {
    for (h in 1:length(groups)) {
        if (g == h) next
        gg = groups[[g]]
        hh = groups[[h]]
        
        if (mode == "original") {
            # Starting with two groups of treatments/cells gg and hh
            #      this calculates the euclidean distance between all points using the PC1 and PC2 from
            #      PCA. Then it does a t test between the group of distances inside each group with
            #      distances between the groups gg->hh
            # Note: this should work fairly well but may be subject to pathological conditions. For example
            #       if the points in gg are completely surrounded by points from hh.
            group_comp = sapply(rownames(gspca[[1]]), 
                                function (p) {
                                  x=as.matrix(dist(t(rbind(gspca[[1]][p,], gspca[[2]][p,])))); 
                                  l=try(t.test(x[gg, gg][upper.tri(x[gg, gg])], x[gg, hh])$p.value); 
                                  if (is.numeric(l)) return (l); 
                                  return (NA)
                                  })
        }
        else if (mode == "distance") {
          group_comp = sapply(rownames(gspca[[1]]), 
                              function (p) {
                                x=as.matrix(dist(t(rbind(gspca[[1]][p,], gspca[[2]][p,])))); 
                                ggd = x[gg,gg][upper.tri(x[gg,gg])];
                                ggdhh = x[gg,hh];
                                return(mean(unlist(ggdhh), na.rm=T)-mean(unlist(ggd), na.rm=T))
                              })
        }
        
        for (j in geneset$names) {
          results[[j]][g,h] = group_comp[j]
        }
    }
  }
  return (list(gspca=gspca, distmatrices=results))
  
  
}