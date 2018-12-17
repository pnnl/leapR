# Scripts to do functional enrichment in groups
#
#

#library(GOstats)
#library(GO.db)
#library(hgu95av2.db)
#library(AnnotationDbi)

map_gs_to_entrezgeneid <- function(genes, gene_map) {
	# complicated and probably overblown way to translate a list of gene symbols to entrezgene ids
	return(gene_map[genes[which(genes %in% rownames(gene_map))],]$GeneID)
}

map_gs_to_entrezgene_affy <- function(genes, mode="symbol", gmap=hgu95av2ALIAS2PROBE, emap=hgu95av2ENTREZID) {
	if (mode == "symbol") {
		probeids = unlist(mget(genes, gmap, ifnotfound=NA))
		probeids = probeids[sapply(probeids, function(x) !is.na(x))]
		if (length(probeids) == 0) {
			return(c())	
		}
	}
	else if (mode == "probe") {
		probeids = genes	
	}
	universe = unlist(mget(probeids, envir=emap, ifnotfound=NA))
	universe = universe[sapply(universe, function(x) !is.na(x))]
	return(universe)
}

hypergo_by_list <- function(genelist, gene_map=NA, organism="human", annotation='hgu95av2.db', 
								annotation_db=hgu95av2ENTREZID, alias_map=hgu95av2ALIAS2PROBE, 
								map_mode="symbol", minsize=10, ontology="BP", pvalueCutoff=0.001, conditional=TRUE, 
								big_universe=TRUE) {
	# just a wrapper around the hypergo by cluster to
	#      run a single list of genes
	dumclust = rep(1, length(genelist))
	names(dumclust) = genelist
	return(hypergo_by_cluster(dumclust, gene_map, organism, annotation, annotation_db, alias_map,
			map_mode, minsize=0, ontology, pvalueCutoff, conditional, big_universe))
}

hypergo_by_cluster <- function(cutclusters, gene_map=NA, organism="human", annotation='hgu95av2.db', 
								annotation_db=hgu95av2ENTREZID, alias_map=hgu95av2ALIAS2PROBE, 
								map_mode="symbol", minsize=10, ontology="BP", pvalueCutoff=0.001, conditional=TRUE, 
								big_universe=TRUE) {
	# take a treecut from hclust and perform hyperGOtest on each
	#      cluster in turn- uses a gene_map to map gene symbols to entrezgeneids
	
	# shortcut
	if (organism == "mouse") {
		# if this fails: 
		# source("http://bioconductor.org/biocLite.R")
		# biocLite("mouse4302.db")
		# to download the annotation file
		
		# Affymetrix chip mapping
		library(mouse4302.db)
		annotation = "mouse4302.db"
		annotation_db = mouse4302ENTREZID
		alias_map = mouse4302ALIAS2PROBE
	} else if (organism == "mouse_moe430a") {
		# also Affy chip mapping
		library(moe430a.db)
		annotation = "moe430a.db"
		annotation_db = moe430aENTREZID
		alias_map = moe430aALIAS2PROBE
	} else if (organism == "mouse_agilent") {
		# Agilent 4x44k mapping
		library(mgug4122a.db)
		annotation = "mgug4122a.db"
		annotation_db = mgug4122aENTREZID
		alias_map = mgug4122aALIAS2PROBE
		
	} else if (organism == "human_agilent") {
		library()
		biocLite("")
		
	}
	
	# we can add more annotations
	
	# map gene symbols to entrezgeneids
	#universe = map_gs_to_entrezgeneid(names(cutclusters), gene_map)
	
	if (big_universe == TRUE) {
		# get them all
		universe = mget(keys(annotation_db), envir=annotation_db)	
	} else if (big_universe == FALSE){	
		universe = map_gs_to_entrezgene_affy(names(cutclusters), emap=annotation_db, gmap=alias_map, mode=map_mode)
	} else {
		universe = map_gs_to_entrezgene_affy(big_universe, emap=annotation_db, gmap=alias_map, mode=map_mode)
	}
			
	params = new("GOHyperGParams", geneIds=universe, universeGeneIds=universe,
					annotation=annotation, ontology=ontology, pvalueCutoff=pvalueCutoff,
					conditional=conditional, testDirection="over")
	
	results = list()
	cnames = c()
	
	for (i in 1:max(cutclusters)) {
		cat("Functional enrichment for cluster ", i, "\n")
		paramsclust = params
		
		clustids = map_gs_to_entrezgene_affy(names(which(cutclusters==i)), emap=annotation_db, gmap=alias_map, mode=map_mode)
		#cat(str(clustids), "\n")
		
		clustids[which(clustids %in% universe)]
		if (length(clustids) >= minsize) {
			cat("Calculating hyperGOtest for ", length(clustids), "genes\n")
			
			geneIds(paramsclust) = clustids
			
			hgOver = hyperGTest(paramsclust)
			#testDirection(paramsclust) = "under"
			#hgUnder = hyperGTest(paramsclust)
			
			cat("\tfinished.\n")
			result = c()
			result$cluster = i
			result$over = hgOver
			cnames = c(cnames, paste("cluster_", i, sep=""))
			#result$under = hgUnder
			results = c(results, list(result))
		}
	}
	names(results) = cnames
	return(results)
}

summarize_hypergo_clusters <- function(hypergo_record, pvalue_threshold=0.05) {
	# takes a list of clusters as output above and creates a summary measure
	# in this case a simple count of genes that are in groups below the pvalue_threshold
	# Not sure if this is really fair- it likely overrepresents particular genes
	x = 0
	for (i in 1:length(hypergo_record)) {
		#nm = paste("cluster_", hypergo[[i]]$cluster)
		gc = geneCounts(hypergo_record[[i]]$over)[which(pvalues(hypergo_record[[i]]$over) < pvalue_threshold)]
		
		gids = geneIdsByCategory(hypergo_record[[i]]$over)
		
		genelist = sapply(names(gc), function (g) gids[g])
		genelist = unique(genelist)
		
		y = length(genelist)
		x = x + y
	}	
	return(x)
}

randomize_cluster_membership <- function(ntimes=10, cutclusters=NA, gene_map=NA, organism="human", annotation='hgu95av2.db', 
											annotation_db=hgu95av2ENTREZID, alias_map=hgu95av2ALIAS2PROBE, 
											map_mode="symbol", minsize=10, ontology="BP", pvalueCutoff=0.001, 
											conditional=TRUE, big_universe=TRUE) {
	
    # takes same input as hypergo_by_cluster above. This will just run ntimes random reassortments of
    #    cluster membership, do the functional enrichment, and then do the summary count metric from above.
	results = list()
	for (i in 1:ntimes) {
		randcut = sample(cutclusters, size=length(cutclusters))
		names(randcut) = names(cutclusters)
		
		result = hypergo_by_cluster(cutclusters=cutclusters, gene_map=gene_map, organism=organism, annotation=annotation, 
											annotation_db=annotation_db, alias_map=alias_map, 
											map_mode=map_mode, minsize=minsize, ontology=ontology, pvalueCutoff=pvalueCutoff, 
											conditional=conditional, big_universe=big_universe)
		
		count = summarize_hypergo_clusters(result)
		results = c(results, list(result, count))
		
	}
	return(results)
}

write_function_table <- function(hypergo_record, outfile, cluster_names=NA, pvalue=0.01) {
	# takes a hypergo cluster record produced by
	#       hypergo_by_cluster and writes it to a tdl file
	
	for (i in 1:length(hypergo_record)) {
		if (! is.na(cluster_names)) {
			if (cluster_names == 1) {
				nm = cluster_names[hypergo_record[[i]]$cluster]
			}
			else {
				nm = cluster_names[i]	
			}
		}
		else {
			nm = paste("cluster_", hypergo_record[[i]]$cluster)
		}
		
		cat("cluster\t", nm, "\n", file=outfile, append=TRUE, sep="")
		write.table(summary(hypergo_record[[i]]$over, pvalue=pvalue), file=outfile, sep="\t", quote=FALSE, append=TRUE)
	}
}

enrichment_by_cluster <- function(clusters, annotation, mode="fishers", valuecol=1) {
	# for each cluster in the hclust-type clusters
	# assess the statistical enrichment in the annotation list
	results = c()
	
	for (i in 1:max(clusters)) {
		incluster = names(clusters)[which(clusters == i)]
			
		if (mode == "fishers") {	
			ft = enrichment_by_fishers(incluster, names(clusters), rownames(annotation))
		} else if (mode == "ttest") {
			ft = enrichment_by_ttest(incluster, names(clusters), rownames(annotation), values=annotation, valuecol=valuecol)
		}
		results = c(results, list(ft))
	}
	return(results)
}

enrichment_by_condition <- function(clusters, ratios, condrange=NA) {
	# calculate statistical enrichment of a cluster for a range of conditions
	library(stats)
	
	if (is.na(condrange)) {
		condrange = colnames(ratios)	
	}
	
	sigs = matrix(nrow=max(clusters), ncol=length(condrange))
	means = matrix(nrow=max(clusters), ncol=length(condrange))
	back = matrix(nrow=max(clusters), ncol=length(condrange))
	dev = matrix(nrow=max(clusters), ncol=length(condrange))
	zscores = matrix(nrow=max(clusters), ncol=length(condrange))
	colnames(sigs) = colnames(ratios)
	rownames(sigs) = sapply(1:max(clusters), function(i) paste("cluster_", i, sep=""))
	dimnames(means) = dimnames(sigs)
	dimnames(back) = dimnames(sigs)
	dimnames(dev) = dimnames(sigs)
	dimnames(zscores) = dimnames(sigs) 
	
	for (i in 1:max(clusters)) {
		rows = names(clusters[which(clusters==i)])
    rows = rows[which(rows %in% rownames(ratios))]
		in_cluster = ratios[rows,]
		out_cluster = ratios[! rows %in% names(ratios),]
		if (class(in_cluster) == "matrix") {
			for (j in condrange) {
        if (length(in_cluster) > 2) {
				s = t.test(in_cluster[,j], out_cluster[,j])
				sigs[i,j] = s$p.value
				m = mean(in_cluster[,j])
				means[i,j] = m
				b = mean(out_cluster[,j])
				back[i,j] = b
				d = sd(in_cluster[,j])
				dev[i,j] = d
				z = (m-b)/d
				zscores[i,j] = z
        }
			}
		}
	}
	return(list(mean=means, pvalue=sigs, zscore=zscores, sd=dev, back=back))
}

# Functional enrichment code
#    input: a matrix or data frame of proteins (rows) with different numerical values (columns)
#           a data frame that has a list of proteins in a functional group
# Example input file:
#        Degree  BonPwr  2Step   ARD     Eigenvector     Between
# CBX5    0.401658773     51.80116272     0.951421797     0.692733169     0.086658746     0.014533159
# PSAP    0.359004736     40.90088654     0.946682453     0.670616388     0.068314366     0.009966435
# ...
#
# Example annotation file:
# ACE2    dyerlist        PathoTarget     1.0
# ACHE    dyerlist        PathoTarget     1.0
# ACSS2   dyerlist        PathoTarget     1.0
# ... (this only requires the list of protein names in the first column- the rest is dispensable)
#
# 1. Read in the data frame of proteins and values:
#			input = read.table([filename], sep="\t", header=TRUE, row.names=1)
# 2. Read in the data frame of annotations:
# 			annotations = read.table([filename], sep="\t", row.names=1)
# 3. Do the enrichment for the top 20% of a particular column:
#			results = enrichment_in_top(input, annotations, sort.column="Between", percentage=0.2)
# 4. Do the enrichment for the top 20% of all columns:
#  			results = enrichment_columns(input, annotations, percentage=0.2)
# The results contain a list of:
#      1. $fisher : Fisher's exact test
#      2. $mat : Matrix of counts in groups. Final column (percentage) gives the percentage enrichment of the functional group
#                in the group (to 20% of values), and the background (rest of proteins in the input list)
# 	   3. $fold : fold enrichment of group to background calculated by dividing the in-group-function percentage by the
#				  in-background-function percentage. This provides a summary of results.
#
enrichment_by_fishers <- function(group, background, annotation) {
	# calculate the fishers exact on a group of things versus a background
	#     for those things with a particular annotation (on another list)
	non_group = background[!background %in% group]
	
	group_annot = sum(group %in% annotation)
	non_group_annot = sum(non_group %in% annotation)
	
	group_nonannot = sum(!group %in% annotation)
	non_group_nonannot = sum(!non_group %in% annotation)
	
	test = matrix(c(group_annot, non_group_annot, group_nonannot, non_group_nonannot), nr=2,
					dimnames=list(c("Group", "NonGroup"), c("Annotated", "NonAnnotated")))
	
	per = c(test[1,1]/(test[1,1]+test[1,2]), test[2,1]/(test[2,1]+test[2,2]))
					
	ft = fisher.test(test)
	test = cbind(test, per)
	#print(test)
	#cat(sum(test), length(background), "\n")
	
	fold = per[1]/per[2]
	
	return(list(fisher=ft, mat=test, foldx=fold))
}

enrichment_by_ttest <- function(group, background, annotation, values=NA, valuecol=1) {
	# calculate the fishers exact on a group of things versus a background
	#     for those things with a particular annotation (on another list)
	non_group = background[!background %in% group]
	
	group_list = values[which(annotation %in% group),valuecol]
	non_group_list = values[which(annotation %in% non_group),valuecol]
						
	tt = t.test(group_list, non_group_list)
	foldx = tt$estimate[1]/tt$estimate[2]

	return(list(ttest=tt, mat=list(group_list, non_group_list), foldx=foldx))
}

enrichment_in_top <- function(data, annotation, sort.column="V1",
							percentage=0.2, random=FALSE, decreasing=TRUE, 
							mode="fishers", valuecol=1) {
	# calculate the enrichment of the top portion of a list

	if (random) {
		# randomly order the matrix
		ordered_data = data[sample(rownames(data), nrow(data)),]
	}
	else {
		ordered_data = data[order(data[,sort.column], decreasing=decreasing),]
	}
	
	top_part = rownames(ordered_data[1:(length(rownames(ordered_data))*percentage),])
		
	if (mode == "fishers") {	
		ft = enrichment_by_fishers(top_part, rownames(ordered_data), rownames(annotation))
	} else if (mode == "ttest") {
		ft = enrichment_by_ttest(top_part, rownames(ordered_data), rownames(annotation), values=annotation, valuecol=valuecol)	
	}
	
	return(c(ft, list(top=top_part)))
}

format_enrichment <- function(data) {
	# data is a list of enrichment records (from above)
	group = sapply(data, function (i) i$mat[1,"per"])
	background = sapply(data, function (i) i$mat[2,"per"])
	foldx = sapply(data, function (i) i$foldx)
	pvalue = sapply(data, function (i) i$fisher$p.value)	
	ret = rbind(group, background, foldx, pvalue)
	
	#ret = matrix(group, background, foldx, pvalue, nrow=4, byrow=TRUE)
	
	return(ret)
}

enrichment_columns <- function(data, annotation, percentage=0.2, mode="fishers", valuecol=1, decreasing=T) {
	ret = matrix(nrow=5, ncol=ncol(data))
	colnames(ret) = colnames(data)
	rownames(ret) = c("Group", "NonGroup", "FoldX", "p-value", "N")
	for (c in colnames(data)) {
		ft = enrichment_in_top(data, annotation, sort.column=c, percentage=percentage, 
								mode=mode, valuecol=valuecol, decreasing=decreasing)
		if (mode == "fishers") {
			ret["Group",c] = ft$mat[1,"per"]
			ret["NonGroup",c] = ft$mat[2,"per"]
			ret["FoldX",c] = ft$foldx
			ret["p-value",c] = ft$fisher$p.value
			ret["N",c] = nrow(data)
		} else if (mode == "ttest") {
			ret["Group",c] = ft$ttest$estimate[[1]]
			ret["NonGroup",c] = ft$ttest$estimate[[2]]
			ret["FoldX", c] = ft$foldx
			ret["p-value",c] = ft$ttest$p.value
			ret["N",c] = nrow(data)
		}
	}
	return(ret)
}

random_simple_go_path <- function(golist, path=c()) {
  # for demonstration purposes retrieve a single GO path to root
  #     and discard everything else
  
  # randomly choose a parent
  #cat("0:", golist, "\n")
  golist = golist[which(golist %in% keys(GOBPPARENTS))]
  #cat("1:", golist, "\n")
  if (length(golist) == 0) {
    return(path)
  }
  
  if (length(golist) == 1) thisgo = golist
  else thisgo = sample(as.vector(golist), 1)
  path = c(Term(GOTERM[thisgo]), path)
  #cat(">>>", thisgo, "\n")
  
  # get parents of this parent
  golist = as.vector(as.list(GOBPPARENTS[thisgo])[[1]])
  #cat("2:", golist, "\n")
  
  path = random_simple_go_path(golist, path)
}

random_simple_go_paths <- function(go_table) {
  # operates on a table of ids and go identifers
  # and returns a vector with unique rownames where
  #     each entry is a list containing the longest GO path
  
  ids = unique(go_table[,1])
  
  results = c()
  for (id in ids) {
    res = random_simple_go_path(go_table[which(go_table[,1]==id),2])
    results = c(results, list(list(id=id, results=res)))
  }
  return(results)
}

get_deepest_go_path <- function(golist) {
  # operating on a list of GOIDs (as from the biomaRt getBM table) 
  #  returns the longest ancestor path
  # Actually returns an unstructured set of ancestors
  
  # ignore ids not in GOBP
  howlongs = sapply(golist, function (id) try(length(as.list(GOBPANCESTOR[id])[[1]]),silent=T))
  
  winner = golist[which.max(howlongs)]
  goids = sapply(as.list(GOBPANCESTOR[winner])[[1]], function (g) Term(GOTERM[g]))
  
  return(c(goids, Term(GOTERM[winner])))
}

get_deepest_go_paths <- function(go_table) {
  # operates on a table of ids and go identifers
  # and returns a vector with unique rownames where
  #     each entry is a list containing the longest GO path
  
  ids = unique(go_table[,1])
  
  results = c()
  for (id in ids) {
    res = get_deepest_go_path(go_table[which(go_table[,1]==id),2])
    results = c(results, list(list(id=id, results=res)))
  }
  return(results)
}

