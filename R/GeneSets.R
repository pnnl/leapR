# Read GSEA gene set file
library(stats)
library(readr)

read_gene_sets <- function(gsfile, gene.labels=NA, gs.size.threshold.min=5, gs.size.threshold.max=15000) {
  # Read input gene set database
      temp <- read_lines(gsfile)
      
      max.Ng <- length(temp)
      temp.size.G <- vector(length = max.Ng, mode = "numeric") 
      for (i in 1:max.Ng) {
          temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
      }

      max.size.G <- max(temp.size.G)      
      gs <- matrix(rep(NA, max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
      temp.names <- vector(length = max.Ng, mode = "character")
      temp.desc <- vector(length = max.Ng, mode = "character")
      gs.count <- 1
      for (i in 1:max.Ng) {
          gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
          gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
          gene.set.name <- gs.line[1] 
          gene.set.desc <- gs.line[2]
          #cat(gene.set.name, "\n")
          #cat(gs.line)
          gene.set.tags <- vector(length = gene.set.size, mode = "character")
          for (j in 1:gene.set.size) {
              gene.set.tags[j] <- gs.line[j + 2]
          } 
          if (is.na(gene.labels)) {
          	  # we just want to load it all (unless we want to filter here?)
          	existing.set = rep(T, length(gene.set.tags))
          } else {
         	existing.set <- is.element(gene.set.tags, gene.labels)
          }
          
          set.size <- length(existing.set[existing.set == T])
          if ((set.size < gs.size.threshold.min) || (set.size > gs.size.threshold.max)) next
          temp.size.G[gs.count] <- set.size
          gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
          temp.names[gs.count] <- gene.set.name
          temp.desc[gs.count] <- gene.set.desc
          gs.count <- gs.count + 1
      } 
      Ng <- gs.count-1
      gs.names <- vector(length = Ng, mode = "character")
      gs.desc <- vector(length = Ng, mode = "character")
      size.G <- vector(length = Ng, mode = "numeric") 
      gs.names <- temp.names[1:Ng]
      gs.desc <- temp.desc[1:Ng] 
      size.G <- temp.size.G[1:Ng]
  	  return(list(names=gs.names, desc=gs.desc, sizes=size.G, matrix=gs))
  }

add_gene_set <- function(geneset, name=NA, desc=NA, genelist=NA) {
  geneset$names = c(geneset$names, name)
  geneset$desc = c(geneset$desc, desc)
  geneset$size = c(geneset$size, length(genelist))
  if (ncol(geneset$matrix)<length(genelist)) {
    geneset$matrix = cbind(geneset$matrix, rep("null", nrow(geneset$matrix)*(length(genelist)-ncol(geneset$matrix))))
  }
  add = genelist
  if ((ncol(geneset$matrix)-length(genelist))>0) add = c(genelist, rep("null", (ncol(geneset$matrix)-length(genelist))))
  
  geneset$matrix = rbind(geneset$matrix, add)
  
  return (geneset)
}

make_gene_set <- function(name, desc, genelist) {
  return(list(names=c(name,"dum"), 
              desc=c(desc,"dum"), 
              size=c(length(genelist),1), 
              matrix=t(data.frame(genelist, rep("null", length(genelist)), ncol=length(genelist)))))
}

extract_gene_set <- function(geneset, setlist) {
  x = which(geneset$names %in% setlist)
  names = geneset$names[x]
  notfound = length(setlist[which(!setlist %in% geneset$names)])
  cat(notfound, "sets not found in geneset\n")
  
  desc = geneset$desc[x]
  size = geneset$sizes[x]
  mat = geneset$matrix[x,]
  
  return(list(names=names, desc=desc, sizes=size, matrix=mat))
}

overlap_gene_sets <- function(geneset, set1, set2) {
  a = get_pathway_information(geneset, set1)$geneset
  b = get_pathway_information(geneset, set2)$geneset

  overlap = a[a %in% b]
  uniquea = a[!a %in% b]
  uniqueb = b[!b %in% a]
  
  return(list(overlap=overlap, set1_unique=uniquea, set2_unique=uniqueb))
  
}

combine_gene_sets <- function(genesetA, genesetB, tagA=NULL, tagB=NULL) {
  # TODO
  
}

get_pathway_information <- function(geneset, path, remove.tags=F) {
  i = which(geneset$names == path)
  thissize = geneset$size[i]
  thisdesc = geneset$desc[i]
  
  # default is to look through all members for effects on the pathway
  grouplist = geneset$matrix[i,1:thissize]
  grouplist = unique(sort(grouplist[which(grouplist != "")]))
  
  if (remove.tags) {
    grouplist = unique(sapply(grouplist, function (n) strsplit(n, "_")[[1]][2]))
  }
  thissize = length(grouplist)
  
  return(list(name=path, size=thissize, description=thisdesc, geneset=grouplist))
}

permute_enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5, 
                                         mapping_column=NULL, abundance_column=NULL, ntimes=100) {
  # Permute Fisher's exact test and provide an FDR for the real thing
  
  real = enrichment_in_groups(geneset, targets=targets, background=background, minsize=minsize,
                              mapping_column=mapping_column, abundance_column=abundance_column,
                              method=method)
  
  # all we need to do is keep track of the number of times we get a better
  #    p-value out of the permutations
  results = rep(0, nrow(real))
  names(results) = rownames(real)
  
  for (i in 1:ntimes) {
    # permute labels
    unreal = enrichment_in_groups(geneset, targets=targets, background=background, minsize=minsize, 
                                  mapping_column=mapping_column, abundance_column=abundance_column, 
                                  method=method, randomize=T)
    
    results[!is.na(unreal[,"pvalue"])] = results[!is.na(unreal[,"pvalue"])] + 
      (unreal[which(!is.na(unreal[,"pvalue"])),"pvalue"]<real[which(!is.na(unreal[,"pvalue"])),"pvalue"])
  }
  # return the percentage of times that each functional group
  #        permutation returns a p-value better than the real
  #        functional group
  return(results/ntimes)
}

enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5,
                                 mapping_column=NULL, abundance_column=NULL, randomize=F) {
  
  resultp = c()
  resultf = c()
  results = data.frame(row.names = geneset$names, 
                       in_path=rep(NA_real_, length(geneset$names)), in_path_names=rep(NA_character_, length(geneset$names)), out_path=rep(NA_real_, length(geneset$names)),
                       in_back=rep(NA_real_, length(geneset$names)), out_back=rep(NA_real_, length(geneset$names)), foldx=rep(NA_real_, length(geneset$names)),
                       pvalue=rep(NA_real_, length(geneset$names)), Adjusted_pvalue=rep(NA_real_, length(geneset$names)), Signed_AdjP=rep(NA_real_, length(geneset$names)), 
                       stringsAsFactors = F)
  
  if (method == "ks") colnames(results)[c(3,5)] = c("MeanPath", "Zscore")
  
  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    grouplist = geneset$matrix[i,1:thissize]
    if (randomize) {
      # choose a random set of genes as this grouplist
      # A disadvantage is that we resample for each functional group rather than
      #   running one set of analyses on a fully scrambled set of functions.
      #   I don't think this should be a huge problem though.
      grouplist = sample(unlist(geneset$matrix), length(grouplist))
    }
    in_back = length(background)
    
    if (method == "fishers") {
      enr = enrichment_by_fishers(targets, background, grouplist)
      p = enr$fisher$p.value
      f = enr$foldx
      mat = enr$mat
      names = enr$in_path_names
      
      results[thisname, ] = list(mat[1,1], names, mat[1,2], mat[2,1], mat[2,2], f, p, NA, NA)
    }
    else if (method == "ks") {    #Kolmogorov-Smirnov test
      # in this case "background" must be the continuous variable from which grouplist can be drawn
      backlist = background
      
      if (is.null(mapping_column)) {
        in_group = background[grouplist[which(grouplist %in% rownames(background))],abundance_column]
        in_group_name = paste(intersect(grouplist, rownames(background)), collapse = ", ")
        backlist = background[,abundance_column]
      }
      else {
        # mapping_column adds the ability to use phospho-type data where the gene name (non-unique) is in the
        #       first column and the rownames are peptide ids
        # unfortunately this means that "background" has to be the whole matrix and abundance_column
        #       has to be specified, which is a bit ugly
        in_group = background[which(background[,mapping_column] %in% grouplist),abundance_column]
        in_group_name = paste(intersect(background[,mapping_column], grouplist), collapse = ", ")
        backlist = background[,abundance_column]
      }
      
      in_path = length(in_group)
      
      
      if (in_path > minsize) {
        in_back = length(backlist)
        
        enr = try(ks.test(in_group, backlist))
        if (class(enr) == "try-error") {
          enr = NA
          p.value = NA
        }
        else {
          p.value = enr$p.value
        }
        
        # this expression of foldx might be subject to some weird pathological conditions
        # e.g. one sample has a background that is always negative, another that's positive
        # may pertain to zscore too (although not sure it should)
        #foldx = mean(in_group, na.rm=T)/mean(background, na.rm=T)
        
        # rank from largest to smallest
        if (is.null(mapping_column)) in_rank = rank(backlist)[grouplist[which(grouplist %in% names(background))]]
        else in_rank = rank(backlist)[which(background[,mapping_column] %in% grouplist)]
        
        foldx = mean(in_rank, na.rm=T)/length(backlist)
        
        zscore = (mean(in_group, na.rm=T)-mean(backlist, na.rm=T))/sd(in_group, na.rm=T)
        
        #padj = enr$p.value*length(geneset$names)
        # c("in_path", "MeanPath", "in_back", "Zscore", "foldx", "pvalue", "Adjusted_pvalue")
        results[thisname, ] = list(in_path, in_group_name, mean(in_group, na.rm=T), in_back, zscore, foldx, p.value, NA, NA)
      }
    }
  }
  results[,"Adjusted_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  results[,"Signed_AdjP"] = results[,"Adjusted_pvalue"]*sign(results[,"in_path"])
  return(results)
}

enrichment_in_relationships <- function(geneset, relationships, idmap=NA, tag=NA, mode="original") {
  # for each category in geneset calculates enrichment of within-group
  #     relationships relative to between-group relationships, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference
  results = matrix(nrow=length(geneset$names), ncol=10)
  rownames(results) = geneset$names
  colnames(results) = c("ingroup_mean", "outgroup_mean", "background_mean", "ingroup_n", 
          "outgroup_n", "background_n", "pvalue", "BH_pvalue", "pvalue_background", 
          "BH_pvalue_background")
  
  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    #cat(thisname, thissize, "\n")
    grouplist = geneset$matrix[i,1:thissize]
    
    if (!is.na(tag)) grouplist = sapply(grouplist, function (n) paste(tag, n, sep="_"))

    if (!is.na(idmap)) {
      grouplist = rownames(idmap)[which(idmap[,1] %in% grouplist)]
    }
    
    
    ingroup_ids = grouplist[which(grouplist %in% rownames(relationships))]
    outgroup_ids = rownames(relationships)[which(!rownames(relationships) %in% grouplist)]

    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup = relationships[ingroup_ids, ingroup_ids][upper.tri(relationships[ingroup_ids, ingroup_ids])]
      
      # pathway members versus all other components
      outgroup = relationships[ingroup_ids, outgroup_ids]
      
      # non-pathway members versus non-pathway members
      outgroup_2 = relationships[outgroup_ids, outgroup_ids][upper.tri(relationships[outgroup_ids, outgroup_ids])]
      
    }
    else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits = strsplit(mode, "-")[[1]]
      tag1 = bits[1]
      tag2 = bits[2]
      
      # filter to give a matrix where rows are tag1 and cols are tag2
      relationshipsa = relationships[grep(tag1, rownames(relationships)), grep(tag2, rownames(relationships))]
      
      relationshipsa = relationshipsa[which(rownames(relationshipsa) %in% ingroup_ids), 
                                        which(colnames(relationshipsa) %in% ingroup_ids)]
     
      # now match up cognates
      match1 = detag(rownames(relationshipsa))[detag(rownames(relationshipsa)) %in% detag(colnames(relationshipsa))]

      ingroup = sapply(match1, function (p) relationshipsa[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])
      
      # in the cognate comparison mode this no longer makes complete sense - that is,
      #    the comparison between groups isn't the same
      outgroup = relationships[ingroup_ids, outgroup_ids]
      
      # non-pathway members versus non-pathway members
      
      
      
      outgroup_2 = relationships[outgroup_ids, outgroup_ids][upper.tri(relationships[outgroup_ids, outgroup_ids])]
      #browser()
    }

    in_mean = mean(unlist(ingroup), na.rm=T)
    out_mean = mean(unlist(outgroup), na.rm=T)
    
    out_mean_2 = mean(unlist(outgroup_2), na.rm=T)
    
    pvalue = NA
    pvalue_2 = NA
    if (length(ingroup)>1) {
      pvalue = try(t.test(ingroup, outgroup)$p.value, silent=T);
      if (class(pvalue)=="try-error") pvalue = NA;
      pvalue_2 = try(t.test(ingroup, outgroup_2)$p.value, silent=T)
      if (class(pvalue_2)=="try-error") pvalue_2 = NA
    }
    
    delta = in_mean - out_mean
    this = c(in_mean, out_mean, out_mean_2, 
             length(ingroup), length(outgroup), length(outgroup_2), 
             pvalue, NA, 
             pvalue_2, NA)
    
    results[thisname,] = this
  }
  
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  results[,"BH_pvalue_background"] = p.adjust(results[,"pvalue_background"], method="BH")
  return(results)
}


ks_enrichment_vector <- function(groupgeneids, backgroundgenes, map_column=0, sort_col=1) {
  # returns a vector that is length(backgroundgenes) long of the
  #     enrichment of the groupgenes at each step in the ordered
  #     list of background genes
  if (map_column == 0) backgroundgenes = rownames(backgroundgenes[order(backgroundgenes[,sort_col]),])
  else backgroundgenes = backgroundgenes[order(backgroundgenes[,sort_col]),map_column]
  
  overallen = sum(groupgeneids %in% backgroundgenes)/length(backgroundgenes)
  
  result = c()
  for (i in 1:length(backgroundgenes)) {
    x = sum(groupgeneids %in% backgroundgenes[1:i])
    
    thisen = x/i
    
    result = c(result, thisen/overallen)
  }
  return(result)
}

difference_enrichment_in_relationships <- function(geneset, relationships1, relationships2, idmap=NA, tag=NA,
                                                   mode="original") {
  # for each category in geneset calculates enrichment of 
  #     relationships in matrix 1 versus those in matrix 2, where
  #     'relationships' are in the form of a square matrix (NxN) of
  #     continuous similarity metrics
  #     Currently uses a two-tailed t-test to assess this difference
  results = matrix(nrow=length(geneset$names), ncol=6)
  rownames(results) = geneset$names
  colnames(results) = c("group_mean_1", "group_mean_2", "ingroup1_n", 
                        "ingroup2_n", "pvalue", "BH_pvalue")
  
  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    #cat(thisname, thissize, "\n")
    grouplist = geneset$matrix[i,1:thissize]
    
    if (!is.na(tag)) grouplist = sapply(grouplist, function (n) paste(tag, n, sep="_"))
    
    if (!is.na(idmap)) {
      grouplist = rownames(idmap)[which(idmap[,1] %in% grouplist)]
    }
    
    # do these separately since the matrices could be different
    ingroup1_ids = grouplist[which(grouplist %in% rownames(relationships1))]
    ingroup2_ids = grouplist[which(grouplist %in% rownames(relationships2))]
    
    # original way to do this is to take all relationships
    if (mode == "original") {
      ingroup1 = relationships1[ingroup1_ids, ingroup1_ids][upper.tri(relationships1[ingroup1_ids, ingroup1_ids])]
      ingroup2 = relationships2[ingroup2_ids, ingroup2_ids][upper.tri(relationships2[ingroup2_ids, ingroup2_ids])]
    }
    else {
      # this does a calculation between cognate components in a multiomics dataset
      # we need to parse the mode and get the tags we want to compare
      # e.g. "cnv-txn", "txn-prot", "prot-phospho"*
      bits = strsplit(mode, "-")[[1]]
      tag1 = bits[1]
      tag2 = bits[2]
      
      # filter to give a matrix where rows are tag1 and cols are tag2
      relationships1a = relationships1[grep(tag1, rownames(relationships1)), grep(tag2, rownames(relationships1))]
      relationships2a = relationships2[grep(tag1, rownames(relationships2)), grep(tag2, rownames(relationships2))]
      
      relationships1a = relationships1a[which(rownames(relationships1a) %in% ingroup1_ids), 
                                        which(colnames(relationships1a) %in% ingroup1_ids)]
      relationships2a = relationships2a[which(rownames(relationships2a) %in% ingroup2_ids), 
                                        which(colnames(relationships2a) %in% ingroup2_ids)]
      
      # now match up cognates
      match1 = detag(rownames(relationships1a))[detag(rownames(relationships1a)) %in% detag(colnames(relationships1a))]
      match2 = detag(rownames(relationships2a))[detag(rownames(relationships2a)) %in% detag(colnames(relationships2a))]
      
      ingroup1 = sapply(match1, function (p) relationships1a[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])
      ingroup2 = sapply(match2, function (p) relationships2a[paste(tag1, p, sep="_"), paste(tag2, p, sep="_")])
      #browser()
    }
    
    in_mean1 = mean(unlist(ingroup1), na.rm=T)
    in_mean2 = mean(unlist(ingroup2), na.rm=T)
    
    pvalue = NA
    
    if (length(ingroup1)>1 && length(ingroup2)>1) {
      pvalue = try(t.test(ingroup1, ingroup2)$p.value, silent=T);
      if (class(pvalue)=="try-error") pvalue = NA;
    }
    
    this = c(in_mean1, in_mean2, 
             length(ingroup1), length(ingroup2), 
             pvalue, NA)
    
    results[thisname,] = this
  }
  
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  return(results)
}

correlation_enrichment <- function(geneset, abundance, mapping_column=NA, tag=NA) {
  allgenes = unique(unlist(as.list(geneset$matrix)))
  
  if (!is.na(tag)) allgenes = sapply(allgenes, function (n) paste(tag, n, sep="_"))
  
  # fixme: add support for an id column
  ids = rownames(abundance)
  cols = 1:ncol(abundance)
  map = NA
  if (!is.na(mapping_column)) {
    allgenes = rownames(abundance)[which(abundance[,1] %in% allgenes)]
    
    # NOTE: assumes that the mapping column is 1 and everything else
    #       is valid data- which may not be the case
    cols = 2:ncol(abundance)
    map = abundance
  }
  allgenes_present = allgenes[which(allgenes %in% ids)]
  allgenes_cor = cor(t(abundance[allgenes_present,cols]), use="p")
  
  return(list(enrichment=enrichment_in_relationships(geneset, allgenes_cor, idmap=map, tag=tag),
              corrmat=allgenes_cor))
}

correlation_comparison_enrichment <- function(geneset, abundance, set1, set2, mapping_column=NA, tag=NA, mode="original") {
  allgenes = unique(unlist(as.list(geneset$matrix)))
  
  if (!is.na(tag)) allgenes = sapply(allgenes, function (n) paste(tag, n, sep="_"))
  
  # fixme: add support for an id column
  ids = rownames(abundance)
  cols = 1:ncol(abundance)
  map = NA
  if (!is.na(mapping_column)) {
    allgenes = rownames(abundance)[which(abundance[,1] %in% allgenes)]
    
    # NOTE: assumes that the mapping column is 1 and everything else
    #       is valid data- which may not be the case
    cols = 2:ncol(abundance)
    map = abundance
  }
  cols1 = colnames(abundance)[which(colnames(abundance) %in% set1)]
  cols2 = colnames(abundance)[which(colnames(abundance) %in% set2)]
  
  allgenes_present = allgenes[which(allgenes %in% ids)]
  allgenes_cor1 = cor(t(abundance[allgenes_present,cols1]), use="p")
  allgenes_cor2 = cor(t(abundance[allgenes_present,cols2]), use="p")
  
  return(enrichment=difference_enrichment_in_relationships(geneset, allgenes_cor1, allgenes_cor2, idmap=map, tag=tag, mode=mode))
}


geneset_pca <- function(geneset, datatable, minsize=5, ncomp=1) {
  out = list()
  for (c in 1:ncomp) {
    results = matrix(nrow=length(geneset$names), ncol=ncol(datatable))
    rownames(results) = geneset$names
    colnames(results) = colnames(datatable)
    
    for (i in 1:length(geneset$names)) {
      thisname = geneset$names[i]
      thissize = geneset$size[i]
      thisdesc = geneset$desc[i]
      grouplist = geneset$matrix[i,1:thissize]
      
      in_group = datatable[grouplist[which(grouplist %in% rownames(datatable))],]
      
      if (nrow(in_group) < minsize) next
      
      prc = prcomp(in_group)
      #return(list(thisname, in_group, prc))
      
      results[thisname,] = prc$rotation[,c]
    }
    out = c(out, list(results))
  }
  return(out)
}


progressive_enrichment <- function(geneset, rankedlist, decreasing=T, ntimes=10, nstep=10, startat=10, endat=0.5) {
	# Progressive enrichment performs a simple Fisher's test on successively more and more
	#      of a ranked list. At each step the list is randomized ntimes and
	#      the percentage of times that a larger enrichment score is obtained is
	#      summed to provide an overall enrichment score
	
  rankedlist = sort(rankedlist, decreasing=decreasing)
  endatx = as.integer(length(rankedlist)*endat)
	presults_p = matrix(nrow=length(geneset$names), ncol=length(seq(startat, endatx, nstep)), dimnames=list(geneset$names, seq(startat, endatx, nstep)))
	presults_foldx = presults_p
	
	x = 1
	for (i in seq(startat, endatx, nstep)) {
    cat("Status", i, "\n")
		result = enrichment_in_groups(geneset, names(rankedlist[1:i]), names(rankedlist))
		
		randmat = sapply(1:ntimes, function (j) {randlist=sample(names(rankedlist), i); enrichment_in_groups(geneset, randlist, names(rankedlist))[,5]>result[,5]})
		presults_p[,x] = rowSums(randmat)/ntimes
		presults_foldx[,x] = result[,5]
		x = x + 1
	}
	return(list(pvalues=presults_p, foldx=presults_foldx))
}

partition_sets_enrichment <- function(partition_sets, geneset, abundance, mapping_column=NULL, abundance_column=NULL,
                                      sample_comparison=NULL, background_comparison=NULL) {
  results = list()
  
  for (i in 1:length(partition_sets$names)) {
    thisname = partition_sets$names[i]
    thissize = partition_sets$size[i]
    description = partition_sets$desc[i]
    grouplist = partition_sets$matrix[i,1:thissize]
    cat(thisname, " : ", description, "\n")
    
    result = enrichment_in_abundance(geneset, abundance[grouplist[which(grouplist %in% rownames(abundance))],],
                                     mapping_column=mapping_column, abundance_column=abundance_column,
                                     sample_comparison=sample_comparison, background_comparison=background_comparison)
    
    results = c(results, list(result))
  }
  names(results) = partition_sets$names
  return(results)
}

partition_sets_ksenrich <- function(partition_sets, geneset, abundance, mapping_column=NULL) {
  results = list()
  
  for (i in 1:length(partition_sets$names)) {
    thisname = partition_sets$names[i]
    thissize = partition_sets$size[i]
    description = partition_sets$desc[i]
    grouplist = partition_sets$matrix[i,1:thissize]
    cat(thisname, " : ", description, "\n")
    
    result = sapply(colnames(abundance), function (c) 
      enrichment_in_groups(geneset, background=abundance[grouplist[which(grouplist %in% rownames(abundance))],c],
                           method="ks", mapping_column=mapping_column)[,"Signed_AdjP"])
    
    
    results = c(results, list(result))
  }
  names(results) = partition_sets$names
  return(results)
}

enrichment_in_abundance <- function(geneset, abundance, mapping_column=NULL, abundance_column=NULL, 
                                    fdr=0, matchset=NULL, longform=F, sample_comparison=NULL, 
                                    background_comparison=NULL, min_p_threshold=NULL, tag=NA, sample_n=NULL) {
	# for each category in geneset calculates the abundance level
	#     of the genes/proteins in the category versus those
  #     not in the category and calculate a pvalue based on a
	#     two sided t test
  # If the sample_comparison variable is set then do a comparison between this
  #     abundance and sample comparison set which is vector of valid column (sample) ids
	results = matrix(nrow=length(geneset$names), ncol=8)
	rownames(results) = geneset$names
	colnames(results) = c("ingroup_mean", "outgroup_mean", "ingroup_n", "outgroup_n", "pvalue", "BH_pvalue", "count_above", "count_below")
	
  if (!is.null(mapping_column)) groupnames = unique(abundance[,mapping_column])
  
	for (i in 1:length(geneset$names)) {
	  thisname = geneset$names[i]
    if (!is.null(matchset) && matchset != thisname) next
		
		thissize = geneset$size[i]
		thisdesc = geneset$desc[i]
    #cat(thisname, thissize, "\n")
		grouplist = geneset$matrix[i,1:thissize]
    if (!is.na(tag)) grouplist = sapply(grouplist, function (n) paste(tag, n, sep="_"))
    
    if (!is.null(mapping_column)) {
      ingroupnames = grouplist[which(grouplist %in% groupnames)]
      outgroupnames = groupnames[which(!groupnames %in% grouplist)]
      
      if (!is.null(sample_n)) {
        if (sample_n > length(ingroupnames) || sample_n > length(outgroupnames)) {
          next  
        }
        #cat(sample_n, length(ingroupnames), length(outgroupnames), "\n")
        ingroupnames = sample(ingroupnames, sample_n)
        outgroupnames = sample(outgroupnames, sample_n)
      }
      
      ingroup = abundance[which(abundance[,mapping_column] %in% ingroupnames), 
                          abundance_column[which(abundance_column %in% colnames(abundance[,2:ncol(abundance)]))]]
      
      if (!is.null(sample_comparison)) outgroup = abundance[which(abundance[,mapping_column] %in% ingroupnames), 
                                                            sample_comparison[which(sample_comparison %in% colnames(abundance[,2:ncol(abundance)]))]]
      else outgroup = abundance[which(abundance[,mapping_column] %in% outgroupnames), abundance_column]
    }
		else {
      # I changed this to make it easier to use- may break old code!
		  # ingroup = abundance[grouplist[which(grouplist %in% names(abundance))]]
		  # outgroup = abundance[which(!names(abundance) %in% grouplist)]
		  ingroup = unlist(abundance[grouplist[which(grouplist %in% rownames(abundance))], 
		                             abundance_column[which(abundance_column %in% colnames(abundance))]])
      
      if (!is.null(sample_comparison)) outgroup = abundance[which(rownames(abundance) %in% grouplist), 
                                                            sample_comparison[which(sample_comparison %in% colnames(abundance))]]
      else outgroup = abundance[which(!rownames(abundance) %in% grouplist), 
                                abundance_column[which(abundance_column %in% colnames(abundance))]]
		}
    #cat(length(ingroup), length(outgroup), "\n")
		in_mean = mean(unlist(ingroup), na.rm=T)
		out_mean = mean(unlist(outgroup), na.rm=T)
		pvalue = NA
		if (length(ingroup)>1) {
			pvalue = try(t.test(unlist(ingroup), unlist(outgroup))$p.value, silent=T);
			if (class(pvalue)=="try-error") pvalue = NA;
		
			# we step through all components to calculate the number
			#    that are above/below the threshold
			if (!is.null(min_p_threshold)) {
			  count_above = 0
			  count_below = 0
			  # this works with sample_comparison
			  for (this_bit in rownames(ingroup)) {
			    petevalue = try(t.test(ingroup[this_bit,], outgroup[this_bit,])$p.value, silent=F);
			    if (class(petevalue)=="try-error" || is.na(petevalue)) { 
			      petevalue = NA;
			    }
			    else {
			      if (petevalue > min_p_threshold) count_above = count_above + 1
			      if (petevalue <= min_p_threshold) count_below = count_below + 1
			    }
			  } 
			  results[thisname,"count_above"] = count_above
			  results[thisname,"count_below"] = count_below
			}
		}
    
    delta = in_mean - out_mean
    if (fdr) {
      background = c()
      abundances = c(ingroup, outgroup)
      for (i in 1:fdr) {
        # randomly sample genes for fdr times
        ingroup = sample(1:length(abundances), length(ingroup))
        outgroup = which(!1:length(abundances) %in% ingroup)
        ingroup = abundances[ingroup]
        outgroup = abundances[outgroup]
        in_mean = mean(ingroup, na.rm=T)
        out_mean = mean(outgroup, na.rm=T)
        delta_r = in_mean - out_mean
        background = c(background, delta_r)
      }
      pvalue = sum(abs(background)>abs(delta))/length(background)
    }
    
		this = c(in_mean, out_mean, length(unlist(ingroup)), length(unlist(outgroup)), pvalue, NA)
    
		results[thisname,1:6] = this
	}
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  if (!is.null(matchset)) {
      results = results[matchset,]
      if (longform==T) results = list(results, ingroup, outgroup)
  }
	return(results)
}

group_membership_matrix <- function(geneset, pathwayname, sources, mode="original") {
	# for the named pathway creates a matrix that indicates membership
	#     from each source where source is a list of gene/protein lists
	members = geneset$matrix[which(geneset$names==pathwayname),1:geneset$sizes[which(geneset$names==pathwayname)]]
	
	results = matrix(nrow=length(members), ncol=length(sources))
	rownames(results) = members
	colnames(results) = names(sources)
	
	for (i in 1:length(sources)) {
		if (mode == "original") results[,i] = members %in% sources[i][[1]]
		else {
			results[members[which(members %in% names(sources[i][[1]]))],i] = sources[i][[1]][members[which(members %in% names(sources[i][[1]]))]]
		}
	}
	return(results)
}

cluster_enrichment <- function(database, clusters, background=NA, sigfilter=0.05) {
	x = length(clusters)-1
  if (is.na(background))  background = clusters[[x+1]]
	this = sapply(1:x, function (i) list(enrichment_in_groups(database, clusters[[i]], background)))
	outlist = list()
	for (i in 1:x) {
		these = this[[i]][which(this[[i]][,7]<sigfilter),]
		if (length(which(this[[i]][,7]<sigfilter)) == 1) {
			these = as.matrix(t(these))
			rownames(these)[[1]] = rownames(this[[i]])[which(this[[i]][,7]<sigfilter)]
		}
		outlist = c(outlist, list(these))
	}
	return(outlist)
}

make_gene_set_matrix <- function(database) {
  # takes a gene set database structure and returns
  #      a matrix that relates genes (rows) to their
  #      pathways (columns). This is useful because
  #      there can be redundance across different pathways
  #      and this structure can be used to figure out
  #      commonalities.
  ids = unique(as.vector(database$matrix))
  ids = ids[which(!ids == "null")]
  
  result = matrix(0, nrow=length(ids), ncol=length(database$names), dimnames=list(ids, database$names))
  
  x = 1
  for (this in database$names) {
    result[,this] = as.numeric(!rownames(result) %in% database$matrix[x,])
    x = x + 1
  }
                  
  return(result)
}

enrichment_redundancy_matrix <- function(geneset, dataframe=NA, enrichment_results=NA, significance_threshold=NA, 
                                         pathway_list=NA, method="jaccard") {
  # takes an enrichment table and the data frame that created it- with the rownames corresponding to gene names
  #    and calculates a pairwise overlap matrix- that is, for each pair of pathways what is the overlap
  #    between the genes represented in those pathways?
  
  if (!is.na(significance_threshold)) enrichment_results = enrichment_results[which(enrichment_results[,7]<significance_threshold),]
  if (!is.na(pathway_list)) {
      enrichment_results = matrix(nrow=length(pathway_list), ncol=1)
      rownames(enrichment_results) = pathway_list
  }
  
  result_matrix = matrix(nrow=nrow(enrichment_results), ncol=nrow(enrichment_results))
  rownames(result_matrix) = rownames(enrichment_results)
  colnames(result_matrix) = rownames(enrichment_results)
  
  for (i in rownames(enrichment_results)) {
    gene_list_i = get_pathway_information(geneset, i)$geneset
    if (!is.na(dataframe)) gene_list_i = gene_list_i[which(gene_list_i %in% rownames(dataframe))]
    
    for (j in rownames(enrichment_results)) {
      if (j==i) overlap = 1.0
      else {
        gene_list_j = get_pathway_information(geneset, j)$geneset
        if (!is.na(dataframe)) gene_list_j = gene_list_j[which(gene_list_j %in% rownames(dataframe))]
        
        # for first pass calculate jaccard index
        total = length(unique(c(gene_list_i, gene_list_j)))
        over = sum(gene_list_i %in% gene_list_j)
        if (method == "jaccard") overlap = over/total
        else overlap = over/length(gene_list_j)
      }
      result_matrix[i,j] = overlap
    }
  }
  return(result_matrix)
}

pairwise_overlap_enrichment <- function(geneset, datamatrix, condition1, condition2, mapping_column=NULL,
                                        subsample_components=NULL) {
  # iterate through a geneset and calculate the enrichment (abundance in condition1 v condition2)
  #      in the overlapping and unique gene sets for every combination
  
  results = list()
  
  overlap = matrix(nrow=length(geneset$names), ncol=length(geneset$names))
  rownames(overlap) = geneset$names
  colnames(overlap) = geneset$names
  
  # this wouldn't work in Python - I think it does in R
  # - that is the resulting matrices are new
  uniquea = unique_count = overlap_count = overlap
  
  for (path1 in geneset$names) {
    cat(path1, "\n")
    for (path2 in geneset$names) {
      if (path1 != path2) {
        this = overlap_gene_sets(geneset, path1, path2)
        
        if (is.null(mapping_column)) {
          theserows_o = rownames(datamatrix)[which(rownames(datamatrix) %in% this$overlap)]
          theserows_a = rownames(datamatrix)[which(rownames(datamatrix) %in% this$set1_unique)]
        }
        else {
          theserows_o = rownames(datamatrix)[which(datamatrix[,mapping_column] %in% this$overlap)]
          theserows_a = rownames(datamatrix)[which(datamatrix[,mapping_column] %in% this$set1_unique)]
        }
        
        if (!is.null(subsample_components)) {
          if (length(theserows_o)>subsample_components) theserows_o = sample(theserows_o, subsample_components)
          else theserows_o = c()
          if (length(theserows_a)>subsample_components) theserows_a = sample(theserows_a, subsample_components)
          else theserows_a = c()
        }
        #return(theserows_b)
        o = a = b = NA
        if (length(theserows_o)>5) o = try(t.test(datamatrix[theserows_o, condition1], datamatrix[theserows_o, condition2])$p.value, silent=T)
        if (length(theserows_a)>5) a = try(t.test(datamatrix[theserows_a, condition1], datamatrix[theserows_a, condition2])$p.value, silent=T)
        if (!class(o)=="try-error") overlap[path1,path2] = o
        if (!class(a)=="try-error") uniquea[path1,path2] = a
        
        unique_count[path1,path2] = length(theserows_a)
        overlap_count[path1,path2] = length(theserows_o)
        
        #cat(path1, path2, o, a, b, "\n", sep="\t")
      }
    }
  }
  return(list(overlap_pmat=overlap, unique_pmat=uniquea, 
              overlap_counts=overlap_count, unique_counts=unique_count))
}

ordered_list_enrichment_from_matrix <- function(geneset, datamatrix, starting_col=NULL, ending_col=NULL, idcolumn=NULL,
                                                threshold=0.02, background_list=NULL, enrichment_threshold=0.05) {
  
  # currently the thresholding is setup to handle fdrs/pvalues where lower is better
  #            but this should be extended
  if (is.null(starting_col)) starting_col = 1
  if (is.null(ending_col)) ending_col = ncol(datamatrix)
  if (is.null(background_list)) {
    if (is.null(idcolumn)) background_list = rownames(datamatrix)
    else background_list = datamatrix[,idcolumn]
  }
  
  thesenames = rownames(datamatrix)
  if (!is.null(idcolumn)) thesenames = datamatrix[,idcolumn]
  
  results = list()
  idlists = list()
  plists = list()
  for (c in colnames(datamatrix)[starting_col:ending_col]) {
      exlist = thesenames[which(datamatrix[,c]<threshold)]
      plist = rownames(datamatrix)[which(datamatrix[,c]<threshold)]
      enriched = enrichment_in_groups(geneset, exlist, background_list)
      results = c(results, list(enriched[which(enriched[,"Adjusted_pvalue"]<enrichment_threshold),]))
      idlists = c(idlists, list(exlist))
      plists = c(plists, list(plist))
  }
  names(results) = colnames(datamatrix)[starting_col:ending_col]
  names(idlists) = colnames(datamatrix)[starting_col:ending_col]
  names(plists) = colnames(datamatrix)[starting_col:ending_col]
  return(list(enrichment=results, idlists=idlists, names=plists))
}


enrichment_by_fishers <- function(group, background, annotation) {
  # calculate the fishers exact on a group of things versus a background
  #     for those things with a particular annotation (on another list)
  non_group = background[!background %in% group]
  
  group_annot = sum(group %in% annotation)
  group_annot_names = paste(intersect(group, annotation), collapse = ", ")
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
  
  return(list(fisher=ft, mat=test, foldx=fold, in_path_names=group_annot_names))
}

