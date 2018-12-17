# Set analysis for examining overlap between sets of DE genes
# JEM 6/15/2010

gene_set_overlap <- function(inmatrix, colsets=NA, filter=1.5) {
	# this will simply examine the overlap in successively more sets
	#      but will ignore which specific sets are overlapping. That is
	# 	   it will identify genes that are DE in ANY N conditions.
	#      
	#  inmatrix : matrix of values to examine
	#  colsets  : a list of vectors describing groups to look for DE in
	#  filter   : what DE do we consider to be DE?
	#
	# To create a colset list, e.g.: 
	#     conds = list(colnames(ratios)[45:59], colnames(ratios)[61:66], colnames(ratios)[77:86])
	outmatrix = matrix(nrow=nrow(inmatrix), ncol=length(colsets))
	rownames(outmatrix) = rownames(inmatrix)
	ids = c()
	for (gene in rownames(inmatrix)) {
		for (i in 1:length(colsets)) {
			if (sum(abs(inmatrix[gene,colsets[[i]]])>filter)) {
				outmatrix[gene,i] = TRUE	
			} else {
				outmatrix[gene,i] = FALSE	
			}
		}
		ids = c(ids, sum(outmatrix[gene,])>0)
	}
	# filter the matrix so it only has shared DE genes
	outmatrix = outmatrix[ids,]
	return(outmatrix)
}

random_gene_set_overlap <- function(inmatrix, colsets=NA, filter=1.5, N=NA, rep=100) {
	# this just returns the number of genes DE in N or more conditions
	#   if N==NA, N defaults to all conditions

	# this returns a list of the number of genes that are DE in each condition
	de_list = c()
	for (cond in colsets) {
		x = sapply(rownames(inmatrix), function (r) sum(abs(inmatrix[r,cond[1]]))>filter)
		de_list = c(de_list, sum(x))
		#cat("Condition", str(cond), "\n")
	}
	
	#cat("here\n")
	#return(de_list)
	results = c()
	
	for (i in 1:rep) {
		# now randomly select genes
		de_names = sapply(de_list, function (dn) sample(rownames(inmatrix), dn))
		
		# count the number of occurrences of each identifier
		de_occ = rle(sort(unlist(de_names)))
				
		if (is.na(N)) {
			N = length(de_list)	
		}
		result = sum(de_occ$lengths >= N)
		results = c(results, result)
	}
	return(results)
}

random_overlap <- function(items, counts=c(10,50,100), N=3,rep=100) {
	# this just returns the number of items that match when randomly drawn
	#     with the numbers indicated in counts
	results = c()
	
	for (i in 1:rep) {
		# now randomly select genes
		bats = c()
		
		for (this in counts) {
			bits = sample(items, this)
			bats = c(bats, bits)
		}
			
		# count the number of occurrences of each item
		de_occ = rle(sort(bats))
						
		result = sum(de_occ$lengths >= N)
		results = c(results, result)
	}
	return(results)
}

pvalue_from_random_gene_set_overlap_db <- function(inmatrix_1, inmatrix_2, names_1, names_2, rep=1000, longform=F) {
	# calculates a pvalue from random iterations for overlap between two sets
	x = sum(names_1 %in% names_2)
	
	pvect = random_gene_set_overlap_db(inmatrix_1, inmatrix_2, length(names_1), length(names_2), rep=rep)
	
	pvalue = sum(pvect>=x)/length(pvect)
	
	if (longform) return(list(x=x, pvalue=pvalue, pvect=pvect, n1=length(names_1), n2=length(names_2)))
	
	return(pvalue)
}

random_gene_set_overlap_db <- function(inmatrix_1, inmatrix_2, n1, n2, rep=100) {
	# This is a variant of the above to test 2 lists of different lengths
	#      to see how much of an overlap there is- and replicate randomly a bunch of times
	
	if (is.matrix(inmatrix_1)) {
		names1 = rownames(inmatrix_1)
	} else names1 = names(inmatrix_1)
	
	if (is.matrix(inmatrix_2)) {
		names2 = rownames(inmatrix_2)
	} else names2 = names(inmatrix_2)
	
	
	results = c()
	
	for (i in 1:rep) {
		sn1 = sample(names1, n1)
		sn2 = sample(names2, n2)
		x = sum(sn1 %in% sn2)
		
		results = c(results, x)
	}
	return(results)
}

random_gene_set_overlap_x3 <- function(list_1, list_2, list_3, n1, n2, n3, rep=100) {
	# This is a variant of the above to test 3 lists of different lengths
	#      to see how much of an overlap there is- and replicate randomly a bunch of times
	
	results = c()
	
	for (i in 1:rep) {
		sn1 = sample(list_1, n1)
		sn2 = sample(list_2, n2)
		sn3 = sample(list_3, n3)
		x1 = sn1[which(sn1 %in% sn2)]
		x = sum(x1 %in% sn3)
		
		results = c(results, x)
	}
	return(results)
}


gene_set_heatmap <- function(inmatrix, colsets=NA, filter=1.5) {
	# this will generate a matrix identical to the one above, but
	#      with the max(abs(expression)) value instead of T/F
	#      
	#  inmatrix : matrix of values to examine
	#  colsets  : a list of vectors describing groups to look for DE in
	#  filter   : what DE do we consider to be DE?
	#
	# To create a colset list, e.g.: 
	#     conds = list(colnames(ratios)[45:59], colnames(ratios)[61:66], colnames(ratios)[77:86])
	outmatrix = matrix(nrow=nrow(inmatrix), ncol=length(colsets))
	rownames(outmatrix) = rownames(inmatrix)
	ids = c()
	for (gene in rownames(inmatrix)) {
		for (i in 1:length(colsets)) {
			if (sum(abs(inmatrix[gene,colsets[[i]]])>filter)) {
				value = inmatrix[gene,colsets[[i]]][which(abs(inmatrix[gene,colsets[[i]]]) == max(abs(inmatrix[gene,colsets[[i]]])))]
				#cat(gene, i, value, "\n")
				outmatrix[gene,i] = value
			} else {
				outmatrix[gene,i] = NA	
			}
		}
		ids = c(ids, sum(outmatrix[gene,])>0)
	}
	# filter the matrix so it only has shared DE genes
	#outmatrix = outmatrix[ids,]
	return(outmatrix)	
}

bottleneck_overlap <- function(set_1, set_2, percentage=0.2) {
	# takes two matrices/data frames that are ordered and returns a
	# matrix of overlapping rows by rowname
	ret_set = set_1[which(rownames(set_1)[1:(nrow(set_1)*percentage)] %in% rownames(set_2)[1:(nrow(set_2)*percentage)]),]
	return(ret_set)
}


plot_expression <- function(ratios, rowname, pfactor=1, xvector=NA, xfactor=0, title=NA, condrange=c(), trendline=FALSE, trendcolor="red", color="grey", legend=FALSE, yr=NA, mode="new") {
	# expand the y margin lower
	par(mar=c(15,4.1,4.2,2.1))
	
	#if (!class(title) == "character") {
	#	if (class(rowname) == "character") {
	#		title = rowname
	#	} else {
	#		title = rownames(ratios)[rowname]
	#	}
	#}
	if (!is.vector(rowname)) {
		rowname = c(rowname)
		if (length(condrange) == 0) {
			condrange = c(1:length(ratios[rowname,]))
		}
	} else {
		if (length(condrange) == 0) {
			condrange = c(1:length(ratios[rowname[1],]))
		}
	}
	
	name= rowname[1]
	if (mode == "new") {
		if (is.na(yr)) {
			if (is.na(xvector)) {
				#cat(name, str(ratios[name, condrange]),"\n")
				# plot the observed values in red with no axes
				plot(ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1)
			} else {
				plot(xvector, ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1)
			}
			
		} else {
			if (is.na(xvector)) {
				plot(ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1, ylim=yr)
			} else {
				plot(xvector, ratios[name,condrange], type="l", col=color, axes=FALSE, xlab="", ylab="expression values", main=title, lwd=1, ylim=yr)
			}
		}
	}
	for (name in rowname[2:(length(rowname))]) {
		if (is.na(xvector)) {
			lines(ratios[name,condrange], type="l", col=color)
		} else {
			lines(xvector, ratios[name,condrange], type="l", col=color)
		}
	}
	if (trendline == TRUE) {
		trend = sapply(condrange, function (c) mean(as.numeric(ratios[rowname,c])))
		if (is.na(xvector)) {
			lines(trend, type="l", col=trendcolor)
		} else {
			lines(xvector, trend, type="l", col=trendcolor)
			}
		}
	# add categories on the x axis
	if (is.na(xvector)) {
		axis(1, labels=names(ratios[name,condrange]), at=1:length(ratios[name,condrange]), las=2)
	} else {
		axis(1, labels=names(ratios[name,condrange]), at=xvector, las=2)
	}
	axis(2)
}

