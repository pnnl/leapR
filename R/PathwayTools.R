#Pathway tools

pathway_activity <- function(data, geneset, method="mean-abundance") {
  # this function provides a pathway summary for each pathway in the gene_sets
  #      (msigdb format) using one of several methods. A matrix is returned
  #      with the same number of columns (conditions) as the input data matrix
  #      and each row being a different pathway
  
  result = matrix(nrow=length(geneset$names), ncol=ncol(data))
  rownames(result) = geneset$names
  colnames(result) = colnames(data)
  
  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    grouplist = geneset$matrix[i,1:thissize]
    
    inmat = data[grouplist[which(grouplist %in% rownames(data))],]
    
    # if this is empty (i.e. no matches) or it only has one member - returning a vector
    #     then we'll leave it blank and not count it as anything
    if (class(inmat) == "matrix") {
      result[thisname,] = colMeans(inmat, na.rm=T)
    }
  }
  return(result)
}