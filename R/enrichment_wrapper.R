#work on wrapper function for enrichment functions

enrichment_wrapper = function(geneset, enrichment_method, ...){
  .enrichment_wrapper(geneset, enrichment_method, ...)
}

.enrichment_wrapper = function(geneset, enrichment_method, abundance=NULL, set1=NULL, set2=NULL, mapping_column=NULL, tag=NA, mode=NULL,
                               relationships1=NULL, relationships2=NULL, idmap=NA, abundance_column=NULL, fdr=0, matchset=NULL, longform=F,
                               sample_comparison=NULL, background_comparison=NULL, min_p_threshold=NULL, sample_n=NULL, targets=NULL,
                               background=NULL, method=NULL, minsize=5, randomize=F, relationships=NULL, dataframe=NA, enrichment_results=NA,
                               significance_threshold=NA, pathway_list=NA, datamatrix=NULL, condition1=NULL, condition2=NULL,
                               subsample_components=NULL, ntimes=100){

  #check that geneset is of correct class
  #if(!inherits(geneset, "geneset_data")) stop("geneset must be of class 'geneset_data")

  #check that enrichment_method is one of the designated methods
  if(!is.element(enrichment_method, c("correlation_comparison_enrichment", "correlation_enrichment", "difference_enrichment_in_relationships", "enrichment_in_abundance", "enrichment_in_groups", "enrichment_in_relationships", "enrichment_redundancy_matrix", "pairwise_overlap_enrichment", "permute_enrichment_in_groups"))) stop("enrichment_method must be one of the methods designated in the function documentation")

  #checking each enrichment method
  if(enrichment_method == "correlation_comparison_enrichment"){

    if(is.null(abundance)|is.null(set1)|is.null(set2)) {stop("'abundance', 'set1' and 'set2' arguments are required")}
    message("'correlation_comparison_enrichment' has the following optional arguments: mapping_column=NA, tag=NA, mode='original'")
    if(is.null(mode)){mode = 'original'}
    if(is.null(mapping_column)){mapping_column = NA}

    result = correlation_comparison_enrichment(geneset=geneset, abundance=abundance, set1=set1, set2=set2, mapping_column=mapping_column, tag=tag, mode=mode)
  }

  else if(enrichment_method == "correlation_enrichment"){

    if(is.null(abundance)){stop("'abundance' argument is required")}
    message("'correlation_enrichment' has the following optional arguments: mapping_column=NA, tag=NA")
    if(is.null(mapping_column)){mapping_column = NA}

    result = correlation_enrichment(geneset=geneset, abundance=abundance, mapping_column=mapping_column, tag=tag)
  }

  else if(enrichment_method == "difference_enrichment_in_relationships"){

    if(is.null(relationships1|is.null(relationships2))){stop("'relationships1' and 'relationships2' arguments are required")}
    message("'difference_enrichment_in_relationships' has the following optional arguments: idmap=NA, tag=NA,mode='original'")
    if(is.null(mode)){mode = 'original'}

    result = difference_enrichment_in_relationships(geneset=geneset, relationships1=relationships1, relationships2=relationships2, idmap=idmap, tag=tag, mode=mode)
  }

  else if(enrichment_method == "enrichment_in_abundance"){

    if(is.null(abundance)){stop("'abundance' argument is required")}
    message("'enrichment_in_abundance' has the following optional arguments: mapping_column=NULL, abundance_column=NULL,fdr=0, matchset=NULL, longform=F, sample_comparison=NULL, background_comparison=NULL, min_p_threshold=NULL, tag=NA, sample_n=NULL")

    result = enrichment_in_abundance(geneset=geneset, abundance=abundance, mapping_column=mapping_column, abundance_column=abundance_column, fdr=fdr, matchset=matchset, longform=longform, sample_comparison=sample_comparison, background_comparison=background_comparison, min_p_threshold=min_p_threshold, tag=tag, sample_n=sample_n)
  }

  else if(enrichment_method == "enrichment_in_groups"){

    message("'enrichment_in_groups' has the following optional arguments: targets=NULL, background=NULL, method='fishers', minsize=5, mapping_column=NULL, abundance_column=NULL, randomize=F")
    if(is.null(method)){method = 'fishers'}

    result = enrichment_in_groups(geneset=geneset, targets=targets, background=background, method=method, minsize=minsize, mapping_column=mapping_column, abundance_column=abundance_column, randomize=randomize)
  }

  else if(enrichment_method == "enrichment_in_relationships"){

    if(is.null(relationships)){stop("'relationships' argument is required")}
    message("'enrichment_in_relationships' has the following optional arguments: idmap=NA, tag=NA, mode='original'")
    if(is.null(mode)){mode = 'original'}

    result = enrichment_in_relationships(geneset=geneset, relationships=relationships, idmap=idmap, tag=tag, mode=mode)
  }

  else if(enrichment_method == "enrichment_redundancy_matrix"){

    message("'enrichment_redundancy_matrix' has the following optional arguments: dataframe=NA, enrichment_results=NA, significance_threshold=NA, pathway_list=NA, method='jaccard'")
    if(is.null(method)){method = 'jaccard'}

    result = enrichment_in_relationships(geneset=geneset, dataframe=dataframe, enrichment_results=enrichment_results, significance_threshold=significance_threshold, pathway_list=pathway_list, method=method)
  }

  else if(enrichment_method == "pairwise_overlap_enrichment"){

    if(is.null(datamatrix)|is.null(condition1)|is.null(condition2)){stop("'datamatrix', 'condition1', 'condition2' arguments are required")}
    message("'pairwise_overlap_enrichment' has the following optional arguments: mapping_column=NULL, subsample_components=NULL")

    result = pairwise_overlap_enrichment(geneset=geneset, datamatrix=datamatrix, condition1=condition1, condition2=condition2, mapping_column=mapping_column, subsample_components=subsample_components)
  }

  else if(enrichment_method == "permute_enrichment_in_groups"){

    message("'permute_enrichment_in_groups' has the following optional arguments: targets=NULL, background=NULL, method='fishers', minsize=5, mapping_column=NULL, abundance_column=NULL, ntimes=100")
    if(is.null(method)){method = 'fishers'}

    result = permute_enrichment_in_groups(geneset=geneset, targets=targets, background=background, method=method, minsize=minsize, mapping_column=mapping_column, abundance_column=abundance_column, ntimes=ntimes)
  }


 return(result)


}

