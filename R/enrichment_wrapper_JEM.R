#' enrichment_wrapper
#'
#' enrichment_wrapper is a wrapper function that consolidates multiple enrichment methods.
#'
#' @param geneset is...
#' @param enrichment_method is...
#' @param datamatrix is...
#' @param id_column is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

enrichment_wrapper = function(geneset, enrichment_method, ...){
  .enrichment_wrapper(geneset, enrichment_method, ...)
}

.enrichment_wrapper = function(geneset, enrichment_method, datamatrix=NULL, id_column=NULL, primary_columns=NULL,
                               secondary_columns=NULL, threshold=NULL, minsize=5, mode=NULL,
                               idmap=NA,  fdr=0,
                               sample_comparison=NULL, min_p_threshold=NULL, sample_n=NULL,
                               enrichment_results=NA,
                               significance_threshold=NA, pathway_list=NA,
                               subsample_components=NULL, ntimes=100, greaterthan=TRUE){

  # JEM: we will remove the 'tag' parameter. This is used to indicate different types of 'omics data and we'll
  #      need to handle that, but will do it differently

  #check that geneset is of correct class
  #if(!inherits(geneset, "geneset_data")) stop("geneset must be of class 'geneset_data")

  #check that enrichment_method is one of the designated methods
  if(!is.element(enrichment_method, c("correlation_comparison_enrichment", "correlation_enrichment",
                                      "difference_enrichment_in_relationships", "enrichment_in_abundance",
                                      "enrichment_by_fishers", "enrichment_by_ks", "enrichment_in_relationships",
                                      "enrichment_redundancy_matrix", "pairwise_overlap_enrichment",
                                      "permute_enrichment_in_groups")))
    stop("enrichment_method must be one of the methods designated in the function documentation")

  #checking each enrichment method
  # JEM: changed here - moved these around so that the most important ones are first
  #      There are two different modes this can be run in and we can just have the user specify in the enrichment_method
  if(enrichment_method == "enrichment_by_ks"){

    message("'enrichment_by_ks' has the following optional arguments: datamatrix=NULL, minsize=5, id_column=NULL, primary_columns=NULL")

    # datamatrix will replace background - but we can do that at this level and not change the underlying code
    # primary_columns will replace abundance_column  - but again, at this level
    # minsize will remain the same
    # id_column will replace mapping_column
    # randomize might change - TBD

    # TODO: add a validation that primary_columns is one column id for this case
    if(!is.null(primary_columns)){
      if(length(primary_columns) > 1){stop("'primary_columns' must be a string of length 1")}
    }

    if(!is.null(primary_columns) & !is.null(datamatrix)){
      # TODO: add a validation that datamatrix has a column named primary_columns
      if(!is.element(primary_columns, colnames(datamatrix))){stop("'primary_columns' must be a column name of 'datamatrix'")}
    }

    # JEM - notes for me: this application takes a background (matrix) and specifies which abundance column for calculating ks enrichment
    result = enrichment_in_groups(geneset=geneset, background=datamatrix,
                                  method="ks", minsize=minsize, mapping_column=id_column,
                                  abundance_column=primary_columns)
  }

  else if(enrichment_method == "enrichment_by_fishers"){

    message("'enrichment_by_fishers' has the following optional arguments: datamatrix=NULL, minsize=5, id_column=NULL, primary_columns=NULL, greaterthan=TRUE")


    # JEM - notes for me: this application takes a list of targets and a list of background genes and applies fisher's exact
    #     if we update this to take a datamatrix and speciify a column then we need to also allow a threshold to be specified
    #     which is probably OK for this application -
    # datamatrix will replace background
    # primary_columns will replace abundance column
    # minsize will remain the same
    # id_column will replace mapping_column
    # randomize isn't implemented for this application

    if(!is.null(primary_columns)){
      # We'll do the list splitting here
      # TODO: add a validation that the primary_columns is one column for this case
      if(length(primary_columns) > 1){stop("'primary_columns' must be a string of length 1")}
    }

    if(!is.null(primary_columns) & !is.null(datamatrix)){
      # TODO: add a validation that the datamatrix has a column with the primary_columns name
      if(!is.element(primary_columns, colnames(datamatrix))){stop("'primary_columns' must be a column name of 'datamatrix'")}
    }

    background = rownames(datamatrix)

    # TODO: we should allow the user to specify if they want it greater than or less than for the threshold
    if(greaterthan == FALSE & !is.null(threshold)){
      targets = rownames(datamatrix[which(datamatrix[,primary_columns]<threshold),])
      }
    else if(greaterthan == TRUE & !is.null(threshold)){
      targets = rownames(datamatrix[which(datamatrix[,primary_columns]>threshold),])
    }

    # TODO: add a validation that targets is non-empty
    if(length(targets) == 0){stop("please adjust 'threshold' argument, there are no target genes captured by current threshold value")}

    # TODO: add a validation that id_column exists in datamatrix
    if(!is.null(id_column) & !is.null(datamatrix)){
      if(!is.element(id_column, colnames(datamatrix))){stop("'id_column' must be a column name of 'datamatrix'")}
    }

    if(!is.null(id_column)) {
      background = unique(datamatrix[background,])
      targets = unique(datamatrix[targets,])
    }

    result = enrichment_in_groups(geneset=geneset, targets=targets, background=background,
                                  method="fishers", minsize=minsize)

  }

  else if(enrichment_method == "enrichment_in_abundance"){

    if(is.null(datamatrix)){stop("'datamatrix' argument is required")}
    message("'enrichment_in_abundance' has the following optional arguments: id_column=NULL, primary_columns=NULL,fdr=0, secondary_columns=NULL,
                                          min_p_threshold=NULL, sample_n=NULL")

    result = enrichment_in_abundance(geneset=geneset, abundance=datamatrix, mapping_column=id_column,
                                     abundance_column=primary_columns, sample_comparison=secondary_columns,
                                     fdr=fdr, min_p_threshold=min_p_threshold, sample_n=sample_n)

    # matchset=matchset, - this is a specific flag that's used for debugging
    # min_p_threshold is a parameter that will change the output - adding a couple of columns

    result = as.data.frame(result)
  }

  else if(enrichment_method == "enrichment_in_relationships"){

    if(is.null(datamatrix)){stop("'datamatrix' argument is required")}
    message("'enrichment_in_relationships' has the following optional arguments: idmap=NA")

    #if(is.null(mode)){mode = 'original'} - this pertains to multiomics comparisons and we'll need to figure out how
    #                                       to recode this -

    result = enrichment_in_relationships(geneset=geneset, relationships=datamatrix, idmap=idmap)
  }

  # JEM- this function needs to be updated. Right now it's very simple and has the user define what set they
  #           want to use for correlation- but we should allow the user to specify a primary_columns parameter
  else if(enrichment_method == "correlation_enrichment"){

    if(is.null(datamatrix)){stop("'datamatrix' argument is required")}
    message("'correlation_enrichment' has the following optional arguments:  id_column=NA, tag=NA")
    if(is.null(id_column)){id_column = NA}

    temp_result = correlation_enrichment(geneset=geneset, abundance=datamatrix, mapping_column=id_column)

    result = as.data.frame(temp_result$enrichment)
    attr(result, "corrmat") = temp_result$corrmat
  }

####IG this function calls "difference_enrichment_in_relationships" which is not an active enrichment method at this point
#
#  else if(enrichment_method == "correlation_comparison_enrichment"){
#
#    if(is.null(datamatrix)|is.null(primary_columns)|is.null(secondary_columns)) {stop("'datamatrix', 'primary_columns' and 'secondary_columns' arguments are required")}
#    message("'correlation_comparison_enrichment' has the following optional arguments: id_column=NA, tag=NA, mode='original'")
#    if(is.null(mode)){mode = 'original'}
#    if(is.null(id_column)){id_column = NA}
#
#    result = correlation_comparison_enrichment(geneset=geneset, abundance=datamatrix, set1=primary_columns, set2=secondary_columns,
#                                               mapping_column=id_column, mode=mode)
#  }
####

  ### JEM - I think we should maybe not have this following function in the wrapper
  ###       the reason being that this takes two data matrices instead of one and there's not
  ###       a graceful way to do this in the context of what we have for the other functions here

  # else if(enrichment_method == "difference_enrichment_in_relationships"){
  #
  #   if(is.null(relationships1|is.null(relationships2))){stop("'relationships1' and 'relationships2' arguments are required")}
  #   message("'difference_enrichment_in_relationships' has the following optional arguments: idmap=NA, tag=NA,mode='original'")
  #   if(is.null(mode)){mode = 'original'}
  #
  #   result = difference_enrichment_in_relationships(geneset=geneset, relationships1=primary_columns, relationships2=secondary_columns,
  #                                                   idmap=idmap, tag=tag, mode=mode)
  # }

  ### JEM - also for this one. We should leave it exposed- but it's something quite different that it's doing

  # else if(enrichment_method == "enrichment_redundancy_matrix"){
  #
  #   message("'enrichment_redundancy_matrix' has the following optional arguments: dataframe=NA, enrichment_results=NA, significance_threshold=NA, pathway_list=NA, method='jaccard'")
  #   if(is.null(method)){method = 'jaccard'}
  #
  #   result = enrichment_in_relationships(geneset=geneset, dataframe=dataframe, enrichment_results=enrichment_results,
  #                                        significance_threshold=significance_threshold, pathway_list=pathway_list, method=method)
  # }

  ### JEM - the following function is useful and has a similar format as other wrapped functions
  ###       but will be very slow, especially for large genesets since it does a pairwise enrichment for all genesets

  else if(enrichment_method == "pairwise_overlap_enrichment"){

    if(is.null(datamatrix)|is.null(primary_columns)|is.null(secondary_columns)){stop("'datamatrix', 'primary_columns', 'secondary_columns' arguments are required")}
    message("'pairwise_overlap_enrichment' has the following optional arguments: id_column=NULL, subsample_components=NULL")

    temp_result = pairwise_overlap_enrichment(geneset=geneset, datamatrix=datamatrix, condition1=primary_columns, condition2=secondary_columns,
                                         mapping_column=id_column, subsample_components=subsample_components)
    result = as.data.frame(temp_result$overlap_pmat)
    attr(result, "unique_pmat") = as.data.frame(temp_result$unique_pmat)
    attr(result, "overlap_counts") = as.data.frame(temp_result$overlap_counts)
    attr(result, "unique_counts") = as.data.frame(temp_result$unique_counts)


  }

  # JEM- we need to do this - but probably should make it work for all the different methods in the same way
  #      Right now keep as a future development.
  # else if(enrichment_method == "permute_enrichment_in_groups"){
  #
  #   message("'permute_enrichment_in_groups' has the following optional arguments: targets=NULL, background=NULL, method='fishers', minsize=5, mapping_column=NULL, abundance_column=NULL, ntimes=100")
  #   if(is.null(method)){method = 'fishers'}
  #
  #   result = permute_enrichment_in_groups(geneset=geneset, targets=targets, background=background, method=method,
  #                                         minsize=minsize, mapping_column=mapping_column, abundance_column=abundance_column, ntimes=ntimes)
  # }


 return(result)


}

