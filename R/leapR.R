#' leapR
#'
#' leapR is a wrapper function that consolidates multiple enrichment methods.
#'
#' @param geneset is a list of four vectors, gene names, gene descriptions, gene sizes and a matrix of genes. It represents .gmt format pathway files.
#' @param enrichment_method is a character string specifying the method of enrichment to be performed, one of: "enrichment_comparison", "enrichment_in_order", "enrichment_in_sets", "enrichment_in_pathway", "correlation_enrichment", "enrichment_in_relationships".
#' @return data frame with results
#' @import Biobase
#' @param ... further arguments
#'
#' 
#' @details 
#' Further arguments and enrichment method optional argument information: \cr
#' \tabular{ll}{
#' eset \tab Is an \code{ExpressionSet} of expression data, with features as rows and \emph{n} sample/conditions as columns. 
#' The \code{Annotation} field ideally describes the data type (i.e. proteomics, phosphoproteomics), the \code{featureData} field describes any mapping identifiers and the \code{phenoData} field describes any 
#' phenotyptic data. This is an required for all active enrichment methods with the exception of 'enrichment_in_relationships'. \cr
#' \cr
#' id_column \tab Is a character string, present in the \code{featureData} slot, that is used to specify a column for identifiers. If missing, the rownames of the ExpressionSet will be used. 
#' \cr
#' primary_columns \tab Is a character vector composed of column names from \code{eset} (either in the `exprs` or in the `featureData`), 
#' that specifies a set of primary columns to calculate enrichment on. 
#' The meaning of this varies according to the enrichment method used - see the descriptions for each method below. 
#' This is an optional argument used with 'enrichment_in_order', 'enrichment_in_sets', and 'enrichment_comparison' methods. \cr
#' \cr
#' secondary_columns \tab Is a character vector of column names. This is an optional argument used with 'enrichment_comparison' methods. \cr
#' \cr
#' threshold \tab Is a numeric value, an optional argument used with 'enrichment_in sets' method which filters out abundance values 
#' either above or below it. \cr
#' \cr
#' greaterthan \tab Is a logical value that defaults to TRUE, it's used with 'enrichment_in_sets' method. 
#' When set to TRUE, genes with abundance data above the \code{threshold} argument are kept. 
#' When set to FALSE genes with abundance data below the \code{threshold} argument are kept. 
#' This is an optional argument used with 'enrichment_in_sets' method. \cr
#' \cr
#' minsize \tab Is a numeric value, an optional argument used with 'enrichment_in_sets' and 'enrichment_in_order". \cr
#' \cr
#' idmap \tab Is...??. This is an optional argument used with 'enrichment_in_relationships' method. \cr
#' \cr
#' fdr \tab A numerical value which specifies how many times to randomly sample genes to calculate an empirical false discovery rate, is an optional argument used with 'enrichment_comparison' method. \cr
#' \cr
#' min_p_threshold \tab Is a numeric value, a lower p-value threshold and is an optional argument used with 'enrichment_comparison' method. \cr
#' \cr
#' sample_n \tab Is a way to subsample the number of components considered for each calculation randomly. This is an optional argument used with 'enrichment_comparison' method. \cr
#' \cr
#' }
#' 
#' \strong{Enrichment Methods:}
#' \cr
#' \cr
#' enrichment_comparison
#' \cr
#' Compares the distribution of abundances between two sets of conditions for each pathway using a t test. For each pathway in \code{geneset}
#' uses a t test to compare the distribution of abundance values/numbers in \code{eset} \code{primary_columns} with those in
#' \code{eset} \code{secondary_columns}. Lower p-values for pathways indicate that the expression of the pathway is
#' significantly different between the set of conditions in primary_columns and the set of conditions in secondary_columns.
#' Optionally, users can specify \code{fdr} which will calculate an empirical p-value by randomizing abdunances
#' \code{fdr} number of times. If the \code{min_p_threshold} is specified the method will only return pathways with an
#' adjusted p-value lower than the specified threshold. If \code{sample_n} is specified the method will subsample the 
#' pathway members to the specified number of components.
#' \cr \cr
#' enrichment_in_order
#' \cr
#' Calculates enrichment of pathways based on a ranked list using the Kologmorov-Smirnov test. 
#' For each pathway in \code{geneset} uses a Kolgmorov-Smirnov test for rank order to test if the distribution
#' of ranked abundance values in the \code{eset} \code{primary_columns} is significant relative to a random
#' distribution. Note that currently \code{primary_columns} only accepts a single column for this method.
#' \cr \cr
#' enrichment_in_sets
#' \cr
#' Calculates enrichment in pathway membership in a list (e.g. highly differential proteins) relative to background using Fisher's exact test.
#' For each pathway in \code{geneset} uses a Fisher's exact test over- or under- representation of a list
#' of components specified. If \code{targets} are specified this must be a vector of identifiers to serve
#' as the target list for comparison. If \code{eset} and \code{primary_columns} are specified then \code{threshold} specifies
#' a threshold value for determining the target list of components to test. Specifying \code{greaterthan} to be False
#' will result in components with values lower than the specified \code{threshold}. If \code{eset} is
#' a data frame or matrix, the background used for calculation will be taken as the rownames of \code{eset}
#' \cr \cr
#' enrichment_in_pathway
#' \cr
#' Compares the distribution of abundances in a pathway with the background distribution of abundances using a t test
#' For each pathway in \code{geneset} calculates the signficance of the difference between the abundances
#' from pathway members versus abundance of non-pathway members in the set of conditions specified by \code{primary_columns}.
#' Optionally, users can specify \code{fdr} which will calculate an empirical p-value by randomizing abdunances
#' \code{fdr} number of times. If the \code{min_p_threshold} is specified the method will only return pathways with an
#' adjusted p-value lower than the specified threshold. If \code{sample_n} is specified the method will subsample the 
#' pathway members to the specified number of components.
#' \cr \cr
#' correlation_enrichment
#' Calculates the enrichment of a pathway based on correlation between pathway members across conditions versus correlation between members not in the pathway.
#' For each pathway in \code{geneset} calculates the pairwise correlation between all pathway members and non-pathway members
#' across the specified \code{primary_columns} conditions in \code{eset}. Note that for large matrices this can take a long
#' time. A p-value is calculated based on comparing the correlation within the members of a pathway with the correlation
#' values between members of the pathway and non-members of the pathway.
#' \cr \cr
#' enrichment_in_relationships
#' Calculates the enrichment of a pathway in specified interactions relative to non-pathway members. For each pathway in \code{geneset}
#' calculates the enrichment in relationships as defined by an adjacency matrix provided in \code{eset}. An adjacency matrix
#' is a square matrix to provide pairwise relationships between components (genes, proteins) as derived from e.g. correlation as
#' \code{correlation_enrichment}.
#'
#' @examples
#'         library(leapR)
#'
#'         # read in the example abundance data
#'         # read in the example transcriptomic data
#'         tdata <- download.file("https://figshare.com/ndownloader/files/55781153",method='libcurl',destfile='transData.rda')
#'         load('transData.rda')
#'         p <- file.remove("transData.rda")
#'
#'         # read in the pathways
#'         data("ncipid")
#'
#'         # read in the patient groups
#'         data("shortlist")
#'         data("longlist")
#'         
#'         # use enrichment_comparison to calculate enrichment in one set of conditions (shortlist) and another
#'         # (longlist)
#'         short_v_long = leapR(geneset=ncipid, enrichment_method='enrichment_comparison', 
#'               eset=tset, primary_columns=shortlist, secondary_columns=longlist)
#'         
#'         # use enrichment_in_sets to calculate the most enriched pathways from the highest abundance proteins
#'         #     from one condition
#'         onept_sets = leapR(geneset=ncipid, enrichment_method='enrichment_in_sets',
#'                eset=tset, primary_columns="TCGA-13-1484", threshold=1.5)
#'                
#'          # use enrichment_in_order to calculate the most enriched pathways from the same condition
#'          #     Note: that this uses the entire set of abundance values and their order - whereas
#'          #     the previous example uses a hard threshold to get a short list of most abundant proteins
#'          #     and calculates enrichment based on set overlap. The results are likely to be similar - but
#'          #     with some notable differences.
#'          onept_order = leapR(geneset=ncipid, enrichment_method='enrichment_in_order',
#'                eset=tset, primary_columns="TCGA-13-1484")
#'                
#'          # use enrichment_in_pathway to calculate the most enriched pathways in a set of conditions
#'          #     based on abundance in the pathway members versus abundance in non-pathway members
#'          short_pathways = leapR(geneset=ncipid, enrichment_method='enrichment_in_pathway',
#'                eset=tset, primary_columns=shortlist)
#'                
#'          # use correlation_enrichment to calculate the most enriched pathways in correlation across
#'          #     the shortlist conditions
#'          short_correlation_pathways = leapR(geneset=ncipid, enrichment_method='correlation_enrichment',
#'                 eset=tset, primary_columns=shortlist)
#'
#'
#' @export
#'

leapR = function(geneset, enrichment_method, ...){
  .enrichment_wrapper(geneset, enrichment_method, ...)
}

.enrichment_wrapper = function(geneset, enrichment_method, eset=NULL, id_column=NULL, primary_columns=NULL,
                               secondary_columns=NULL, threshold=NULL, minsize=5, mode=NULL,
                               idmap=NA, fdr=0, min_p_threshold=NULL, sample_n=NULL,
                               enrichment_results=NA,
                               significance_threshold=NA, pathway_list=NA, background=NULL, targets=NULL,
                               subsample_components=NULL, ntimes=100, greaterthan=TRUE){

  # JEM: we will remove the 'tag' parameter. This is used to indicate different types of 'omics data and we'll
  #      need to handle that, but will do it differently

  #check that geneset is of correct class
  #if(!inherits(geneset, "geneset_data")) stop("geneset must be of class 'geneset_data")

 # if (!is.null(eset))
#    eset = Biobase::exprs(eset)
#  else
#    eset = NULL
  #check that enrichment_method is one of the designated methods
  if(!is.element(enrichment_method, c("correlation_enrichment", "correlation_comparison",
                                      "enrichment_in_relationships", 
                                      "enrichment_in_order", "enrichment_comparison", "enrichment_in_pathway",
                                      "enrichment_in_sets")))
    stop("enrichment_method must be one of the methods designated in the function documentation")
  

  
  #checking each enrichment method
  # JEM: changed here - moved these around so that the most important ones are first
  #      There are two different modes this can be run in and we can just have the user specify in the enrichment_method
  if(enrichment_method == "enrichment_in_order"){

    #message("'enrichment_by_ks' has the following optional arguments: eset=NULL, minsize=5, id_column=NULL, primary_columns=NULL")

    # eset will replace background - but we can do that at this level and 
    # not change the underlying code
    # primary_columns will replace abundance_column  - but again, at this level
    # minsize will remain the same
    # id_column will replace mapping_column
    # randomize might change - TBD

    # TODO: add a validation that primary_columns is one column id for this case
    if(!is.null(primary_columns)){
      if(length(primary_columns) > 1){stop("'primary_columns' must be a string of length 1")}
    }

    if(!is.null(primary_columns) & !is.null(eset)){
      
      ##primary_columsn can be in either the abundance data or the featureData, so here we check
      all_names <- c(colnames(eset), colnames(Biobase::fData(eset)))
      
      # TODO: add a validation that eset has a column named primary_columns
      if(!is.element(primary_columns, all_names)){stop("'primary_columns' must be column name of 'eset'")}
    }

    # JEM - notes for me: this application takes a background (matrix) and specifies which abundance column for calculating ks enrichment
    result = enrichment_in_groups(geneset=geneset, background=eset,
                                  method="ks", minsize=minsize, mapping_column=id_column,
                                  abundance_column=primary_columns)
  }

  else if(enrichment_method == "enrichment_in_sets"){

    #message("'enrichment_in_sets' has the following optional arguments: eset=NULL, minsize=5, id_column=NULL, primary_columns=NULL, greaterthan=TRUE, threshold=NA")


    # JEM - notes for me: this application takes a list of targets and a list of background genes and applies fisher's exact
    #     if we update this to take a eset and speciify a column then we need to also allow a threshold to be specified
    #     which is probably OK for this application -
    # eset will replace background
    # primary_columns will replace abundance column
    # minsize will remain the same
    # id_column will replace mapping_column
    # randomize isn't implemented for this application

    if(!is.null(primary_columns)){
      # We'll do the list splitting here
      # TODO: add a validation that the primary_columns is one column for this case
      if(typeof(primary_columns)=='list' && length(primary_columns) > 1){stop("must specifiy one column for 'primary_columns'")}
      else if(typeof(primary_columns)!="character"){stop("must specify one column for 'primary_columns'")}
    }

    if(!is.null(primary_columns) & !is.null(eset)){
      ##primary_columsn can be in either the abundance data or the featureData, so here we check
      all_names <- c(colnames(eset), colnames(Biobase::fData(eset)))
      
      # TODO: add a validation that the eset has a column with the primary_columns name
      if(!is.element(primary_columns, all_names)){stop("'primary_columns' must be a column name of 'eset'")}
    }
    
    # if we've been passed background instead of a eset then use it as the background
    if (is.null(background)) {
      background = eset #rownames(eset)
      backlist <- rownames(background)
    }else if(is.list(background)) {backlist  = background }
    else{ backlist <- rownames(background) }
    
    # TODO: we should allow the user to specify if they want it greater than or less than for the threshold
    if(greaterthan == FALSE & !is.null(threshold)){
      if(primary_columns%in%colnames(eset)) {
        targets = rownames(eset[which(Biobase::exprs(eset)[,primary_columns]<threshold),])
      }else{
        targets = rownames(Biobase::fData(eset)[which(Biobase::fData(eset)[,primary_columns]<threshold),])
      }
    }
    else if(greaterthan == TRUE & !is.null(threshold)){
      if(primary_columns%in%colnames(eset)) {
        targets = rownames(eset[which(Biobase::exprs(eset)[,primary_columns]>threshold),])
      }else{
        targets = rownames(Biobase::fData(eset)[which(Biobase::fData(eset)[,primary_columns]>threshold),])
      }
    }

    if(length(targets) == 0){stop("please adjust 'threshold' argument, there are no target genes captured by current threshold value")}

    if(!is.null(id_column) & !is.null(eset)){
      if(!is.element(id_column, colnames(Biobase::fData(eset)))){stop("'id_column' must be a column name of 'fData(eset)'")}
    }

    if(!is.null(id_column)) {
      background = unique(Biobase::fData(eset)[backlist,id_column])
      targets = unique(Biobase::fData(eset)[targets,id_column])
    }

    result = enrichment_in_groups(geneset=geneset, targets=targets, background=background,
                                  method="fishers", minsize=minsize)

  }

  else if(enrichment_method == "enrichment_in_pathway"){

    if(is.null(eset)){stop("'eset' argument is required")}
    
    #JEM removing these messages since they clutter
    #message("'enrichment_in_pathway' has the following optional arguments: id_column=NULL, primary_columns=NULL,fdr=0, secondary_columns=NULL,
    #                                      min_p_threshold=NULL, sample_n=NULL")

    result = enrichment_in_abundance(geneset=geneset, eset=eset, mapping_column=id_column,
                                     abundance_column=primary_columns, sample_comparison=NULL,
                                     fdr=fdr, min_p_threshold=min_p_threshold, sample_n=sample_n)

    # matchset=matchset, - this is a specific flag that's used for debugging
    # min_p_threshold is a parameter that will change the output - adding a couple of columns

    result = as.data.frame(result)
  }

  else if(enrichment_method == "enrichment_comparison"){
    
    if(is.null(eset)){stop("'eset' argument is required")}
    if(is.null(secondary_columns)){stop("'secondary_columns' argument is required")}
    if(is.null(primary_columns)){stop("'primary_columns' argument is required")}
    
    #JEM removing these messages since they clutter
    #message("'enrichment_comparison' has the following optional arguments: id_column=NULL, primary_columns=NULL,fdr=0, secondary_columns=NULL,
    #                                      min_p_threshold=NULL, sample_n=NULL")
    
    result = enrichment_in_abundance(geneset=geneset, eset=eset, mapping_column=id_column,
                                     abundance_column=primary_columns, sample_comparison=secondary_columns,
                                     fdr=fdr, min_p_threshold=min_p_threshold, sample_n=sample_n)
    
    # matchset=matchset, - this is a specific flag that's used for debugging
    # min_p_threshold is a parameter that will change the output - adding a couple of columns
    
    result = as.data.frame(result)
  }
  
  else if(enrichment_method == "enrichment_in_relationships"){

    if(is.null(eset)){stop("'eset' argument is required")}
    #message("'enrichment_in_relationships' has the following optional arguments: idmap=NA")

    #if(is.null(mode)){mode = 'original'} - this pertains to multiomics comparisons and we'll need to figure out how
    #                                       to recode this -

    result = enrichment_in_relationships(geneset=geneset, relationships=Biobase::exprs(eset), idmap=NULL)
  }

  # JEM- this function needs to be updated. Right now it's very simple and has the user define what set they
  #           want to use for correlation- but we should allow the user to specify a primary_columns parameter
  else if(enrichment_method == "correlation_enrichment"){

    if(is.null(eset)){stop("'eset' argument is required")}
    #message("'correlation_enrichment' has the following optional arguments:  id_column=NA, tag=NA")
    if(is.null(id_column)){id_column = NA}

    temp_result = correlation_enrichment(geneset = geneset, eset = eset, mapping_column = id_column)

    result = as.data.frame(temp_result$enrichment)
    attr(result, "corrmat") = temp_result$corrmat
  }
  
  else if (enrichment_method == "correlation_comparison") {
    if (is.null(eset)){stop("'eset' argument is required")}
    if(is.null(id_column)){id_column = NA}
    
    temp_result = correlation_comparison_enrichment(geneset=geneset, eset=eset, 
                                               set1=primary_columns, set2=secondary_columns,
                                               mapping_column=id_column)
    result = as.data.frame(temp_result)
  }
  
  return(result)
}

