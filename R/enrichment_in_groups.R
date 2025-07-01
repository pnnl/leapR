#' enrichment_in_groups
#'
#' Calculate the enrichment in pathways using Fisher's exact or Kolgmorov-Smirnov test, using either
#' the primary_columns to identify feature or the targets list.
#' # access through leapr wrapper
#' @import Biobase
#' @param geneset geneset to use for enrichment
#' @param targets targets to use for enrichmenet
#' @param background `ExpressionSet` describing background to use
#' @param method method to use for statistical test, options are 'fishers' or 'ks'
#' @param minsize minimum size of set
#' @param mapping_column column name of mapping identifiers
#' @param abundance_column columns mapping abundance, either in the `exprs` matrix or `featureData`
#' @param randomize true/false whether to randomize
#' @param silence_try_errors true/false to silence errors
#' @import stats
#' @return data frame with enrichment results
#' 
enrichment_in_groups <- function(geneset, targets=NULL, background=NULL, method="fishers", minsize=5,
                                 mapping_column=NULL, abundance_column=NULL, randomize=FALSE,
                                 silence_try_errors=TRUE) {
  
  resultp = c()
  resultf = c()
  results = data.frame(row.names = geneset$names,
                       ingroup_n=rep(NA_real_, length(geneset$names)), ingroupnames=rep(NA_character_, length(geneset$names)), 
                       ingroup_mean=rep(NA_real_, length(geneset$names)), outgroup_n=rep(NA_real_, length(geneset$names)), 
                       outgroup_mean=rep(NA_real_, length(geneset$names)), score=rep(NA_real_, length(geneset$names)), oddsratio=rep(NA_real_, length(geneset$names)), 
                       pvalue=rep(NA_real_, length(geneset$names)), BH_pvalue=rep(NA_real_, length(geneset$names)), 
                       SignedBH_pvalue=rep(NA_real_, length(geneset$names)), background_n=rep(NA_real_, length(geneset$names)),
                       background_mean=rep(NA_real_, length(geneset$names)), stringsAsFactors = FALSE)
  
  
  for (i in 1:length(geneset$names)) {
    thisname = geneset$names[i]
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    grouplist = setdiff(geneset$matrix[i,1:thissize],'null')
    if (randomize) {
      # choose a random set of genes as this grouplist
      # A disadvantage is that we resample for each functional group rather than
      #   running one set of analyses on a fully scrambled set of functions.
      #   I don't think this should be a huge problem though.
      grouplist = sample(unlist(geneset$matrix), length(grouplist))
    }
    
    #here backlist is actually a list of feature names
    #todo fix this!!
    if (length(background) == 1) { #klugey way to see if its an ExpressionSet
      if (!is.null(mapping_column)) {
        backlist <- Biobase::fData(background)[,mapping_column] |>
          unique()
      }
      backlist <- rownames(background)
    }else{
      backlist <- names(background)
    }
    in_back = length(backlist)
    
    if (method == "fishers") {
      
      enr = enrichment_by_fishers(targets, backlist, grouplist)
      p = enr$fisher$p.value
      f = enr$foldx
      mat = enr$mat
      names = enr$in_path_names
      
      results[thisname,"ingroup_n"] = mat[1,1]
      results[thisname,"ingroupnames"] = names
      results[thisname,"outgroup_n"] = mat[1,2]
      results[thisname,"pvalue"] = p
      results[thisname,"oddsratio"] = f
      results[thisname,"background_n"] = mat[2,2]
      # putting this here is a bit of a kludge (since it's not actually a mean, it's a count)
      results[thisname,"background_mean"] = mat[2,1]
    }
    else if (method == "ks") {    #Kolmogorov-Smirnov test
      # in this case "background" must be the continuous variable from which grouplist can be drawn
      #backlist = background
      
      ##here backlist is a list of values, with namesa s feature names
      
      if (is.null(mapping_column)) { #use rownames here
        in_group_name = paste(intersect(grouplist, rownames(background)), collapse = ", ")
        
        if(abundance_column %in% colnames(background)){ ##use the exprs
          in_group = Biobase::exprs(background)[grouplist[which(grouplist %in% rownames(background))],abundance_column]
          backlist = Biobase::exprs(background)[,abundance_column]
        }else{ #mapping_column must be a featureData item
          in_group = Biobase::fData(background)[grouplist[which(grouplist %in% rownames(background))],abundance_column]
          backlist = Biobase::fData(background)[,abundance_column]
        }
      }
      else {
        # mapping_column adds the ability to use phospho-type data where the gene name (non-unique) is in the
        #       first column and the rownames are peptide ids
        # unfortunately this means that "background" has to be the whole matrix and abundance_column
        #       has to be specified, which is a bit ugly
        in_group_name = paste(intersect(backlist, grouplist), collapse = ", ")
        if (abundance_column %in% colnames(background)) { #we are using exprs fro values
          in_group = Biobase::exprs(background)[which(Biobase::fData(background)[,mapping_column] %in% grouplist),abundance_column]
          backlist = Biobase::exprs(background)[,abundance_column]
        }else{
          in_group = Biobase::fData(background)[which(Biobase::fData(background)[,mapping_column] %in% grouplist),abundance_column]
          backlist = Biobase::fData(background)[,abundance_column]          
        }
      }
      
      in_path = length(in_group)

      if ((in_path > minsize) & (any(!is.na(in_path)))) {
        in_back = length(backlist)

        # enr = try(ks.test(in_group, backlist), silent=silence_try_errors)
        # if (class(enr) == "try-error") {
        #   enr = NA
        #   p.value = NA
        #   #browser()
        # }
        # else {
        #   p.value = enr$p.value
        # }
        # The above block of code was replaced by the tryCatch block below to handle errors and warnings more elegantly.
        # The if(class(enr)) statement causes an error which doesn't let the rest of the code run.
        # Proposed change by Harkirat Sohi:
        enr <- NA
        enr <- tryCatch(
          {
            suppressWarnings(ks.test(in_group, backlist))
          },
          error=function(e) {
            if (!silence_try_errors) message('An Error Occurred in KS test calculation')
            return(NA)
          }
        )
        if(length(enr)>1)
          p.value <- enr$p.value
        else
          p.value < -NA
     
        
        # this expression of foldx might be subject to some weird pathological conditions
        # e.g. one sample has a background that is always negative, another that's positive
        # may pertain to zscore too (although not sure it should)
        #foldx = mean(in_group, na.rm=T)/mean(background, na.rm=T)
        
        # rank from largest to smallest
        # NOTE: By default, the function 'rank' outputs the position in the list from smallest to largest.
        # that is, rank(backlist) has the most negative values at the top, with the most positive at the bottom.
        # To get the positive entries at the top, negative at the bottom, we use rank(-backlist) instead.
        if (is.null(mapping_column)) in_rank = rank(-backlist)[which(rownames(background) %in% grouplist)]
        else in_rank = rank(-backlist)[which(names(backlist) %in% grouplist)]
        
        # this will give not a fold enrichment - but a score that ranges from 1 (most in top)
        #      to -1 (most in bottom).
        foldx = 1 - ((mean(in_rank, na.rm = TRUE)/length(backlist))/0.5)
        
        zscore = (mean(in_group, na.rm = TRUE) - mean(backlist, na.rm = TRUE))/sd(in_group, na.rm = TRUE)
        #browser()
        
        results[thisname,"ingroup_n"] = in_path
        results[thisname,"ingroupnames"] = in_group_name
        results[thisname,"ingroup_mean"] = mean(in_group, na.rm = TRUE)
        results[thisname,"outgroup_n"] = in_back
        results[thisname,"outgroup_mean"] = mean(backlist, na.rm = TRUE)
        results[thisname,"zscore"] = zscore
        results[thisname,"pvalue"] = p.value
        results[thisname,"oddsratio"] = foldx
      }
    }
  }
  results[,"BH_pvalue"] = p.adjust(results[,"pvalue"], method="BH")
  results[,"SignedBH_pvalue"] = results[,"BH_pvalue"]*sign(results[,"ingroup_mean"])
  return(results)
}
