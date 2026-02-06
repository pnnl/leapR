#' enrichment_in_groups
#'
#' Calculate the enrichment in pathways using Fisher's exact or
#' Kolmogorov-Smirnov test, using either the abundance column to identify
#' feature or the targets list. access through leapr wrapper
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats sd
#' @importFrom stats p.adjust
#'
#' @param geneset geneset to use for enrichment
#' @param targets targets to use for enrichment
#' @param background `SummarizedExperiment` describing background to use
#' @param assay_name is the name of the assay to use from the background
#' @param method method to use for statistical test, options are
#' 'fishers', 'ks', 'ztest', or 'chisq'. Remember that KS test assumes normality, so it would be good
#' to log your data before calling. NOTE: if you do not call `suppressWarnings` then
#' the KS test will warn you about ties.
#' @param minsize minimum size of set
#' @param mapping_column column name of mapping identifiers
#' @param abundance_column columns mapping abundance, either in the `assay`
#' matrix or `rowData`
#' @param randomize true/false whether to randomize
#' @param silence_try_errors true/false to silence errors
#' @return data frame with enrichment results
#'
enrichment_in_groups <- function(geneset,
                                 targets = c(),
                                 background = NULL,
                                 assay_name = NULL,
                                 method = "fishers",
                                 minsize = 5,
                                 mapping_column = NULL,
                                 log_transformed = FALSE,
                                 abundance_column = NULL,
                                 randomize = FALSE,
                                 silence_try_errors = TRUE) {

  #required input is the geneset, background and list of targets
  stopifnot(is(geneset,'geneset_data'),
            is(background,'SummarizedExperiment'))

  ##loop through every pathway for enrichment test
  results2 <- vapply(seq_len(length(geneset$names)), function(i){
    thisname <- geneset$names[i]
    thissize <- geneset$size[i]
    thisdesc <- geneset$desc[i]
    grouplist <- setdiff(geneset$matrix[i, seq_len(thissize)], c("null",""))
    if (randomize) {
      # choose a random set of genes as this grouplist
      # A disadvantage is that we resample for each functional group rather than
      #   running one set of analyses on a fully scrambled set of functions.
      grouplist <- sample(unlist(geneset$matrix), length(grouplist))
    }

    # grouplist and backlist
    # are both list of feature names of foreground/background
     if (!is.null(mapping_column)) {
        backlist <- SummarizedExperiment::rowData(background)[, mapping_column,
                                                                drop = TRUE] |>
        unique()
     }else{
      backlist <- rownames(background)
     }
    backlist <- unlist(backlist)
    in_back <- length(backlist)

    if (method == "fishers") {
      stopifnot(length(targets) > 0)
      enr <- enrichment_by_fishers(targets, backlist, grouplist)
      p <- enr$fisher$p.value
      f <- enr$foldx
      mat <- enr$mat
      names <- enr$in_path_names
      res <- c(ingroup_n = mat[1, 1], ingroupnames = names, ingroup_mean = NA,
               outgroup_n = mat[1, 2], outgroup_mean = NA, zscore = NA,
               oddsratio = f, pvalue = p,BH_pvalue = NA, SignedBH_pvalue = NA,
               background_n = mat[2,2], background_mean = mat[2,1])

    } else if (method %in% c('ks','ztest','chisq')) { ##try out one of our 
      #three continuous tests
       # in this case "background" must be the continuous variable
       #lets ensure that the data are normally distributed

      #these are the indices of the background
      group_ind <- which(backlist %in% grouplist)
      in_group_name <- paste(intersect(backlist, grouplist), collapse = ", ")

      if (any(abundance_column %in% colnames(background))) {
       #  in_group <- SummarizedExperiment::assay(background)[group_ind,
      #                                                       abundance_column,
      #                                                       drop = TRUE]
         backvals <- SummarizedExperiment::assay(background)[,
                                                             abundance_column,
                                       drop = TRUE]
      } else { ##we are working with rowData
       # in_group <- SummarizedExperiment::rowData(background)[group_ind,
      #                                                        abundance_column,
      #                                drop = TRUE]
        backvals <- SummarizedExperiment::rowData(background)[,
                                                              abundance_column,
                                      drop = TRUE]
      }

      ##make normal - this is new and might change results
     backvals <- (backvals - mean(backvals, na.rm=T))/sd(backvals,na.rm=T)
      
      in_group <- backvals[group_ind]
      ##remove NA vals from in group
      in_group <- in_group[!is.na(in_group)]
      in_group_mean <- mean(in_group)
      
      names(backvals) <- backlist#[-group_ind]
      in_back <- length(backvals)
      
      outgroup_mean = mean(backvals[-group_ind], na.rm = T)
      
      in_path <- length(in_group) #how many left after na.rm
      if ((in_path > minsize) & (any(!is.na(in_path))) &
         !all(is.na(in_group))) {
        if (method == 'ks') { 
              in_back <- length(backlist)
              enr <- NA
              enr <- tryCatch(
          {
            #suppressWarnings(
            ks.test(in_group, backvals[-group_ind], exact = FALSE)#)
          })
          if (length(enr) > 1) {
            p.value <- enr$p.value
          } else {
            p.value < -NA
          }
          # this expression of foldx might be subject to some weird pathological
          # conditions e.g. one sample has a background that is always negative,
          # another that's positive may pertain to zscore too (although not sure
          # it should) rank from largest to smallest
          if (is.null(mapping_column)) {
            in_rank <- rank(-backvals)[which(rownames(background) %in%
                    #  backvals[which(rownames(background) %in%                  
                                       grouplist)]
          } else {
            in_rank <- rank(-backvals)[which(names(backvals) %in% grouplist)]
          }
          # this will give not a fold enrichment - but a score that ranges
          # from 1 (most in top)
          #      to -1 (most in bottom).
          foldx <- 1 - ((mean(in_rank) / length(backvals)) / 0.5)
        } else if (method == 'ztest') {
        ##add z test here
            foldx <- sqrt(thissize) * in_group_mean #z test
            p.value <- 2 * pnorm(-abs(foldx)) #2 sided pvalue
        } else if (method == 'chisq') {
         #add chisq test here
          foldx <- vapply(in_group, function(x) (x - in_group_mean)^2,numeric(1))
          foldx <- (sum(foldx) - (thissize - 1))/(2*(thissize - 1))
          p.value <- 2 * pnorm(-abs(foldx)) # 2 sided pvalue
        } 
        zscore <- (in_group_mean - outgroup_mean)
        zscore <- zscore / sd(in_group)

       
      res <- c(ingroup_n = in_path, ingroupnames = in_group_name,
               ingroup_mean = in_group_mean,
               outgroup_n = in_back, outgroup_mean = outgroup_mean,
               zscore = zscore,
               oddsratio = foldx, pvalue = p.value,
               BH_pvalue = NA, SignedBH_pvalue = NA,
               background_n = NA, background_mean = NA)
      }  else{
        res <- c(ingroup_n = in_path, ingroupnames = in_group_name,
                 ingroup_mean = in_group_mean,
                 outgroup_n = in_back, outgroup_mean = outgroup_mean,
                 zscore = NA,
                 oddsratio = NA, pvalue = NA,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)
      } 
    } 
    else {
        res <- c(ingroup_n = in_path, ingroupnames = in_group_name,
                   ingroup_mean = NA,
                   outgroup_n = in_back, outgroup_mean = mean(backvals,
                                                              na.rm = TRUE),
                   zscore = NA,
                 oddsratio = NA, pvalue = NA,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)
      }
  
    
    return(res)
  
  }, c(ingroup_n = numeric(1), ingroupnames = "", ingroup_mean = numeric(1),
       outgroup_n = numeric(1), outgroup_mean = numeric(1),
       zscore = numeric(1),
       oddsratio = numeric(1), pvalue = numeric(1),
       BH_pvalue = numeric(1), SignedBH_pvalue = numeric(1),
       background_n = numeric(1), background_mean = numeric(1)))
  results <- t(results2)
  rownames(results) <- geneset$names
  results <- as.data.frame(results)

  return(results)
}
