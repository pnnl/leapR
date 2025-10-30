#' enrichment_in_groups
#'
#' Calculate the enrichment in pathways using Fisher's exact or
#' Kolmogorov-Smirnov test, using either
#' the primary_columns to identify feature or the targets list.
#' # access through leapr wrapper
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
#' 'fishers' or 'ks'
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

    # here backlist will be a list of feature names
     if (!is.null(mapping_column)) {
        backlist <- SummarizedExperiment::rowData(background)[, mapping_column,
                                                                drop = TRUE] |>
        unique()
     }else{
      backlist <- rownames(background)
     }

    in_back <- length(backlist)

    if (method == "fishers") {
      stopifnot(length(targets)>0)
      enr <- enrichment_by_fishers(targets, backlist, grouplist)
      p <- enr$fisher$p.value
      f <- enr$foldx
      mat <- enr$mat
      names <- enr$in_path_names
      res <- c(ingroup_n = mat[1, 1], ingroupnames = names, ingroup_mean = NA,
               outgroup_n = mat[1, 2], outgroup_mean = NA, zscore = NA,
               oddsratio = f, pvalue = p,BH_pvalue = NA, SignedBH_pvalue = NA,
               background_n = mat[2,2], background_mean = mat[2,1])

    } else if (method == "ks") { # Kolmogorov-Smirnov test
         # in this case "background" must be the continuous variable
      ## here backlist is a list of values, with namesa s feature names
      if (is.null(mapping_column)) { # use rownames here
          backlist <- rownames(SummarizedExperiment::assay(background,
                                                           assay_name))
        if (any(abundance_column %in% colnames(background))) { ## use the exprs
            in_group <- SummarizedExperiment::assay(background,
                                                    assay_name)
            backlist <- SummarizedExperiment::assay(background,assay_name)
        } else { # mapping_column must be a featureData item
            in_group <- SummarizedExperiment::rowData(background)
            backlist <- SummarizedExperiment::rowData(background)

        }
        # reassign backlist to values not strings
      } else {
        # mapping_column adds the ability to use phospho-type data
        # where the gene name (non-unique) is in the
        backlist <- SummarizedExperiment::rowData(background)
        backlist <- backlist[, mapping_column, drop = TRUE]

        if (abundance_column %in% colnames(background)) { s
            in_group <- SummarizedExperiment::assay(background, assay_name)
            backlist <- SummarizedExperiment::assay(background, assay_name)
        } else {
          in_group <- SummarizedExperiment::rowData(background)
          backlist <- SummarizedExperiment::rowData(background)
        }
      }

      in_group <- in_group[grouplist[which(grouplist %in% backlist)],
                           abundance_column, drop = FALSE]
      backlist <- backlist[,abundance_column, drop = FALSE]
      in_group_name <- paste(intersect(backlist, grouplist), collapse = ", ")

      in_path <- length(in_group)
      if ((in_path > minsize) & (any(!is.na(in_path))) &
          !all(is.na(in_group))) {
        in_back <- length(backlist)
        enr <- NA
        enr <- tryCatch(
          {
            #suppressWarnings(
            ks.test(in_group, backlist)#)
          },
          error = function(e) {
            return(NA)
          }
        )
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
          in_rank <- rank(-backlist)[which(rownames(background) %in%
                                             grouplist)]
        } else {
          in_rank <- rank(-backlist)[which(names(backlist) %in% grouplist)]
        }

        # this will give not a fold enrichment - but a score that ranges
        # from 1 (most in top)
        #      to -1 (most in bottom).
        foldx <- 1 - ((mean(in_rank, na.rm = TRUE) / length(backlist)) / 0.5)

        zscore <- (mean(in_group, na.rm = TRUE) - mean(backlist, na.rm = TRUE))
        zscore <- zscore / sd(in_group, na.rm = TRUE)
        res <- c(ingroup_n = in_path, ingroupnames = in_group_name,
                 ingroup_mean = mean(in_group, na.rm = TRUE),
                 outgroup_n = in_back, outgroup_mean = mean(backlist,
                                                            na.rm = TRUE),
                 zscore = zscore,
                 oddsratio = foldx, pvalue = p.value,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)

      }else{
          res <- c(ingroup_n = in_path, ingroupnames = in_group_name,
                   ingroup_mean = mean(in_group, na.rm = TRUE),
                   outgroup_n = in_back, outgroup_mean = mean(backlist,
                                                              na.rm = TRUE),
                   zscore = NA,
                 oddsratio = NA, pvalue = NA,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)
      }
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
