#' enrichment_in_groups
#'
#' Calculate the enrichment in pathways using Fisher's exact or Kolmogorov-Smirnov test, using either
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
#' @param method method to use for statistical test, options are 'fishers' or 'ks'
#' @param minsize minimum size of set
#' @param mapping_column column name of mapping identifiers
#' @param abundance_column columns mapping abundance, either in the `assay` matrix or `rowData`
#' @param randomize true/false whether to randomize
#' @param silence_try_errors true/false to silence errors
#' @return data frame with enrichment results
#'
enrichment_in_groups <- function(geneset, targets = NULL, background = NULL, assay_name = NULL, method = "fishers", minsize = 5,
                                 mapping_column = NULL, abundance_column = NULL, randomize = FALSE,
                                 silence_try_errors = TRUE) {

  ##loop through every pathway for enrichment test
  results2 <- vapply(seq_along(1:length(geneset$names)), function(i){
    thisname <- geneset$names[i]
    thissize <- geneset$size[i]
    thisdesc <- geneset$desc[i]
    grouplist <- setdiff(geneset$matrix[i, 1:thissize], c("null",""))
    if (randomize) {
      # choose a random set of genes as this grouplist
      # A disadvantage is that we resample for each functional group rather than
      #   running one set of analyses on a fully scrambled set of functions.
      #   I don't think this should be a huge problem though.
      grouplist <- sample(unlist(geneset$matrix), length(grouplist))
    }

    # here backlist is actually a list of feature names
    if (length(background) == 1) { # klugey way to see if its an ExpressionSet
      if (!is.null(mapping_column)) {
        backlist <- SummarizedExperiment::rowData(background)[, mapping_column, drop = TRUE] |>
          unique()
      }
      backlist <- rownames(background)
    } else {
      backlist <- names(background)
    }

    in_back <- length(backlist)

    if (method == "fishers") {
      enr <- enrichment_by_fishers(targets, backlist, grouplist)
      p <- enr$fisher$p.value
      f <- enr$foldx
      mat <- enr$mat
      names <- enr$in_path_names

    #  res <- c(ingroup_n = mat[1, 1], ingroupnames = names, outgroup_n = mat[1, 2],
    #              ingroup_mean = mean(in_group, na.rm = TRUE), background_mean = mat[2,1], pvalue = p,
    #              oddratio = foldx, zscore = zscore)
      res <- c(ingroup_n = mat[1, 1], ingroupnames = names, ingroup_mean = NA,
               outgroup_n = mat[1, 2], outgroup_mean = NA, zscore = NA,
               oddsratio = f, pvalue = p,
               BH_pvalue = NA, SignedBH_pvalue = NA,
               background_n = mat[2,2], background_mean = mat[2,1])

    } else if (method == "ks") { # Kolmogorov-Smirnov test
      # in this case "background" must be the continuous variable from which grouplist can be drawn
      # backlist = background

      ## here backlist is a list of values, with namesa s feature names
      if (is.null(mapping_column)) { # use rownames here
        backlist <- rownames(SummarizedExperiment::assay(background, assay_name))
        in_group_name <- paste(intersect(grouplist, backlist), collapse = ", ")

        if (any(abundance_column %in% colnames(background))) { ## use the exprs
          in_group <- SummarizedExperiment::assay(background, assay_name)[grouplist[which(grouplist %in% backlist)], abundance_column, drop = FALSE]
          backlist <- SummarizedExperiment::assay(background, assay_name)[, abundance_column, drop = FALSE] # reassign backlist to values not strings
        } else { # mapping_column must be a featureData item
          in_group <- SummarizedExperiment::rowData(background)[grouplist[which(grouplist %in% backlist)], abundance_column, drop = FALSE]
          backlist <- SummarizedExperiment::rowData(background)[, abundance_column, drop = FALSE] # reassign backlist to values not strings
        }
      } else {
        # mapping_column adds the ability to use phospho-type data where the gene name (non-unique) is in the
        #       first column and the rownames are peptide ids
        # unfortunately this means that "background" has to be the whole matrix and abundance_column
        #       has to be specified, which is a bit ugly
        backlist <- SummarizedExperiment::rowData(background)[, mapping_column, drop = TRUE]

        in_group_name <- paste(intersect(backlist, grouplist), collapse = ", ")

        if (abundance_column %in% colnames(background)) { # we are using exprs fro values
          in_group <- SummarizedExperiment::assay(background, assay_name)[which(backlist %in% grouplist), abundance_column, drop = FALSE]
          backlist <- SummarizedExperiment::assay(background, assay_name)[, abundance_column, drop = FALSE]
        } else {
          in_group <- SummarizedExperiment::rowData(background)[which(backlist %in% grouplist), abundance_column, drop = FALSE]
          backlist <- SummarizedExperiment::rowData(background)[, abundance_column, drop = FALSE]
        }
      }

      in_path <- length(in_group)
      if ((in_path > minsize) & (any(!is.na(in_path))) & !all(is.na(in_group))) {
        in_back <- length(backlist)
        # The above block of code was replaced by the tryCatch block below to handle errors and warnings more elegantly.
        # The if(class(enr)) statement causes an error which doesn't let the rest of the code run.
        # Proposed change by Harkirat Sohi:
        enr <- NA
        enr <- tryCatch(
          {
            #suppressWarnings(
            ks.test(in_group, backlist)#)
          },
          error = function(e) {
            if (!silence_try_errors) message("An Error Occurred in KS test calculation")
            return(NA)
          }
        )
        if (length(enr) > 1) {
          p.value <- enr$p.value
        } else {
          p.value < -NA
        }
        # this expression of foldx might be subject to some weird pathological conditions
        # e.g. one sample has a background that is always negative, another that's positive
        # may pertain to zscore too (although not sure it should)
        # rank from largest to smallest
        # NOTE: By default, the function 'rank' outputs the position in the list from smallest to largest.
        # that is, rank(backlist) has the most negative values at the top, with the most positive at the bottom.
        # To get the positive entries at the top, negative at the bottom, we use rank(-backlist) instead.
        if (is.null(mapping_column)) {
          in_rank <- rank(-backlist)[which(rownames(background) %in% grouplist)]
        } else {
          in_rank <- rank(-backlist)[which(names(backlist) %in% grouplist)]
        }

        # this will give not a fold enrichment - but a score that ranges from 1 (most in top)
        #      to -1 (most in bottom).
        foldx <- 1 - ((mean(in_rank, na.rm = TRUE) / length(backlist)) / 0.5)

        zscore <- (mean(in_group, na.rm = TRUE) - mean(backlist, na.rm = TRUE)) / sd(in_group, na.rm = TRUE)
        res <- c(ingroup_n = in_path, ingroupnames = in_group_name, ingroup_mean = mean(in_group, na.rm = TRUE),
                 outgroup_n = in_back, outgroup_mean = mean(backlist, na.rm = TRUE), zscore = zscore,
                 oddsratio = foldx, pvalue = p.value,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)

      }else{
        res <- c(ingroup_n = in_path, ingroupnames = in_group_name, ingroup_mean = mean(in_group, na.rm = TRUE),
                 outgroup_n = in_back, outgroup_mean = mean(backlist, na.rm = TRUE), zscore = NA,
                 oddsratio = NA, pvalue = NA,
                 BH_pvalue = NA, SignedBH_pvalue = NA,
                 background_n = NA, background_mean = NA)
      }
    }
    return(res)

  }, c(ingroup_n = numeric(1), ingroupnames = "", ingroup_mean = numeric(1),
       outgroup_n = numeric(1), outgroup_mean = numeric(1), zscore = numeric(1),
       oddsratio = numeric(1), pvalue = numeric(1),
       BH_pvalue = numeric(1), SignedBH_pvalue = numeric(1),
       background_n = numeric(1), background_mean = numeric(1)))
  results <- t(results2)
  rownames(results) <- geneset$names
  results <- as.data.frame(results)

  return(results)
}
