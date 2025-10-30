#' enrichment_in_abundance
#'
#' Enrichment in abundance calculates enrichment in pathways by the difference
#' in abundance of the pathway members.
#' @importFrom stats sd
#' @importFrom stats p.adjust
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment SummarizedExperiment

#' @param geneset Gene set to calculate enrichment
#' @param eset Molecular abundance data in `SummarizedExperiment` format
#' @param assay_name Name of assay to compare
#' @param mapping_column Column to use to map identifiers
#' @param abundance_column  Columns to use to quantify abundance
#' @param fdr number of times to sample for FDR value
#' @param matchset Name of a set to use for enrichment
#' @param sample_comparison list of samples to use as comparison. if missing
#' background (eset) is used
#' @param min_p_threshold  Only include p-values lower than this
#' @param sample_n size of sample to use
#' @param silence_try_errors set to true to silence try errors
#' @return data frame of enrichment result

enrichment_in_abundance <- function(geneset,eset, assay_name,
            mapping_column = NULL, abundance_column = NULL, fdr = 0,
           matchset = NULL, sample_comparison = NULL,
           min_p_threshold = NULL, sample_n = NULL, silence_try_errors = TRUE) {

    #required input is the geneset, background and an abundance column
    stopifnot(is(geneset,'geneset_data'),
              is(eset,'SummarizedExperiment'),
              !is.null(abundance_column))

    if (!is.null(mapping_column)) {
      groupnames <- SummarizedExperiment::rowData(eset)[, mapping_column,
                                                        drop = TRUE]
    } else {
      groupnames <- rownames(eset)
    }

    if (!is.null(matchset)) {
      geneset <- geneset[which(geneset$names %in% matchset)]
    }
    results2 <- vapply(seq_len(length(geneset$names)), function(i){

      thisname <- geneset$names[i]
      thissize <- geneset$size[i]
      thisdesc <- geneset$desc[i]
      grouplist <- setdiff(geneset$matrix[i, seq_len(thissize),
                                          drop = FALSE],"")
      ingroupnames <- grouplist[which(grouplist %in% groupnames)]
      outgroupnames <- groupnames[which(!groupnames %in% grouplist)] |> unique()

      if (!is.null(sample_n)) {
        if (sample_n <= length(ingroupnames) ||
          sample_n <= length(outgroupnames)) {
          ingroupnames <- sample(ingroupnames, sample_n)
          outgroupnames <- sample(outgroupnames, sample_n)
        }
      }
      ingroup <- SummarizedExperiment::assay(eset, assay_name)[
        which(groupnames %in% ingroupnames),
        abundance_column[which(abundance_column %in% colnames(eset))]
      ]

      if (!is.null(sample_comparison)) {
        outgroup <- SummarizedExperiment::assay(eset, assay_name)[
          which(groupnames %in% ingroupnames),
          sample_comparison[which(sample_comparison %in% colnames(eset)),
                            drop = FALSE]
        ]
      } else {
        outgroup <- SummarizedExperiment::assay(eset, assay_name)
        outgroup <- outgroup[which(groupnames %in% outgroupnames),
                             abundance_column, drop = FALSE]
      }
      in_mean <- mean(unlist(ingroup), na.rm = TRUE)
      out_mean <- mean(unlist(outgroup), na.rm = TRUE)
      pvalue <- NA
      if (length(ingroup) > 1) {
        pvalue <- try(t.test(unlist(ingroup), unlist(outgroup))$p.value,
          silent =
            silence_try_errors
        )
        if (is(pvalue, "try-error")) {
          pvalue <- NA
        }
      }

      delta <- in_mean - out_mean
      if (fdr) {
        background <- c()
        abundances <- c(ingroup, outgroup)
        for (i in seq_len(fdr)) { #add this to unit tests when you can
          # randomly sample genes for fdr times
          ingroup <- sample(seq_len(length(abundances)), length(ingroup))
          outgroup <- which(!seq_len(length(abundances)) %in% ingroup)
          ingroup <- abundances[ingroup]
          outgroup <- abundances[outgroup]
          in_mean <- mean(ingroup, na.rm = TRUE)
          out_mean <- mean(outgroup, na.rm = TRUE)
          delta_r <- in_mean - out_mean
          background <- c(background, delta_r)
        }
        pvalue <- sum(abs(background) > abs(delta)) / length(background)
      }

      ingroupnames <- paste(ingroupnames, collapse = ", ")
      # question : do we want to calculate an oddsratio for this too?
      # answer: yes, but for now we'll use the mean of the ingroup compared
      # with the
      #        distribution of the background as a zscore
      zscore <- (out_mean - in_mean) / sd(unlist(outgroup), na.rm = TRUE)

      res <- c(ingroup_n = length(unlist(ingroup)),
               ingroupnames = ingroupnames, ingroup_mean = in_mean,
               outgroup_n = length(unlist(outgroup)),
               outgroup_mean = out_mean, zscore = zscore,
               oddsratio = delta, pvalue = pvalue,
               BH_pvalue = NA, SignedBH_pvalue = NA,
               background_n = NA, background_mean = NA)

      return(res)

    },c(ingroup_n = numeric(1), ingroupnames = "",
        ingroup_mean = numeric(1),
        outgroup_n = numeric(1), outgroup_mean = numeric(1),
        zscore = numeric(1),
        oddsratio = numeric(1), pvalue = numeric(1),
        BH_pvalue = numeric(1), SignedBH_pvalue = numeric(1),
        background_n = numeric(1), background_mean = numeric(1)))
    # update
    results <- t(results2)
    rownames(results) <- geneset$names
    results <- as.data.frame(results)

    if (!is.null(min_p_threshold)) {
      results <- results[results[,'pvalue'] < min_p_threshold, ]
      return(results)
    } else {
      return(results)
    }
  }
