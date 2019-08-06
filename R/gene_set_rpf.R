#'gene_set_rpf
#'
#'
#'

gene_set_rpf <- function(data, survivaldata, timevar, statusvar, geneset,
                         load.saved.file=T, only.load.files=F, gsprefix="",
                         min.geneset.size=6, max.steps=30, num.repeats=0,
                         pre.pvalue.threshold=NULL, validation.data=NULL, make.validata=T,
                         collect.validation.pvalues=F, p.value.threshold=0.05, randomize.order=F,
                         randomize.genes=F, match.set.size=F, training.set.ratio=0.5,
                         extra.features=NULL, generator.mode=F, save.cv.pvalues=F,
                         translation_table=NULL, signature.level=20, sig.vote=F, ...) {
  # for each gene set construct a separate rpf model

  results = list()
  datas = list()
  raw_cv_pvalues = list()
  pvalues = matrix(nrow=length(geneset$names), ncol=max.steps)
  patientids = colnames(data)[which(colnames(data) %in% rownames(survivaldata))]

  if (!is.null(validation.data)) {
    pvalues_valid = matrix(nrow=length(geneset$names), ncol=max.steps)
    rownames(pvalues_valid) = geneset$names

    collect_pvalues_training = matrix(ncol=max.steps)
    collect_pvalues_validation = matrix(ncol=max.steps)

    model_stats = matrix(nrow=length(geneset$names), ncol=5)
    rownames(model_stats) = geneset$names
    colnames(model_stats) = c("Significant-Testing", "Significant-Validation",
                              "Significant-Matching", "Number-Models-Testing",
                              "Number-Models-Validation")

    validata = validation.data
    # this assumes that the validation data can be referenced from the same survival data
    # FIXME
    if (make.validata) validdata = rpf_coxph_formatdata(validation.data, survivaldata, timevar, statusvar,
                                                        extra.features=extra.features)
  }
  rownames(pvalues) = geneset$names

  runorder = 1:length(geneset$names)

  if (randomize.order) runorder = sample(1:length(geneset$names))

  for (i in runorder) {
    thisname = geneset$names[i]
    thissize = geneset$size[i]
    thisdesc = geneset$desc[i]
    grouplist = geneset$matrix[i,1:thissize]

    cat("\n\n\n\n", thisname, "\n")
    pvalues_vx = NULL

    # implement a translator that can take many to one relationships (e.g. phosphopeptides to gene symbols in pathway)
    #           and return the correct identifiers
    # This should be a matrix/data frame with the pathway identifiers (e.g. gene symbols) in the first column
    #         and the data identifiers (e.g. phosphopeptide ids) as rownames, which should be unique
    # NOTE: this could easily be used for incorporation through interactions like PPIs!!!
    if (!is.null(translation_table)) {
      translated_group = rownames(translation_table)[which(translation_table[,1] %in% grouplist)]
      grouplist = c(grouplist, translated_group)
    }

    if (match.set.size && randomize.genes) x = length(grouplist)
    else x = length(which(rownames(data) %in% grouplist))

    # Pull a random set of genes from the experimental data
    if (randomize.genes) {
      grouplist = sample(rownames(data), x)
      cat("Randomizing gene list...\n")
    }

    if (x <= min.geneset.size) {
      # add a blank for this pathway since it doesn't
      #     make sense to fit a model
      results = c(results, list(NULL))
      datas = c(datas, list(NULL))
      raw_cv_pvalues = c(raw_cv_pvalues, list(NULL))
      result = NULL
    }
    else {
      fprefix = paste(gsprefix, thisname, sep="")
      fname = paste(fprefix, "_rpfmodel.data", sep="")
      result = NULL
      formatdata = NULL
      if (file.exists(fname) && load.saved.file) {
        cat("Loading saved file", fname, "\n")
        load(fname)
        # sanity check on file type- this checks to see if the
        #     saved file is from a single run- which can happen
        #     when running multiple times if you forget your num.repeats
        #     on one of the runs and causes an ugly problem
        if (num.repeats > 0 && "call" %in% names(result)) {
          cat("Bad run for cross-validation: ", thisname, "\n")
          result = NULL
          formatdata = NULL
        }
      }
      if (is.null(result) &! only.load.files) {
        pfname = paste(gsprefix, "pre.", thisname, "_rpfmodel.data", sep="")
        cat("Fitting", thisname, "with", x, "variables\n")

        setdata = data[which(rownames(data) %in% grouplist),]


        formatdata = rpf_coxph_formatdata(setdata, survivaldata, timevar, statusvar,
                                          extra.features=extra.features)
        cat(dim(formatdata[[1]]$extra.features), extra.features, "\n")
        save(setdata, formatdata, file=pfname)

        result = list()

        groupids = grouplist[which(grouplist %in% rownames(data))]

        if (num.repeats > 0) {
          cat("...running CV for", thisname, "\n")
          if (is.null(pre.pvalue.threshold))
            result = run.cv(formatdata[[1]], formatdata[[2]], max.steps=max.steps,
                            num.repeats=num.repeats, prefix=fprefix,
                            training.set.ratio=training.set.ratio, ...)
          else {
            for (c in 1:num.repeats) {
              training_samples = sample(patientids,
                                        training.set.ratio*ncol(data))
              testing_samples = patientids[which(!patientids %in% training_samples)]

              thisformatdata = as.data.frame(cbind(formatdata[[1]]$time[which(rownames(formatdata[[1]]$x) %in% training_samples)],
                                                   formatdata[[1]]$status[which(rownames(formatdata[[1]]$x) %in% training_samples)],
                                                   formatdata[[1]]$x[training_samples, groupids]))
              colnames(thisformatdata)[1] = timevar
              colnames(thisformatdata)[2] = statusvar

              plist = coxph_expression(thisformatdata,
                                       statusvar, timevar, 3)

              siggroupids = sample(rownames(plist)[which(plist[,2]<pre.pvalue.threshold)])

              trainformatdata = rpf_coxph_formatdata(data[siggroupids, training_samples],
                                                     survivaldata, timevar, statusvar,
                                                     extra.features=extra.features)

              cat("...selecting", length(siggroupids), "significant components for use in model\n")

              res = itemset.coxpath(trainformatdata[[1]], trainformatdata[[2]],
                                    max.steps=max.steps, ...)
              result[[length(result)+1]] = list(fit=res, training.set=patientids %in% training_samples)
            }
          }
        }
        else {
          result = itemset.coxpath(formatdata[[1]], formatdata[[2]], max.steps=max.steps, ...)
        }

        save(formatdata, result, file=fname)
      }
      if (!is.null(result)) {
        if (num.repeats > 0) {
          pvalues_run = NULL
          pvaluefname = paste(gsprefix, thisname, "_rpfmodelcvpvalues.data", sep="")

          overlap = colnames(data)[which(colnames(data) %in% rownames(survivaldata))]
          signatures = matrix(ncol=length(overlap), nrow=num.repeats)
          colnames(signatures) = overlap

          if (file.exists(pvaluefname)) {
            load(pvaluefname)
            cat("Loading saved CV pvalues file", pvaluefname, "\n")
            if (!is.null(validation.data)) {
              pvalues_run = combined$testing
              pvalues_vrun = combined$validation
            }
          }
          else if (!is.null(result)) {
            cat("Calculating pvalue matrix for ", thisname, "\n")
            if (!is.null(validation.data)) cat("Calculating validation pvalues for ", thisname, "\n")
            pvalues_run = matrix(nrow=num.repeats, ncol=max.steps)
            pvalues_vrun = matrix(nrow=num.repeats, ncol=max.steps)
            for (j in 1:num.repeats) {
              pvaluex = sapply(1:max.steps, function (l) {
                r = try(plot_model_kaplan_meier(result[[j]][[1]], formatdata[[1]],
                                                level=l, training.set=result[[j]][[2]], plotthis=F)$testing, silent=T)
                #r = try(plot_model_kaplan_meier(result[[j]][[1]], formatdata[[1]],
                #                                level=l, training.set=result[[j]][[2]], plotthis=F)$training, silent=T)
                if (class(r) == "try-error") r = NA
                return(r)
              })
              # this is the vector of pvalues for each cross-validation run
              pvalues_run[j,] = pvaluex

              #create a signature at the specified level to incorporate into an integrated prediction
              cat("Creating signature", j,"at level", signature.level, "\n")

              testing.set = overlap[which(!result[[j]][[2]])]

              siglevel = min(signature.level, length(result[[j]][[1]]$b.predictor))
              #browser()
              signatures[j,testing.set] = build_predictor(result[[j]][[1]], formatdata[[1]], level=siglevel, browser=F)[testing.set]

              if (!is.null(validation.data)) {
                pvalue_vx = sapply(1:max.steps, function (l) {
                  r = try(plot_model_kaplan_meier(result[[j]][[1]], validdata[[1]], level=l, plotthis=F), silent=T)
                  if (class(r) == "try-error") {cat ("error\n"); r = NA}
                  return(r)
                })
                #pvalues_valid[thisname,] = pvalue_vx
                pvalues_vrun[j,] = pvalue_vx
              }
            }
            if (!is.null(validation.data)) {
              combined = list(testing=pvalues_run, validation=pvalues_vrun)
              save(combined, file=pvaluefname)
            }
            else save(pvalues_run, file=pvaluefname)
            # calculate integrated p.value from signatures
            if (sig.vote)
            {
              signatures_rank = t(apply(signatures, 1, rank))
              signatures_rank[which(is.na(signatures))] = NA
              signatures = signatures_rank
            }
            signatures_int = colMeans(signatures, na.rm=T)

            intsigp = try(plot_model_kaplan_meier(1, list(time=formatdata[[1]][,"time"],
                                                          status=formatdata[[1]][,"status"]), 1,
                                                  prediction=signatures_int, plotthis=F))
            if (!is.numeric(intsigp)) f=0
            cat("Integrated p-value:", intsigp, "\n")
            #browser()
          }
          else {
            pvalues_run = NULL
          }
          # now we get the lowest mean pvalue- note that averaging pvalues may not
          #     be that reasonable, but I think it should be OK for ranking
          #     results

          # one issue is that there are places where there are missing values since
          #     the runs never got to that stage so we need to filter out the one-
          #     and two- 'hit'
          #cat("...Pvalues for ", thisname, "\n")

          #pvaluex = 10^sapply(1:max.steps, function (l) {
          #    if (sum(!is.na(pvalues_run[,l]))<3) return(NA)
          #    mean(log(pvalues_run[,l], 10), na.rm=T)
          #    })
          # The problem with using the mean is that it can be skewed by single
          #     very low pvalues. Better to use the median
          pvaluex = 10^sapply(1:max.steps, function (l) {
            if (sum(!is.na(pvalues_run[,l]))<3) return(NA)
            median(log(pvalues_run[,l], 10), na.rm=T)
          })
          pvalues[thisname,] = pvaluex

          if (!is.null(validation.data)) {
            pvaluevx = 10^sapply(1:max.steps, function (l) {
              if (sum(!is.na(pvalues_vrun[,l]))<3) return(NA)
              median(log(pvalues_vrun[,l], 10), na.rm=T)
            })
            pvalues_valid[thisname,] = pvaluevx

            model_stats[thisname,] = c(sum(pvalues_run<p.value.threshold, na.rm=T),
                                       sum(pvalues_vrun<p.value.threshold, na.rm=T),
                                       sum(pvalues_vrun[which(pvalues_run<p.value.threshold)]<p.value.threshold, na.rm=T),
                                       sum(!is.na(pvalues_run)), sum(!is.na(pvalues_vrun)))

            if (collect.validation.pvalues) {
              collect_pvalues_training = rbind(collect_pvalues_training, pvalues_run)
              collect_pvalues_validation = rbind(collect_pvalues_validation, pvalues_vrun)
            }
          }
        }
        else {
          # we just need to extract a list of the pvalues
          pvaluex = sapply(1:max.steps, function (l) {
            r = try(plot_model_kaplan_meier(result, formatdata[[1]], level=l, plotthis=F), silent=T)
            if (class(r) == "try-error") r = NA
            return(r)
          })
          pvalues[thisname,] = pvaluex
          if (!is.null(validation.data)) {
            pvalue_vx = sapply(1:max.steps, function (l) {
              r = try(plot_model_kaplan_meier(result, validdata[[1]], level=l, plotthis=F), silent=T)
              if (class(r) == "try-error") {cat(r, "\n"); r = NA}
              return(r)
            })
            pvalues_valid[thisname,] = pvalue_vx
          }
        }
      }
      # in generator mode we're just generating the results files- not
      # storing them in memory for the final analysis
      if (generator.mode) {
        formatdata = NULL
        result = NULL
        pvalues_run = NULL
      }
      if (save.cv.pvalues == TRUE) {
        #save the individual pvalues in a data structure
        raw_cv_pvalues = c(raw_cv_pvalues, list(pvalues_run))
      }
      datas = c(datas, list(formatdata))
      results = c(results, list(result))
    }
  }
  gs = geneset$names[1]
  names(results) = gs
  names(datas) = gs

  if (save.cv.pvalues == TRUE) names(raw_cv_pvalues) = geneset$names

  bestmodel = sapply(rownames(pvalues), function (r) which.min(pvalues[r,]))
  bestp = unlist(sapply(rownames(pvalues), function (r)
    pvalues[r,which.min(pvalues[r,])]))
  names(bestp) = sapply(names(bestp), function (n) strsplit(n, "\\.")[[1]][1])
  bestp = cbind(names(bestp), sapply(names(bestp), function (i) geneset$desc[which(geneset$names==i)]), bestp)

  if (is.null(validation.data)) return(list(fits=results, data=datas,
                                            pvalues=pvalues, bestp=bestp,
                                            bestmodels=bestmodel, raw_cv_pvalues=raw_cv_pvalues))
  return(list(fits=results, data=datas, pvalues=pvalues, bestp=bestp, bestmodels=bestmodel,
              pvalues_valid=pvalues_valid, model_stats=model_stats,
              training_matrix=collect_pvalues_training, validation_matrix=collect_pvalues_validation))
}
