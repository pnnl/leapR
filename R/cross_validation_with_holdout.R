#'cross_validation_with_holdout
#'
#'
#'

cross_validation_with_holdout = function(data, survivaldata, timevar, statusvar,
                                         consensus_holdout=0.5, validation_holdout=0.3333,
                                         validation_dataset=NA,
                                         foldx=5, consensus_foldx=25, consensus_maxsteps=100,
                                         consensus_level=40, consensus_count_filter=4,
                                         load.saved.file=T, ...) {
  # controller function to take a (largish) dataset and build a consensus
  # signature using a cross-validation approach _on_a_subset_of_the_total_data
  # then apply the consensus to the held out set as a more appropriate validation

  results = list()

  if (!is.na(validation_dataset)) {
    # we can supply a different validation dataset to serve here
    #   in which case we don't hold anything out from the consensus building CV
    validation_holdout = 0
  }

  namebase = paste("hocv.", as.character(sys.call()[2]), ".subset-", validation_holdout, '.data', sep = '')

  for (i in 1:foldx) {
    prefix = paste("hold",i, sep="")
    fname = paste(prefix, namebase, sep="")

    cv.load.saved.file = F
    if(load.saved.file && file.exists(fname)) {
      load(fname)
      # we only use the runcv function's load saved file if
      #     we are loading the (hopefully corresponding) save
      #     file of our own
      cv.load.saved.file = T

      cat("Loading saved file", fname, "\n")
    }
    else {
      subset = sample(colnames(data), ncol(data)*(1-validation_holdout))
      holdout = data[,which(!colnames(data) %in% subset)]

      formatted_data = rpf_coxph_formatdata(data[,subset], survivaldata, timevar, statusvar)
      if (!is.na(validation_dataset)) {
        #NOTE: currently this requires the validation dataset to be covered by the survival data
        #      that is, for the survival data to be applicable to both training and validation data
        formatted_holdout = rpf_coxph_formatdata(validation_dataset, survivaldata, timevar, statusvar)
      }
      else {
        formatted_holdout = rpf_coxph_formatdata(holdout, survivaldata, timevar, statusvar)
      }
      save(subset, holdout, formatted_data, formatted_holdout, file=fname)
    }

    fitcv = run.cv(formatted_data, training.set.ratio=consensus_holdout,
                   num.repeats=consensus_foldx, trace=3, max.steps=consensus_maxsteps,
                   prefix=paste("hold",i,sep=""), load.saved.file=cv.load.saved.file, ...)

    # now construct a consensus signature
    #browser()
    consensus = consensus_coefficients(fitcv, consensus_level, filter_count=consensus_count_filter,
                                       return.all=F)

    # now apply it to the oritinal dataset and the holdout to
    #     assess performance
    training_pvalue = plot_model_kaplan_meier(fitcv[[1]][[1]],
                                              formatted_data[[1]], consensus_level,
                                              signature=consensus, plotthis=F)

    holdout_pvalue = plot_model_kaplan_meier(fitcv[[1]][[1]],
                                             formatted_holdout[[1]], consensus_level,
                                             signature=consensus, plotthis=F)

    result = list(training_data=formatted_data, training_fitcv=fitcv, consensus_signature=consensus,
                  holdout_data=formatted_holdout, training_pvalue=training_pvalue,
                  holdout_pvalue=holdout_pvalue)

    results = c(results, list(result))
  }
  return(results)
}
