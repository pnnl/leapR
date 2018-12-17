# Survival analysis
library(survival)

coxph_expression <- function(survivalmatrix, statusvar, timevar, startcolumn) {
  # function to build a coxph model from each expression pattern
  # in a datamatrix.
  
  results = matrix(nrow=ncol(survivalmatrix[,startcolumn:ncol(survivalmatrix)]), ncol=2)
  rownames(results) = colnames(survivalmatrix[,startcolumn:ncol(survivalmatrix)])
  
  for (col in colnames(survivalmatrix)[startcolumn:ncol(survivalmatrix)]) {
    thisform = paste("Surv(",timevar,",",statusvar,") ~", col, sep="")
    model = try(coxph(as.formula(thisform), data=survivalmatrix), silent=F)
    if (class(model) != "try-error") {
      results[col,1] = summary(model)$waldtest[3]
      results[col,2] = summary(model)$logtest[3]
    }
  }
  return(results)
}

coxph_subset <- function(survivalmatrix, statusvar="status", timevar="time", columns, rows, verbose=F) {
  # function to build a coxph model from an input data matrix given a range of columns and rows
  # the orientation of the matrix is columns=proteins, rows=examples
  ids = columns
  results = NA
  
  thisform = paste("Surv(",timevar,",",statusvar,") ~", paste(ids, collapse= "+"), sep="")
  if (verbose) cat(thisform, "\n")
  
  model = try(coxph(as.formula(thisform), data=survivalmatrix[rows,]), silent=T)
  
  if (class(model) != "try-error") {
    results = list(waldp=summary(model)$waldtest[3], logp=summary(model)$logtest[3], model=model)
  }
  return(results)
}

coxph_random_subset <- function(survivalmatrix, statusvar, timevar, startcolumn, size=5, ntimes=50) {
  # function to build a coxph model from each expression pattern
  # in a datamatrix.
  
  results = matrix(nrow=ntimes, ncol=2)
  ids = colnames(survivalmatrix)[startcolumn:ncol(survivalmatrix)]
  
  for (i in 1:ntimes) {
    randsub = sample(ids, size)
    thisform = paste("Surv(",timevar,",",statusvar,") ~", paste(randsub, collapse= "+"), sep="")
    model = try(coxph(as.formula(thisform), data=survivalmatrix))
    if (class(model) != "try-error") {
      results[i,1] = summary(model)$waldtest[3]
      results[i,2] = summary(model)$logtest[3]
    }
  }
  return(results)
}

average_redundant_rows <- function(data) {
  thesenames = unique(rownames(data))
  uniquedata = matrix(0, ncol=ncol(data), nrow=length(thesenames))
  for(i in 1:length(thesenames)) {
    lines = which(rownames(data) == thesenames[i])
    uniquedata[i,] = sapply(colnames(data), function (i) mean(data[lines,i], na.rm=T))
  }
  rownames(uniquedata) = thesenames
  colnames(uniquedata) = colnames(data)
  return(uniquedata)
}

rpf_coxph_formatdata <- function(data, survdata, timevar, statusvar, extra.features=NULL) {
  # format for suvival data analysis using the regularisation path following algorithm
  namemat = as.data.frame(cbind(1:nrow(data), rownames(data)))
  colnames(namemat) = c("ID", "Gene.Symbol")
  
  matchcols = colnames(data)[which(colnames(data) %in% rownames(survdata))]
  
  efeatures = NULL
  if (!is.null(extra.features)) {
    efeatures = survdata[matchcols,extra.features]
  }
  
  stimes = as.numeric(survdata[matchcols,timevar])
 
  status = as.numeric(survdata[matchcols,statusvar])
  
  # some of the status or time values are missing and thus the data can't be used- screen these out
  x = which(!is.na(status)&!is.na(stimes))
  
  matchcols = matchcols[x]
  
  stimes = stimes[x]
  status = status[x]
  efeatures = efeatures[x,]
  
  mat = I(t(data[,matchcols]))
  
  if (is.null(extra.features)) {
    data.struct = list(data.frame(x=mat, 
                                      time=stimes, 
                                      status=status),
                                  namemat)
  } 
  else {
    data.struct = list(list(x=mat, 
                              time=stimes, 
                              status=status,
                              extra.features=efeatures),
                       namemat)
  }
  
  return(data.struct) 
}

consensus_coefficients = function(fitcv, level, return.all=T, filter_count=NULL, meanweight=NULL, subsample=NULL) {
  # get a list of the lengths of each model- then we can screen for which have the
  # appropriate number of levels (some are shorter than max)
  if (!is.null(subsample)) {
    usethese = sample(1:length(fitcv), subsample)
  }
  else {
    usethese = 1:length(fitcv)
  }
  
  fitlengths = sapply(usethese, function (i) length(fitcv[[i]][[1]]$b.corrector))
  
  passed = usethese[which(fitlengths>=level)]
  if (sum(fitlengths>=level)==0) return(NA)
  
  #varsrle = rle(sort(unlist(sapply(1:length(fitcv), function (i) names(fitcv[[i]][[1]]$b.corrector[[level]])))))
  varsrle = try(rle(sort(unlist(sapply(usethese[which(fitlengths>=level)], 
                                       function (i) names(fitcv[[i]][[1]]$b.corrector[[level]]))))))
  
  if (class(varsrle) == "try-error") browser()
  
  #browser()
  results = matrix(nrow=length(fitlengths), ncol=length(varsrle$values))
  colnames(results) = varsrle$values
  
  #browser()
  for (i in 1:length(usethese)) {
      x = usethese[i]
      
      if (level > length(fitcv[[x]][[1]]$b.corrector)) next()
          
      for (j in 1:length(fitcv[[x]][[1]]$b.corrector[[level]])) {
          val = fitcv[[x]][[1]]$b.corrector[[level]][j]
          results[i,names(val)] = val
      } 
  }
  meansig = colMeans(results, na.rm=T)
  names(meansig) = colnames(results)
  
  if (!is.null(meanweight)) {
    meansig = meansig*varsrle$lengths
  }
  
  if (!is.null(filter_count)) {
    meansig = meansig[which(varsrle$lengths>filter_count)]
  }

  #browser()
  
  if (return.all) {
    results_annotated = results
    colnames(results_annotated) = line.to.gene.symbol(colnames(results), fitcv[[level]][[1]]$gene.symbols)
    meansig_annotated = meansig
    names(meansig_annotated) = line.to.gene.symbol(names(meansig_annotated), fitcv[[level]][[1]]$gene.symbols)
    return(list(mat=results, meansig=meansig, counts=varsrle$lengths, 
                meansig_annotated=meansig_annotated, mat_annotated=results_annotated))
  }
  
  return(meansig)
}

consensus_model = function(model_list, return.all=T, filter_count=NULL) {
  # similar to the consensus_coefficients above but works with a
  #         list of models and outputs helpful things :)
  
  varsrle = rle(sort(unlist(sapply(1:length(model_list), function (i) names(model_list[[i]])))))
  
  results = rep(list(NA), length(varsrle$values))
  names(results) = varsrle$values
  
  for (i in 1:length(model_list)) {
      for (j in names(model_list[[i]])) {
          val = model_list[[i]][j]
          results[[j]] = c(results[[j]], val)
      }
    } 
  
  meancoefficients = sapply(names(results), function (n) mean(results[[n]], na.rm=T))
  stddev = sapply(names(results), function (n) sd(results[[n]], na.rm=T))
  counts = sapply(names(results), function (n) sum(!is.na(results[[n]])))
  return(list(mean=meancoefficients, sd=stddev, counts=counts))
}


cv_selected_variable_analysis = function(fitcv, level) {
  varsrle = rle(sort(unlist(sapply(1:length(fitcv), function (i) names(fitcv[[i]][[1]]$b.corrector[[level]])))))
  
  results = matrix(0, nrow=length(fitcv[[1]][[2]]), ncol=length(varsrle$values))
  colnames(results) = varsrle$values
  rownames(results) = names(fitcv[[1]][[2]])
  for (i in 1:length(fitcv)) {
    for (j in names(fitcv[[i]][[1]]$b.corrector[[level]])) {
      results[,j] = results[,j] + as.integer(fitcv[[i]][[2]])
    } 
  }
  return(results)
}

build_predictor = function(fit, data, level=NULL, signature=NULL, genesymbols_from=TRUE, genesymbols_to=NULL, browser=F) {
  if (is.null(signature)) {
    coeffs = fit$b.corrector[[level]]
    #coeffs = fit$b.predictor[[level]]
    ids = names(coeffs)
  }
  else {
    coeffs = signature
    ids = names(signature)
  }
  if (is.null(coeffs)) return(NA)
  
  inters = ids[grep("\\*", ids)]
  inters_coeffs = coeffs[inters]
  singles = ids[grep("\\*", ids, invert=TRUE)]
  singles_coeffs = coeffs[singles]
  
  # first process single variables
  idx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][2])
  udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])
  
  if (!is.null(genesymbols_from)) {
    # these are _levels_ for some odd reason
    #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
    #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
    idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
    
    names(singles_coeffs) = idx
    names(udx) = idx
    idx = idx[which(idx %in% colnames(data[,1]))]
    udx = udx[which(idx %in% colnames(data[,1]))]
  }
  #browser()
  udx[which(udx == "dn")] = -1
  udx[which(udx == "up")] = 1
  
  sigd = data$x[,idx]
  sig = sapply(rownames(sigd), function (r) sum(sigd[r,]*singles_coeffs[idx]*as.integer(udx[idx]), na.rm=T))

  # now process interactions
  for (intid in inters) {
    inters_coeff = inters_coeffs[intid]
    ids = strsplit(intid, "\\*")[[1]]
    
    idx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][2])
    idxp = idx
    udx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][1])
    
    if (!is.null(genesymbols_from)) {
      #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
      #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
      idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
      names(udx) = idx
      idx = idx[which(idx %in% colnames(data[,1]))]
      udx = udx[which(idx %in% colnames(data[,1]))]
    }
    udx[which(udx == "dn")] = -1
    udx[which(udx == "up")] = 1
    
    if (length(idx)<2) next
    sigd = data$x[,idx]
    pdt = rep(1, nrow(sigd))
    
    # calculate column products
    for (i in 1:ncol(sigd)) {
      x = sigd[,i]
      x[which(is.na(x))] = 0
      pdt = pdt * x * as.integer(udx[i])
    }
    # add to model
    sig = sig + (pdt*inters_coeff)
  }
  if (browser) browser()
  return(sig)
}

build_predictor_save = function(fit, data, level=NULL, signature=NULL, genesymbols_from=TRUE, genesymbols_to=NULL, browser=F) {
  if (is.null(signature)) {
    coeffs = fit$b.corrector[[level]]
    #coeffs = fit$b.predictor[[level]]
    ids = names(coeffs)
  }
  else {
    coeffs = signature
    ids = names(signature)
  }
  inters = ids[grep("\\*", ids)]
  singles = ids[grep("\\*", ids, invert=TRUE)]
  
  # first process single variables
  idx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][2])
  udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])
  
  if (!is.null(genesymbols_from)) {
    # these are _levels_ for some odd reason
    #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
    #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
    idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
    
  }
  browser()
  udx[which(udx == "dn")] = 1
  udx[which(udx == "up")] = -1
  
  sigd = data$x[,as.integer(idx)]
  sig = sapply(rownames(sigd), function (r) sum(sigd[r,]*coeffs[singles]*as.integer(udx), na.rm=T))
  
  # now process interactions
  for (intid in inters) {
    coeff = coeffs[intid]
    ids = strsplit(intid, "\\*")[[1]]
    
    idx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][2])
    udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])
    
    if (!is.null(genesymbols_from)) {
      #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
      #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
      idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
      idx = which(colnames(data[,1]) %in% idx)
    }
    udx[which(udx == "dn")] = 1
    udx[which(udx == "up")] = -1
    
    sigd = data$x[,as.integer(idx)]
    pdt = rep(1, nrow(sigd))
    
    # calculate column products
    for (i in 1:ncol(sigd)) {
      x = sigd[,i]
      x[which(is.na(x))] = 0
      pdt = pdt * x * as.integer(udx[i])
    }
    
    # add to model
    sig = sig + (pdt*coeff)
  }
  
  if (browser) browser()
  return(sig)
}

build_predictor_simple = function(fit, data, level=NULL, signature=NULL, genesymbols_from=TRUE, genesymbols_to=NULL, browser=F) {
  if (is.null(signature)) {
    coeffs = fit$b.corrector[[level]]
    #coeffs = fit$b.predictor[[level]]
    ids = names(coeffs)
  }
  else {
    coeffs = signature
    ids = names(signature)
  }
  if (is.null(coeffs)) return(NA)
  
  inters = ids[grep("\\*", ids)]
  inters_coeffs = coeffs[inters]
  singles = ids[grep("\\*", ids, invert=TRUE)]
  singles_coeffs = coeffs[singles]
  
  # first process single variables
  idx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][2])
  udx = sapply(singles, function (id) strsplit(id, "\\.")[[1]][1])
  
  if (!is.null(genesymbols_from)) {
    # these are _levels_ for some odd reason
    #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
    #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
    idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
    
    names(singles_coeffs) = idx
    names(udx) = idx
    idx = idx[which(idx %in% rownames(data))]
    udx = udx[which(idx %in% rownames(data))]
  }
  #browser()
  udx[which(udx == "dn")] = -1
  udx[which(udx == "up")] = 1
  
  sigd = t(data[idx,])
  sig = sapply(rownames(sigd), function (r) sum(sigd[r,]*singles_coeffs[idx]*as.integer(udx[idx]), na.rm=T))
  
  # now process interactions
  for (intid in inters) {
    inters_coeff = inters_coeffs[intid]
    ids = strsplit(intid, "\\*")[[1]]
    
    idx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][2])
    idxp = idx
    udx = sapply(ids, function (id) strsplit(id, "\\.")[[1]][1])
    
    if (!is.null(genesymbols_from)) {
      #idx = as.character(genesymbols_from[as.numeric(idx), "Gene.Symbol"])
      #idx = which(genesymbols_to[, "Gene.Symbol"] %in% idx)
      idx = as.character(fit$gene.symbols[as.numeric(idx), "Gene.Symbol"])
      names(udx) = idx
      idx = idx[which(idx %in% rownames(data))]
      udx = udx[which(idx %in% rownames(data))]
    }
    udx[which(udx == "dn")] = -1
    udx[which(udx == "up")] = 1
    
    if (length(idx)<2) next
    sigd = t(data[idx,])
    pdt = rep(1, nrow(sigd))
    
    # calculate column products
    for (i in 1:ncol(sigd)) {
      x = sigd[,i]
      x[which(is.na(x))] = 0
      pdt = pdt * x * as.integer(udx[i])
    }
    # add to model
    sig = sig + (pdt*inters_coeff)
  }
  if (browser) browser()
  return(sig)
}



plot_model_kaplan_meier = function(fit, data, level, signature=NULL, prediction=NULL,
                                   training.set=NULL, group.ratios=c(0.5,0.5), plotthis=T, extra=NULL,
                                   return.survival=F, genesymbols_from=T, ...) {
  if (is.null(prediction)) prediction = try(build_predictor(fit, data, level, signature, genesymbols_from))
  
  if (class(prediction) == "try-error") {
    if (return.survival) return(list(pvalue=NA, surv=NA, logrank=NA, 
                                     prediction=prediction, raw_pred=rep(NA, ncol(data$x))))
      return(NA)
    
  }
  
  raw_pred = prediction
  
  prediction = (prediction > prediction[order(prediction)]
                [floor(length(data$time)*group.ratios[1])]) +
               (prediction > prediction[order(prediction)]
                 [ceiling(length(data$time)*group.ratios[2])])
  
  if (is.null(training.set)) {
    prediction = (prediction>=1)+(prediction>=2)
    
    if (is.null(extra)) {
      #browser()
      surv <- try(survfit(Surv(data$time[prediction != 1], 
                        data$status[prediction != 1] == 1)~prediction[prediction != 1]), silent=T);
      #if (!is.numeric(surv)) surv = NA
    }
    else {
      surv <- survfit(Surv(data$time[prediction != 1], 
                           data$status[prediction != 1] == 1)~prediction[prediction != 1]+extra[prediction !=1]);
    }
    logrank = NA
    p.val = NA
    try({
      if (is.null(extra)) {
        logrank <- survdiff(Surv(data$time[prediction != 1], 
                           data$status[prediction != 1] == 1)~prediction[prediction != 1]);
      }
      else { 
        logrank <- survdiff(Surv(data$time[prediction != 1], 
                                 data$status[prediction != 1] == 1)~prediction[prediction != 1]+extra[prediction !=1]);
      }
        p.val = pchisq(logrank$chisq, 1, lower.tail=FALSE)
    })
    if (plotthis) plot(surv, conf.int=F, col=c("blue", "red"), main=paste("Signature p-value:", signif(p.val, 5)), ...)
    if (!return.survival) return(p.val)
    return(list(pvalue=p.val, surv=surv, logrank=logrank, prediction=prediction, raw_pred=raw_pred))
  }
  training = prediction[which(training.set)]
  testing = prediction[which(!training.set)]
  training.time = data$time[which(training.set)]
  testing.time = data$time[which(!training.set)]
  training.status = data$status[which(training.set)]
  testing.status = data$status[which(!training.set)]
  
  training = (training>=1)+(training>=2)
  testing = (testing>=1)+(testing>=2)
  
  #return(list(training.time, training.status, training))
#   training.surv = survfit(Surv(training.time,
#                                training.status ==1)~training[training != 1])
#   #return(list(testing, testing.time, testing.status))
#   testing.surv = survfit(Surv(testing.time,
#                               testing.status == 1)~testing[testing != 1])
#   
#   training.logrank = survdiff(Surv(training.time,
#                                    training.status ==1)~training[training != 1])
#   testing.logrank = survdiff(Surv(testing.time,
#                                   testing.status == 1)~testing[testing != 1])
#   
  training.surv = survfit(Surv(training.time,
                               training.status ==1)~training)
  #return(list(testing, testing.time, testing.status))
  testing.surv = survfit(Surv(testing.time,
                              testing.status == 1)~testing)
  
  training.logrank = survdiff(Surv(training.time,
                                   training.status ==1)~training)
  testing.logrank = survdiff(Surv(testing.time,
                                  testing.status == 1)~testing)
  
  training.p.val = pchisq(training.logrank$chisq, 1, lower.tail=FALSE)
  testing.p.val = pchisq(testing.logrank$chisq, 1, lower.tail=FALSE)
  if (plotthis) plot(training.surv, conf.int=F, col=c("blue", "red"), main=paste("Training p-value:", signif(training.p.val, 5)), ...)
  if (plotthis) plot(testing.surv, conf.int=F, col=c("blue", "red"), main=paste("Testing p-value:", signif(testing.p.val, 5)), ...)
  return(list(training=training.p.val, testing=testing.p.val, prediction=raw_pred))
}

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

cross_validation_matrix = function(validation_dataset, survivaldata, timevar, statusvar, 
                                    fitcv, data, consensus_foldx=25, consensus_maxsteps=100, 
                                   meanweight=NULL, ...) {
  # controller function to take a (largish) dataset and build a consensus
  # signature using a cross-validation approach _on_a_subset_of_the_total_data
  # then apply the consensus to the held out set as a more appropriate validation
  
  results_training = matrix(ncol=consensus_maxsteps, nrow=consensus_foldx)
  results_testing = matrix(ncol=consensus_maxsteps, nrow=consensus_foldx)
  
  formatted_data = rpf_coxph_formatdata(data, survivaldata, timevar, statusvar)
  formatted_holdout = rpf_coxph_formatdata(validation_dataset, survivaldata, timevar, statusvar) 
    
  for (level in 1:consensus_maxsteps) {
    # now construct a consensus signature
    #browser()
    for (filt in 1:consensus_foldx) {
      
      consensus = consensus_coefficients(fitcv, level, filter_count=filt,
                                         return.all=F, subsample=consensus_foldx,
                                         meanweight=meanweight)
      
      if (length(consensus)==0) next()
      
      #browser()
      # now apply it to the oritinal dataset and the holdout to
      #     assess performance
      training_pvalue = plot_model_kaplan_meier(fitcv[[1]][[1]], 
                                                formatted_data[[1]], level,
                                                signature=consensus, plotthis=F)
      
      holdout_pvalue = plot_model_kaplan_meier(fitcv[[1]][[1]], 
                                               formatted_holdout[[1]], level,
                                               signature=consensus, plotthis=F)
      results_training[filt,level] = training_pvalue
      results_testing[filt,level] = holdout_pvalue
      
    }
  }
  return(list(training=results_training, testing=results_testing))
}


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

gene_set_coxph <- function(data, survivaldata, timevar, statusvar, geneset,
                          min.geneset.size=6, num.repeats=10, randomize.order=F, 
                          randomize.genes=F, match.set.size=F, 
                          extra.features=NULL, sigthresh=0.05,
                          translation_table=NULL, training.set.ratio=0.75,
                          sum.sig=F, sig.vote=F) {
  # for each gene set construct a separate coxph model
  
  results = list()
  
  overlap = rownames(survivaldata)[which(rownames(survivaldata) %in% colnames(data))]
  formatdata = cbind(survivaldata[overlap,statusvar], survivaldata[overlap,timevar], t(data[,overlap]))
  colnames(formatdata)[1:2] = c(statusvar, timevar)
  
  patientids = colnames(data)[which(colnames(data) %in% rownames(survivaldata))]
  
  # STUPID KLUDGE
  formatdata[which(is.na(formatdata))] = 0
  colnames(formatdata) = c(statusvar, timevar, rownames(data))
  
  runorder = 1:length(geneset$names)
  if (randomize.order) runorder = sample(1:length(geneset$names))
   
  for (i in runorder) {
    thisname = geneset$names[i]
    thissize = geneset$sizes[i]
    thisdesc = geneset$desc[i]
    grouplist = geneset$matrix[i,1:thissize]
    
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
      result = NULL
    }
    else {
      result = NULL
      cat("Fitting", thisname, "with", x, "variables\n")
      
      groupids = grouplist[which(grouplist %in% rownames(data))]
      
      if (num.repeats > 0) {
        cat("...running CV for", thisname, "\n")
        result = list()
        signatures = matrix(ncol=ncol(data), nrow=num.repeats)
        colnames(signatures) = colnames(data)
        
        models = matrix(nrow=num.repeats, ncol=length(groupids))
        colnames(models) = groupids
               
        for (c in 1:num.repeats) {
          training_samples = sample(patientids, 
                                    training.set.ratio*ncol(data))
          testing_samples = patientids[which(!patientids %in% training_samples)]
          
          plist = coxph_expression(as.data.frame(formatdata[training_samples,c(statusvar, timevar, groupids)]), statusvar, timevar, 3)
          
          siggroupids = sample(rownames(plist)[which(plist[,2]<sigthresh)])
          
          cat("...selecting", length(siggroupids), "significant components for use in model\n")
          
          res = coxph_subset(as.data.frame(formatdata), statusvar, timevar, siggroupids, training_samples, verbose=F)
          
          if (is.na(res)) next
          
          if (!sum.sig)
            # now apply the model to the testing set
            signature = sapply(testing_samples, function (n) 
                    sum(res$model$coefficients * formatdata[n,names(res$model$coefficients)]))
         
          else 
            # this just tests the thing to see if abundance values alone can make
            # a predictive model simply summed
            signature = sapply(testing_samples, function (n) 
                mean(1 * formatdata[n,names(res$model$coefficients)]))
          
          signatures[c,testing_samples] = signature
          models[c,names(res$model$coefficients)] = res$model$coefficients
          
          #return(signature)
          cvp = try(plot_model_kaplan_meier(1, list(time=formatdata[testing_samples,timevar], 
                                                status=formatdata[testing_samples,statusvar]), 1, 
                                        prediction=signature, plotthis=F))
          
          if (class(cvp) == "try-error") cvp=NA
          
          this = list(training=res, testing_p=cvp, testing_signature=signature)
          result = c(result, list(this))
          
          cat(cvp, "\n")
          
        }
      }
      
      if (sig.vote) 
      {
        signatures_rank = t(apply(signatures, 1, rank))
        signatures_rank[which(is.na(signatures))] = NA
        signatures = signatures_rank
      }
      signatures_int = colMeans(signatures, na.rm=T)
      
      intsigp = plot_model_kaplan_meier(1, list(time=formatdata[,timevar], 
                                                status=formatdata[,statusvar]), 1, 
                                        prediction=signatures_int, plotthis=F)
      cat("Integrated p-value:", intsigp, "\n")
      
      this = list(results=result, integrated.p=intsigp, signatures=signatures, models=models)
      
      results = c(results, list(this))
    }
  }
  names(results) = geneset$names
  
  return(results)
}

pathway_model_predictor_matrix = function(pathfit, validation_data=NULL) {
  # this builds a matrix of predictions of pathway models from the data
  np = length(pathfit$fits)
  result = matrix(ncol=nrow(pathfit$data[[1]][[1]]), nrow=np)
  colnames(result) = rownames(pathfit$data[[1]][[1]])
  rownames(result) = names(pathfit$fits)
  
  pvalues = rep(1, length(names(pathfit$fits)))
  names(pvalues) = names(pathfit$fits)
  
  for (i in 1:np) {
    thisfit = pathfit$fits[[i]]
    
    pathname = names(pathfit$bestmodels[i])
    bestx = which.min(pathfit$pvalues[i,])
    
    if (any(sapply(1:ncol(pathfit$pvalues), function (c) !is.na(min(pathfit$pvalues[i,c]))))) {
      #cat(pathname, " - ", bestx, "\n")
      if (!is.null(validation_data)) {
            sig = build_predictor(thisfit, validation_data[[1]], bestx, 
                                  genesymbols_from=thisfit$gene.symbols, genesymbols_to=validation_data[[2]])
            model_pvalue = plot_model_kaplan_meier(thisfit, 
                                         validation_data[[1]], 0,
                                         prediction=sig, plotthis=F)
            pvalues[pathname] = model_pvalue
      }
      sig = build_predictor(thisfit, pathfit$data[[i]][[1]], bestx) 
      result[pathname,] = sig
    }
  }
  if (!is.null(validation_data)) {
    return(list(matrix=result, pvalues=pvalues))
  }
  return(result)
}

print.model <- function(fit, level=NULL, pvalues=NULL) {
  if (is.null(level)) {
    level = which.min(pvalues)
  }
  model = fit$b.corrector[[level]]
  names(model) = line.to.gene.symbol(names(model), fit$gene.symbols)
  
  return(model)
}

cv_significance_summary <- function(mintable, results_struct, pthresh=0.05, return.percentage=T) {
  # mintable is a matrix of pathways (rows) by cv runs (cols) with the min pvalue
  #         from the cv in each cell
  # results is the results structure from the LEAP run which includes which tumors
  #         were included in the training and testing sets
  # The output is a table of pathways by tumors that provides the counts of the times
  #         where the tumor was included in a cv testing set that gave the pathway
  #         as significant by the pthresh
  results = matrix(nrow=nrow(mintable), ncol=length(results_struct$fits[[1]][[2]]$training.set))  
  rownames(results) = rownames(mintable)
  colnames(results) = rownames(results_struct$data[[1]][[1]])
  
  counts = results
  
  rmatrix = sapply(1:ncol(mintable), function (i) (!results_struct$fits[[i]][[2]]$training.set))
  
  for (path in rownames(mintable)) {
    sigvec = matrix(mintable[path,]<pthresh, ncol=ncol(mintable), nrow=ncol(results))
    if (return.percentage == T) results[path,] = rowSums(sigvec * rmatrix)/rowSums(rmatrix)
    else results[path,] = rowSums(sigvec * rmatrix)
    counts[path,] = rowSums(rmatrix)
  }
  return(list(results, counts))
}

clinical_data_table_1 = function(patientlist, clinical_data) {
  # we presuppose certain column names
  a = mean(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  b = sd(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  c = min(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  d = max(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  e = sum(clinical_data[patientlist,"PlatinumStatus"]=="Sensitive")
  f = sum(clinical_data[patientlist,"PlatinumStatus"]=="Resistant")
  
  g = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IIA", "IIB", "IIC", "IID", "II"))
  h = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IIIA", "IIIB", "IIIC", "IIID", "III"))
  i = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IVA", "IVB", "IVC", "IVD", "IV"))
  
  j = sum(clinical_data[patientlist,"tumor_grade"] == "G2")
  k = sum(clinical_data[patientlist,"tumor_grade"] == "G3")
  
  l = sum(clinical_data[patientlist,"tumor_residual_disease"] %in% c("1-10 mm", "No Macroscopic disease"))
  m = length(patientlist)-l
  
  n = sum(clinical_data[patientlist,"vital_status"]=="DECEASED")
  o = length(patientlist)-n
  
  p = median(clinical_data[patientlist,"days.to.death.or.last_followup"])
  q = sd(clinical_data[patientlist,"days.to.death.or.last_followup"])
  r = median(as.numeric(clinical_data[patientlist,"days_to_tumor_prog_or_recur"]), na.rm=T)
  s = sd(as.numeric(clinical_data[patientlist,"days_to_tumor_prog_or_recur"]), na.rm=T)
  
  return(c(paste(format(a, digits=1), " (", format(b, digits=3),")", sep=""),paste(c,d,sep="-"),
           "",
           paste(g, " (", format(g*100/length(patientlist), digits=2), "%)", sep=""),
           paste(h, " (", format(h*100/length(patientlist), digits=3), "%)", sep=""),
           paste(i, " (", format(i*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(j, " (", format(j*100/length(patientlist), digits=3), "%)", sep=""),
           paste(k, " (", format(k*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(l, " (", format(l*100/length(patientlist), digits=3), "%)", sep=""),
           paste(m, " (", format(m*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(e, " (", format(e*100/length(patientlist), digits=3), "%)", sep=""),
           paste(f, " (", format(f*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(o, " (", format(o*100/length(patientlist), digits=3), "%)", sep=""),
           paste(n, " (", format(n*100/length(patientlist), digits=3), "%)", sep=""),
           paste(r, " (", format(s, digits=4), ")", sep=""),
           paste(p, " (", format(q, digits=4), ")", sep="")
           ))
}

survival_formodel = function(genename, leapstruct, dataset, consensusmat) {
  winner = which.min(pvalues_complete[genename,])
  surv_struct = plot_model_kaplan_meier(leapResults$fits[[genename]],
                                        dataset[[1]], winner,
                                        return.survival=T)
  consensus_sig = rbind(consensusmat,
                        surv_struct$prediction)
  consensus_pred = sapply(colnames(consensus_sig), 
                          function (c) sum(consensus_sig[,c]))
  return(plot_model_kaplan_meier(NA, dataset[[1]], NA,
                                 prediction=consensus_pred, return.survival=T))
  
}

