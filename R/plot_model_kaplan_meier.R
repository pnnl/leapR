#'plot_model_kaplan_meier
#'
#'
#'
#'

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
