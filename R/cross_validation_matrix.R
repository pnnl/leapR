#' cross_validation_matrix
#'
#' cross_validation_matrix function description is...
#'
#' @param validation_dataset is...
#' @param survivaldata is...
#' @param timevar is...
#' @param statusvar is...
#' @param fitcv is...
#' @param data is...
#' @param consensus_foldx default is 25
#' @param consensus_maxsteps default is 100
#' @param meanweight default is NULL
#' @param ... is...
#'
#' @examples
#' dontrun{
#'
#'
#' }
#'
#' @export
#'

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
