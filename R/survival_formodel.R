#'survival_formodel
#'
#'
#'


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
