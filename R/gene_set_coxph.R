#'gene_set_coxph
#'
#'
#'

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
