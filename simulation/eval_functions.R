## @knitr metrics

eta <- new_metric("eta", "eta",
                  metric = function(model, out) out$eta)

sigma2 <- new_metric("sigma2", "sigma2",
                     metric = function(model, out) out$sigma2)

modelerror <- new_metric("me", "model error",
                     metric = function(model, out) out$model_error)

tpr <- new_metric("tpr", "True Positive Rate",
                  metric = function(model, out) {
                    length(intersect(out$nonzero_names, model$causal))/length(model$causal)
                  })

"%ni%" <- Negate("%in%")

fpr <- new_metric("fpr", "False Positive Rate",
                  metric = function(model, out){
                    # these are the terms which the model identified as zero
                    modelIdentifyZero <- setdiff(colnames(model$x), out$nonzero_names)
                    FPR <- sum(out$nonzero_names %ni% model$causal) / (sum(out$nonzero_names %ni% model$causal) + sum(modelIdentifyZero %in% model$not_causal))
                    FPR
                  })

nactive <- new_metric("nactive", "Number of Active Variables",
                      metric = function(model, out) {
                        length(out$nonzero_names)
                      })
