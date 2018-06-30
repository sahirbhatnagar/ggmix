## @knitr metrics

eta <- new_metric("eta", "eta",
                  metric = function(model, out) out$eta)

sigma2 <- new_metric("sigma2", "sigma2",
                     metric = function(model, out) out$sigma2)

modelerror <- new_metric("me", "model error",
                     metric = function(model, out) out$model_error)

prederror <- new_metric("prederror", "Prediction error",
                         metric = function(model, out) out$prediction_error)

tpr <- new_metric("tpr", "True Positive Rate",
                  metric = function(model, out) {
                    length(intersect(out$nonzero_names, out$causal))/length(out$causal)
                  })

"%ni%" <- Negate("%in%")

fpr <- new_metric("fpr", "False Positive Rate",
                  metric = function(model, out){
                    FP <- length(setdiff(out$nonzero_names, out$causal)) # false positives
                    TN <- length(out$not_causal) # True negatives
                    FPR <- FP / (FP + TN)
                    FPR
                  })

nactive <- new_metric("nactive", "Number of Active Variables",
                      metric = function(model, out) {
                        length(out$nonzero_names)
                      })


correct_sparsity <- new_metric("correct_sparsity", "Correct Sparsity",
                               metric = function(model, out){
  causal <- out$causal
  not_causal <- out$not_causal
  active <- out$nonzero_names
  p <- out$p

  correct_nonzeros <- sum(active %in% causal)
  correct_zeros <- length(setdiff(not_causal, active))
  #correct sparsity
  (1 / p) * (correct_nonzeros + correct_zeros)
})


mse <- new_metric("mse", "Test Set MSE",
                  metric = function(model, out) {
                    as.numeric(crossprod(out$yhat - out$y) / (length(out$y)))
                  })
