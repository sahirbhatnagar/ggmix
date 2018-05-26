## @knitr metrics

his_loss <- new_metric("hisloss", "His loss function",
                       metric = function(model, out) {
                         return((model$mu - out$fit)^2)
                       })

her_loss <- new_metric("herloss", "Her loss function",
                       metric = function(model, out) {
                         return(abs(model$mu - out$fit))
                       })


sqrerr <- new_metric("sqrerr", "squared error",
                     metric = function(model, out) {
                       as.numeric(crossprod(out$beta - model$beta))
                     })

# best_sqrerr <- new_metric("best_sqrerr", "best squared error",
#                           metric = function(model, out) {
#                             min(colMeans(as.matrix(out$beta - model$beta)^2))
#                           })

# nnz <- new_metric("nnz", "number of nonzeros",
#                   metric = function(model, out) {
#                     colSums(as.matrix(out$beta) != 0)
#                   })

eta <- new_metric("eta", "eta",
                  metric = function(model, out) out$eta)

sigma2 <- new_metric("sigma2", "sigma2",
                     metric = function(model, out) out$sigma2)

rmse <- new_metric("rmse", "Root mean squared error",
                   metric = function(model, out) {
                     as.numeric(sqrt(crossprod(out$y - out$yhat)))
                   })

r2 <- new_metric("r2", "R squared",
                 metric = function(model, out) {
                   cor(out$y,out$yhat)^2
                 })

muy <- new_metric("muy", "mean of draw",
                  metric = function(model, out) {
                    mean(model$y)
                  })


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
