# SVM (with ksvm) ----
# Cross validation for tuning ----

cv_svm = function(data, params, fold=5, verbose=F) {
  stopifnot(class(data) == "CVST.data" && class(params) == "CVST.params")
  nParams = length(params)
  dimnames = list(as.character(1:fold), names(params))
  results = matrix(0, fold, nParams, dimnames=dimnames)
  validationIndex = createFolds(data$y, k = fold)
  
  for (ind in 1:nParams) {
    p = params[[ind]]
    if (p$kernel == "Mallows_kernel") {
      kmat = Mallows_kernel(bandwidth = p$bandwidth, data$x)
      curK = as.kernelMatrix(kmat)
    } else if (p$kernel == "Kendall_kernel") {
      kmat = Kendall_kernel(data$x)
      curK <- as.kernelMatrix(kmat)
    } 
    for (f in 1:fold) {
      curTest = getSubset(data, validationIndex[[f]])
      curTrainK = as.kernelMatrix(curK[-validationIndex[[f]], ])
      if (p$kernel == "Mallows_kernel") {
        curmodel = svm(curTrainK, data$y[-validationIndex[[f]]], kernel = "precomputed",
                       type = "C-classification", cost = p$C, scale = F)
      } else if (p$kernel == "Kendall_kernel") {
        curmodel = svm(curTrainK, data$y[-validationIndex[[f]]], kernel = "precomputed", 
                       type = "C-classification", cost = p$C, scale = F)
      }
      curTestK = as.kernelMatrix(curK[validationIndex[[f]], ])
      curpred  = predict(curmodel, curTestK)
      results[f, ind] = mean(curTest$y != curpred)
    }
    if (verbose) {
      cat(names(params)[ind], "(", mean(results[, ind]), ")\n")
    }
  }
  winner = which.min(apply(results, 2, mean))
  if (length(winner) == 0) {
    return(NULL)
  } else {
    return(params[winner])
  }
}