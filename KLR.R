# Kerenl Logistic ----
# .. Train & Test ----

Klogreg = function (data, kernel, kernelMat, labels, lambda, tol, maxiter) 
{ 
  y = (labels != levels(labels)[1]) + 0
  N = nrow(kernelMat)
  alpha = rep(1/N, N)
  iter = 1
  while (TRUE) {
    Kalpha = as.vector(kernelMat %*% alpha)
    spec = 1 + exp(-Kalpha)
    pi = 1/spec
    diagW = pi * (1 - pi)
    e = (y - pi)/diagW
    q = Kalpha + e
    theSol = try(solve(kernelMat + lambda * Diagonal(x = 1/diagW), q))
    if (class(theSol) == "try-error") {
      break
    }
    alphan = as.vector(theSol)
    if (any(is.nan(alphan)) || all(abs(alphan - alpha) <= tol)) {
      break
    } else if (iter > maxiter) {
      cat("klogreg:maxiter!")
      break
    } else {
      alpha = alphan
      iter = iter + 1
    }
  }
  yLevels = levels(labels)
  return(list(data = data, kernel = kernel, alpha = as.vector(alpha), pi = pi, yLevels = yLevels))
}

Klogreg_predict = function (klogreg, newData) 
{
  if (klogreg$kernel == 'Mallows_kernel') {
    K = newData %*% klogreg$alpha
  } else if (klogreg$kernel == "Kendall_kernel") {
    kernel = eval(parse(text = klogreg$kernel))
    K = kernelMult(kernel, newData$x, klogreg$data, klogreg$alpha)
  }
  pi = 1/(1 + exp(-as.vector(K)))
  result = (pi >= 0.5) + 0
  result = factor(result, c("0", "1"), klogreg$yLevels, ordered=FALSE)
  output = list(pi=pi, result=result)
  return(output)
}

# .. Cross Validation for tuning ----
cv_klr = function(data, params, fold=5, verbose=F) {
  stopifnot(class(data) == "CVST.data" && class(params) == "CVST.params")
  nParams = length(params)
  dimnames = list(as.character(1:fold), names(params))
  results = matrix(0, fold, nParams, dimnames=dimnames)
  validationIndex = createFolds(data$y, k = fold)
  
  for (ind in 1:nParams) {
    p = params[[ind]]
    if (p$kernel == "Mallows_kernel") {
      kernel = eval(parse(text = p$kernel))
      curK <- as.kernelMatrix(kernel(bandwidth = p$bandwidth, data$x))
    } else if (p$kernel == "Kendall_kernel") {
      kernel = eval(parse(text = p$kernel))
      curK <- as.kernelMatrix(kernel(data$x))
    }
    for (f in 1:fold) {
      curTrain = getSubset(data, -validationIndex[[f]])
      curTest = getSubset(data, validationIndex[[f]])
      curTrainK = Matrix(curK[-validationIndex[[f]],-validationIndex[[f]]])
      curmodel = Klogreg(curTrain$x, p$kernel, curTrainK, curTrain$y, 
                         getN(curTrain) * p$lambda, p$tol, p$maxiter)
      if (p$kernel == "Mallows_kernel") {
        curTestK = Mallows_kernel(bandwidth = p$bandwidth, curTest$x, curTrain$x)
        curpred = Klogreg_predict(curmodel, curTestK)
      } else if (p$kernel == "Kendall_kernel") {
        curTestK = curTest
        curpred = Klogreg_predict(curmodel, curTestK)
      }
      results[f, ind] = mean(curTest$y != curpred$result)
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