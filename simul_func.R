# function for simulation----

simul_func = function (data, kernel, bandwidth, C, lambda, tol, maxiter, rank = T){
  
  # Mallows result ----
  M_piResult_klr = matrix(0, length(data$test$y), 1); M_pred_klr_tmp = matrix(0, length(data$test$y), 1);
  M_piResult_svm = matrix(0, length(data$test$y), 1);  
  M_MISMat_klr = 0; M_MISMat_svm = 0
  
  # Kendall result ----
  K_piResult_klr = matrix(0, length(data$test$y), 1); K_pred_klr_tmp = matrix(0, length(data$test$y), 1);
  K_piResult_svm = matrix(0, length(data$test$y), 1);  
  K_MISMat_klr = 0; K_MISMat_svm = 0
  
  # CVST.params for Mallows ----
  M_klr_par = constructParams(kernel = kernel[1], bandwidth = bandwidth, lambda=lambda, tol=tol, maxiter=maxiter) 
  M_svm_par = constructParams(kernel = kernel[1], bandwidth = bandwidth, C = C)  
  
  # CVST.params for Kendall ----
  K_klr_par = constructParams(kernel = kernel[2], lambda=lambda, tol=tol, maxiter=maxiter) 
  K_svm_par = constructParams(kernel = kernel[2], C = C) 
  
  # Data transformation ----
  if(rank == T) data$train$x = floor(t(apply(data$train$x, 1, rank)))
  if(rank == T) data$test$x  = floor(t(apply(data$test$x, 1, rank)))
  
  # CVST.Train & Test set ----
  Train = constructData(data$train$x, data$train$y)
  Test = constructData(data$test$x, data$test$y)
  rm(list = "data")
  
  # choose (C or lambda) output
  M_svm_C = 0; M_klr_lambda = 0
  K_svm_C = 0; K_klr_lambda = 0
  
  #.. Mallows tuning ----
  
  #..... KLR Prediction ----
  M_KLR = paste("Mallows KLR tuning...:", date(), "\n"); cat(M_KLR)
  M_opt_klr = cv_klr(Train, M_klr_par)
  
  M_klrK = as.kernelMatrix(Mallows_kernel(bandwidth = M_opt_klr[[1]]$bandwidth, Train$x))
  M_klr_TrainK = Matrix(M_klrK)
  M_klr_model = Klogreg(Train$x, kernel[1], M_klr_TrainK, Train$y, 
                        M_opt_klr[[1]]$lambda, M_opt_klr[[1]]$tol, M_opt_klr[[1]]$maxiter)
  M_klr_TestK = Mallows_kernel(bandwidth = M_opt_klr[[1]]$bandwidth, Test$x, Train$x)
  M_pred_klr_tmp = Klogreg_predict(M_klr_model, M_klr_TestK)
  M_pred_klr = M_pred_klr_tmp$result
  M_piResult_klr = M_pred_klr_tmp$pi
  M_MISMat_klr = mean(M_pred_klr != Test$y)
  M_klr_lambda = M_opt_klr[[1]]$lambda
  
  rm(list = c("M_opt_klr", "M_klrK", "M_klr_TrainK", "M_klr_model", "M_klr_TestK", "M_pred_klr_tmp"))
  gc()
  
  #..... SVM Prediction ----
  M_SVM = paste("Mallows SVM tuning...:", date(), "\n"); cat(M_SVM)
  M_opt_svm = cv_svm(Train, M_svm_par)
  
  M_svmK = Mallows_kernel(bandwidth = M_opt_svm[[1]]$bandwidth, Train$x)
  M_svm_TrainK = as.kernelMatrix(M_svmK)
  M_svm_model = svm(M_svm_TrainK, Train$y, kernel = "precomputed",
                    type = "C-classification", cost = M_opt_svm[[1]]$C, scale = F, probability = T)
  M_svm_TestK = as.kernelMatrix(Mallows_kernel(bandwidth = M_opt_svm[[1]]$bandwidth, Test$x, Train$x))
  M_pred_svm = predict(M_svm_model, M_svm_TestK)
  W = t(M_svm_model$coefs) %*% M_svm_model$SV
  fx = t(W %*% t(M_svm_TestK)) - M_svm_model$rho
  M_piResult_svm = 1 / ( 1 + exp(fx) )
  M_MISMat_svm = mean(M_pred_svm != Test$y)
  M_svm_C = M_opt_svm[[1]]$C
  
  rm(list = c("M_opt_svm", "M_svmK", "M_svm_TrainK", "M_svm_model", "M_svm_TestK", "W", "fx"))
  gc()
  
  #.. Kendall tuning ----
  
  #..... KLR Prediction ----
  K_KLR = paste("Kendall KLR tuning...:", date(), "\n"); cat(K_KLR)
  K_opt_klr = cv_klr(Train, K_klr_par) 
  
  K_klrK = as.kernelMatrix(Kendall_kernel(Train$x))
  K_klr_TrainK = Matrix(K_klrK)
  K_klr_model = Klogreg(Train$x, kernel[2], K_klr_TrainK, Train$y, 
                        K_opt_klr[[1]]$lambda, K_opt_klr[[1]]$tol, K_opt_klr[[1]]$maxiter)
  K_klr_TestK = Test
  K_pred_klr_tmp = Klogreg_predict(K_klr_model, K_klr_TestK)
  K_pred_klr = K_pred_klr_tmp$result
  K_piResult_klr = K_pred_klr_tmp$pi
  K_MISMat_klr = mean(K_pred_klr != Test$y)
  K_klr_lambda = K_opt_klr[[1]]$lambda
  
  rm(list = c("K_opt_klr", "K_klrK", "K_klr_TrainK", "K_klr_model", "K_klr_TestK", "K_pred_klr_tmp"))
  gc()
  
  #..... SVM Prediction ----
  K_SVM = paste("Kendall SVM tuning...:", date(), "\n"); cat(K_SVM)
  K_opt_svm = cv_svm(Train, K_svm_par) 
  
  K_svmK = Kendall_kernel(Train$x)
  K_svm_TrainK = as.kernelMatrix(K_svmK)
  K_svm_model = svm(K_svm_TrainK, Train$y, kernel = "precomputed", 
                    type = "C-classification", cost = K_opt_svm[[1]]$C, scale = F, probability = T)
  K_svm_TestK = as.kernelMatrix(Kendall_kernel(Test$x, Train$x))
  K_pred_svm = predict(K_svm_model, K_svm_TestK)
  W = t(K_svm_model$coefs) %*% K_svm_model$SV
  fx = t(W %*% t(K_svm_TestK)) - K_svm_model$rho
  K_piResult_svm = 1 / ( 1 + exp(fx) )
  K_MISMat_svm = mean(K_pred_svm != Test$y)
  K_svm_C = K_opt_svm[[1]]$C
  
  rm(list = c("K_opt_svm", "K_svmK", "K_svm_TrainK", "K_svm_model", "K_svm_TestK", "W", "fx"))
  gc()
  
  # Mallows result ----
  # .. KLR Result ----
  M_Output_klr  = list(M_pred_klr = M_pred_klr, M_MISMat_klr = M_MISMat_klr, M_piResult_klr = M_piResult_klr, M_klr_lambda = M_klr_lambda)
  
  # .. SVM Result ----
  M_Output_svm  = list(M_pred_svm = M_pred_svm, M_MISMat_svm = M_MISMat_svm, M_piResult_svm = M_piResult_svm, M_svm_C = M_svm_C)
  
  # Kendall result ----
  # .. KLR Result ----
  K_Output_klr  = list(K_pred_klr = K_pred_klr, K_MISMat_klr = K_MISMat_klr, K_piResult_klr = K_piResult_klr, K_klr_lambda = K_klr_lambda)
  
  # .. SVM Result ----
  K_Output_svm  = list(K_pred_svm = K_pred_svm, K_MISMat_svm = K_MISMat_svm, K_piResult_svm = K_piResult_svm, K_svm_C = K_svm_C)
  
  Test_y = Test$y
  
  Total_Output = list(K_Output_klr = K_Output_klr, K_Output_svm = K_Output_svm,
                      M_Output_klr = M_Output_klr, M_Output_svm = M_Output_svm, Test_y = Test_y)
  return(Total_Output)
}