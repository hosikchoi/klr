# function for simulation repitation----

simul_rep = function(n, p, D, rho, r, iter){
  Start = paste("Starting the", iter, "times iteration :", date(), "\n")
  cat(Start)
  call = match.call()
  
  result_list = list()
  for(rep in 1:iter) {
    set.seed(rep)
    simul_dat = list(train = simul_data(n[1], p, rho, D, r), test = simul_data(n[2], p, rho, D, r)) 
    
    result = simul_func(data = simul_dat, kernel = c('Mallows_kernel', 'Kendall_kernel'), 
                        bandwidth = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
                        C = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                        lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), 
                        tol = 10e-6, maxiter = 100, rank = T)
    
    ith.message = paste(rep, "result calculation is complete\n")
    cat(ith.message)
    result_list[[rep]] = result
  }
  
  rest = which(sapply(result_list, function(x) is.null(x$K_Output_klr$K_piResult_klr)) == T)
  
  if(length(rest) > 0){
    cat(paste(length(rest), "times iteration restart"))
    
    for(repp in 1:length(rest)){
      set.seed(rest[repp])
      simul_dat = list(train = simul_data(n[1], p, rho, D, r), test = simul_data(n[2], p, rho, D, r)) 
      
      result_list[[rest[repp]]] = simul_func(data = simul_dat, kernel = c('Mallows_kernel', 'Kendall_kernel'), 
                                             bandwidth = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1), 
                                             C = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
                                             lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), 
                                             tol = 10e-6, maxiter = 100, rank = T)
      
      ith.message = paste(repp, "result calculation is complete\n")
      cat(ith.message)
    }
    
  }
  # Mallows result ----
  # .. KLR Result ----
  M_pred_klr = lapply(result_list, function(x) x$M_Output_klr$M_pred_klr)
  M_MISMat_klr = as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$M_Output_klr$M_MISMat_klr)), iter, 1)
  colnames(M_MISMat_klr) = NULL; row.names(M_MISMat_klr) = NULL;
  M_piResult_klr = lapply(result_list, function(x) x$M_Output_klr$M_piResult_klr)
  M_Average_klr = round(colMeans(M_MISMat_klr), 4); M_SE_klr = round(sd(M_MISMat_klr)/sqrt(nrow(M_MISMat_klr)), 4)
  M_MeanSE_klr  = paste0(M_Average_klr, ' (', M_SE_klr, ')')
  M_klr_lambda =  as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$M_Output_klr$M_klr_lambda)), iter, 1)
  colnames(M_klr_lambda) = NULL; row.names(M_klr_lambda) = NULL;
  M_Output_klr  = list(M_pred_klr = M_pred_klr, M_MISMat_klr = M_MISMat_klr, M_piResult_klr = M_piResult_klr, 
                       M_MeanSE_klr = M_MeanSE_klr, M_klr_lambda = M_klr_lambda)
  
  # .. SVM Result ----
  M_pred_svm = lapply(result_list, function(x) x$M_Output_svm$M_pred_svm)
  M_MISMat_svm = as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$M_Output_svm$M_MISMat_svm)), iter, 1)
  colnames(M_MISMat_svm) = NULL; row.names(M_MISMat_svm) = NULL;
  M_piResult_svm = lapply(result_list, function(x) x$M_Output_svm$M_piResult_svm)
  M_Average_svm = round(colMeans(M_MISMat_svm), 4); M_SE_svm = round(sd(M_MISMat_svm[,1])/sqrt(nrow(M_MISMat_svm)), 4)
  M_MeanSE_svm  = paste0(M_Average_svm, ' (', M_SE_svm, ')')
  M_svm_C =  as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$M_Output_svm$M_svm_C)), iter, 1)
  colnames(M_svm_C) = NULL; row.names(M_svm_C) = NULL;
  M_Output_svm  = list(M_pred_svm = M_pred_svm, M_MISMat_svm = M_MISMat_svm, M_piResult_svm = M_piResult_svm, 
                       M_MeanSE_svm = M_MeanSE_svm, M_svm_C = M_svm_C)
  
  # Kendall result ----
  # .. KLR Result ----
  K_pred_klr = lapply(result_list, function(x) x$K_Output_klr$K_pred_klr)
  K_MISMat_klr = as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$K_Output_klr$K_MISMat_klr)), iter, 1)
  colnames(K_MISMat_klr) = NULL; row.names(K_MISMat_klr) = NULL;
  K_piResult_klr = lapply(result_list, function(x) x$K_Output_klr$K_piResult_klr)
  K_Average_klr = round(colMeans(K_MISMat_klr), 4); K_SE_klr = round(sd(K_MISMat_klr)/sqrt(nrow(K_MISMat_klr)), 4)
  K_MeanSE_klr = paste0(K_Average_klr, ' (', K_SE_klr, ')')
  K_klr_lambda =  as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$K_Output_klr$K_klr_lambda)), iter, 1)
  colnames(K_klr_lambda) = NULL; row.names(K_klr_lambda) = NULL;
  K_Output_klr  = list(K_pred_klr = K_pred_klr, K_MISMat_klr = K_MISMat_klr, K_piResult_klr = K_piResult_klr, 
                       K_MeanSE_klr = K_MeanSE_klr, K_klr_lambda = K_klr_lambda)
  
  # .. SVM Result ----
  K_pred_svm = lapply(result_list, function(x) x$K_Output_svm$K_pred_svm)
  K_MISMat_svm = as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$K_Output_svm$K_MISMat_svm)), iter, 1)
  colnames(K_MISMat_svm) = NULL; row.names(K_MISMat_svm) = NULL;
  K_piResult_svm = lapply(result_list, function(x) x$K_Output_svm$K_piResult_svm)
  K_Average_svm = round(colMeans(K_MISMat_svm), 4); K_SE_svm = round(sd(K_MISMat_svm[,1])/sqrt(nrow(K_MISMat_svm)), 4)
  K_MeanSE_svm  = paste0(K_Average_svm, ' (', K_SE_svm, ')')
  K_svm_C =  as.matrix(do.call(rbind.data.frame, lapply(result_list, function(x) x$K_Output_svm$K_svm_C)), iter, 1)
  colnames(K_svm_C) = NULL; row.names(K_svm_C) = NULL;
  K_Output_svm  = list(K_pred_svm = K_pred_svm, K_MISMat_svm = K_MISMat_svm, K_piResult_svm = K_piResult_svm, K_MeanSE_svm = K_MeanSE_svm, K_svm_C = K_svm_C)
  
  Test_y = lapply(result_list, function(x) x$Test_y)
  
  Close = paste("Finishing the", iter, "times iteration :", date(), "\n")
  cat(Close)
  StartClose = c(Start, Close)
  
  Total_Output = list(StartClose = StartClose, call = call, K_Output_klr = K_Output_klr, K_Output_svm = K_Output_svm,
                      M_Output_klr = M_Output_klr, M_Output_svm = M_Output_svm, Test_y = Test_y)
  return(Total_Output)
}
