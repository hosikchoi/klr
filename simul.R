setwd("/home/ayoung/svm")
dyn.load("Rsvm.so")
library(CVST); library(kernrank); library(MASS); library(caret); library(pROC);
source('SVM (2).R'); source('KLR.R'); source('Kernel.R'); source('simul_func.R');
source('SimulData.R'); source('svm.R'); source('test_func.R'); source('simul_rep.R');
i=1
j=1
k=1
iter=200
# n = c(100, 10000);
p = c(10, 50, 250); D = c(1, 1.5, 2, 2.5, 2.75, 3); rho = c(0, 0.5, -0.5); r = c(0.1, 0.3, 0.5); 

simul_result = simul_rep(n = c(100, 10000), p[i], D[6], rho[j], r[k], iter)
save(simul_result, file = paste0('result_', p[i], '_', rho[j], '_', r[k], '_', iter, '_', n[1], '.Rdata' ))

test_MAE = test_func(simul_result, "MAE")
test_MSE = test_func(simul_result, "MSE")
test_LogL = test_func(simul_result, "LogL")
test_CalB = test_func(simul_result, "CalB")
test_AUC = test_func(simul_result, "AUC")
test_deviance = test_func(simul_result, "deviance")

test_result = list(test_MAE=test_MAE, test_MSE=test_MSE, test_LogL=test_LogL, 
                   test_CalB=test_CalB, test_AUC=test_AUC, test_deviance=test_deviance)
save(test_result, file = paste0('test_result_', p[i], '_', rho[j], '_', r[k], '_', iter, '_', n[1], '.Rdata'))
