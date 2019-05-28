# Data generation function for simulation----

simul_data = function(n, p, rho, D, r){
  mu_vec1 = c(rep(D/(2*sqrt(p*r)), p*r), rep(0, p - p*r))
  mu_vec2 = -c(rep(D/(2*sqrt(p*r)), p*r), rep(0, p - p*r))
  
  sig_mat = diag(1, p, p)
  sig_mat[1:p*r, 1:p*r] = rho^abs(outer(1:p, 1:p, '-'))[1:p*r, 1:p*r]
  
  x1 = mvrnorm(n/2, mu_vec1, sig_mat, T)
  x2 = mvrnorm(n/2, mu_vec2, sig_mat, T)
  df1 = cbind(y = 1, x1); df2 = cbind(y = 0, x2)
  df = rbind(df1, df2)
  Df = df[sample(nrow(df)),]
  simul_dat = list(x = Df[,-1], y = as.factor(Df[,1])) 
  return(simul_dat)
}