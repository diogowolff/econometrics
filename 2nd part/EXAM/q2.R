set.seed(1337)

library(evd)

X_mat_500 = matrix(rlnorm(100*500, 5, 3), nrow = 500)
X_mat_1000 = matrix(rlnorm(100*1000, 5, 3), nrow = 1000)

eps0_mat_500 = matrix(rgumbel(100*500), nrow = 500)
eps1_mat_500 = matrix(rgumbel(100*500), nrow = 500)

eps0_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)
eps1_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)


Y_mat_500 = ifelse(0.8 + 0.7*X_mat_500 + eps1_mat_500 >= eps0_mat_500, 1, 0)
Y_mat_1000 = ifelse(0.8 + 0.7*X_mat_1000 + eps1_mat_1000 >= eps0_mat_1000, 1, 0)

# c)




data_y_mean_500 = colMeans(Y_mat_500)
data_cov_500 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_500), 
                      as.data.frame(X_mat_500))

alpha_gmm_minimizer_500 = function(param, col_index) {
  data_x = matrix(rep(X_mat_500[, col_index], 30), nrow = 500)
  sim_eps0 = matrix(rgumbel(30*500), nrow = 500)
  sim_eps1 = matrix(rgumbel(30*500), nrow = 500)
  
  estimated_choice = ifelse(param[1] + param[2]*data_x + sim_eps0 >= sim_eps1, 1, 0)
  
  
  
  mean_of_sims = colMeans(estimated_choice)
  
  mean_of_means = mean(mean_of_sims)
  
  
  
  
  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(data_x))
  
  mean_of_covs = mean(cov_of_sims_with_x)
  
  
  moment_diff_vec = c(data_y_mean_500[col_index] - mean_of_means, 
                      data_cov_500[col_index] - mean_of_covs)
  
  t(moment_diff_vec) %*% moment_diff_vec
}

results = purrr::map(1:100, ~ optim(c(0.5, 0.5), alpha_gmm_minimizer_500, col_index = .x)$par)
alpha_est_500 = t(matrix(unlist(results), nrow = 2))

colMeans(alpha_est_500)
