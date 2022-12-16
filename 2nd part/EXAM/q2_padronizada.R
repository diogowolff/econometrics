set.seed(1337)

library(evd)
library(furrr)

######
# Q2 #
######

### (a)

# Check LaTeX pdf document.

### (b)

# Creating matrices containing 100 vectors of incomes, 
# with size 500 and 1000 each, respectively, generated from the 
# log-normal distribution with location 5 and scale sqrt(3).
X_mat_500 = matrix(rlnorm(100*500, 5, sqrt(3)), nrow = 500) 
X_mat_1000 = matrix(rlnorm(100*1000, 5, sqrt(3)), nrow = 1000)

# Creating matrices containing 100 vectors of alternative-specific shocks, 
# with size 500 and 1000 each, respectively, generated from the 
# gumbel distribution.

eps0_mat_500 = matrix(rgumbel(100*500), nrow = 500)
eps1_mat_500 = matrix(rgumbel(100*500), nrow = 500)

eps0_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)
eps1_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)

# Generating 500 choice vectors and 1000 choice vectors, respectively,
# using data generated above. 
Y_mat_500 = ifelse(0.8 + 0.7*X_mat_500 + eps1_mat_500 >= eps0_mat_500, 1, 0)
Y_mat_1000 = ifelse(0.8 + 0.7*X_mat_1000 + eps1_mat_1000 >= eps0_mat_1000, 1, 0)

# Generating the requested datasets:
datasets_500 <- list()
for (i in 1:100) {
datasets_500[[i]] <- data.frame("x" = X_mat_500[,i], "d" = Y_mat_500[,i])
}

# 100 datasets with 500 generated observations of x and d, as requested:
View(datasets_500)

datasets_1000 <- list()
for (i in 1:100) {
datasets_1000[[i]] <- data.frame("x" = X_mat_1000[,i], "d" = Y_mat_1000[,i])
}

# 100 datasets with 1000 generated observations of x and d, as requested:
View(datasets_1000)

### (c)

plan(multisession) # For computational performance.

## SIZE 500 DATASETS ESTIMATES, S = 30:

# Generating the "actual" moments, using the data generated in (b). 
data_y_mean_500 = colMeans(Y_mat_500) # The "actual" means.
data_cov_500 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_500), 
                      as.data.frame(X_mat_500)) # The "actual" covariances.

# Defining the SMM estimator function for datasets of size 500, with S = 30.
alpha_gmm_minimizer_30_500 = function(param, col_index) {
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

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-500 datasets
# generated in (a), for S = 30: 
results_30_500 = future_map(1:100, ~ optim(c(0.5, 0.5), 
                                    alpha_gmm_minimizer_30_500, col_index = .x,
                                    method = 'BFGS',
                                    control = list('maxit' = '1000')),
                                   .options = furrr_options(seed = T))

# Checking convergence status:
convergence_30_500 = purrr::map(1:100, ~ results_30_500[[.x]]$convergence) 

# Results: 
alpha_est_30_500 = t(matrix(unlist(purrr::map(1:100, ~ results_30_500[[.x]]$par)), nrow = 2))

colMeans(alpha_est_30_500) # Mean of alpha_0 and alpha_1 estimates, respectively.
var(alpha_est_30_500[,1]) # Variance of alpha_0 estimates.
var(alpha_est_30_500[,2]) # Variance of alpha_1 estimates.

mean((alpha_est_30_500[,1]-0.8)^2) # MSQE of alpha_0 estimates.
mean((alpha_est_30_500[,2]-0.7)^2) # MSQE of alpha_1 estimates.

## SIZE 1000 DATASETS ESTIMATES, S = 30:

# Generating the "actual" moments, using the data generated in (b). 
data_y_mean_1000 = colMeans(Y_mat_1000)
data_cov_1000 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_1000), 
                       as.data.frame(X_mat_1000))

# Defining the SMM estimator function for datasets of size 1000, with S = 30.
alpha_gmm_minimizer_30_1000 = function(param, col_index) {
  
  data_x = matrix(rep(X_mat_1000[, col_index], 30), nrow = 1000)
  sim_eps0 = matrix(rgumbel(30*1000), nrow = 1000)
  sim_eps1 = matrix(rgumbel(30*1000), nrow = 1000)
  
  estimated_choice = ifelse(param[1] + param[2]*data_x + sim_eps0 >= sim_eps1, 1, 0)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(data_x))
  
  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_1000[col_index] - mean_of_means, 
                      data_cov_1000[col_index] - mean_of_covs)
  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-1000 datasets
# generated in (b), for S = 30: 
results_30_1000 = future_map(1:100, ~ optim(c(0.5, 0.5), 
                                            alpha_gmm_minimizer_30_1000, col_index = .x,
                                            method = 'BFGS',
                                            control = list('maxit' = '1000')),
                                            .options = furrr_options(seed = T))

# Checking convergence status: 
convergence_30_1000 = purrr::map(1:100, ~ results_30_1000[[.x]]$convergence)

# Results:
alpha_est_30_1000 = t(matrix(unlist(purrr::map(1:100, ~ results_30_1000[[.x]]$par)), nrow = 2))

colMeans(alpha_est_30_1000) # Mean of alpha_0 and alpha_1 estimates, respectively.
var(alpha_est_30_1000[,1]) # Variance of alpha_0 estimates.
var(alpha_est_30_1000[,2]) # Variance of alpha_1 estimates. 

mean((alpha_est_30_1000[,1]-0.8)^2) # MSQE of alpha_0 estimates.
mean((alpha_est_30_1000[,2]-0.7)^2) # MSQE of alpha_1 estimates.

### (d)

## SIZE 500 DATASETS ESTIMATES, S = 100:

# Defining the SMM estimator function for datasets of size 500, with S = 100.
alpha_gmm_minimizer_100_500 = function(param, col_index) {
  data_x = matrix(rep(X_mat_500[, col_index], 100), nrow = 500)
  sim_eps0 = matrix(rgumbel(100*500), nrow = 500)
  sim_eps1 = matrix(rgumbel(100*500), nrow = 500)
  
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

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-500 datasets
# generated in (b), for S = 100: 
results_100_500 = future_map(1:100, ~ optim(c(0.5, 0.5), 
                                          alpha_gmm_minimizer_100_500, col_index = .x,
                                          method = 'BFGS',
                                          control = list('maxit' = '1000')),
                                          .options = furrr_options(seed = T))

# Checking convergence status:
convergence_100_500 = purrr::map(1:100, ~ results_100_500[[.x]]$convergence)

# Results: 
alpha_est_100_500 = t(matrix(unlist(purrr::map(1:100, ~ results_100_500[[.x]]$par)), nrow = 2))

colMeans(alpha_est_100_500) # Mean of alpha_0 and alpha_1 estimates, respectively.
var(alpha_est_100_500[,1]) # Variance of alpha_0 estimates.
var(alpha_est_100_500[,2]) # Variance of alpha_1 estimates.

mean((alpha_est_100_500[,1]-0.8)^2) # MSQE of alpha_0 estimates.
mean((alpha_est_100_500[,2]-0.7)^2) # MSQE of alpha_1 estimates.

## SIZE 1000 DATASETS ESTIMATES, S = 100:

# Defining the SMM estimator function for datasets of size 1000, with S = 100.
alpha_gmm_minimizer_100_1000 = function(param, col_index) {
  data_x = matrix(rep(X_mat_1000[, col_index], 100), nrow = 1000)
  sim_eps0 = matrix(rgumbel(100*1000), nrow = 1000)
  sim_eps1 = matrix(rgumbel(100*1000), nrow = 1000)
  
  estimated_choice = ifelse(param[1] + param[2]*data_x + sim_eps0 >= sim_eps1, 1, 0)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(data_x))
  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_1000[col_index] - mean_of_means, 
                      data_cov_1000[col_index] - mean_of_covs)
  
  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-1000 datasets
# generated in (b), for S = 100:
results_100_1000 = future_map(1:100, ~ optim(c(0.5, 0.5), 
                                          alpha_gmm_minimizer_100_1000, col_index = .x,
                                          method = 'BFGS',
                                          control = list('maxit' = '1000')), 
                                          .options = furrr_options(seed = T))

# Checking convergence status:
convergence_100_1000 = purrr::map(1:100, ~ results_100_1000[[.x]]$convergence)

# Results: 
alpha_est_100_1000 = t(matrix(unlist(purrr::map(1:100, ~ results_100_1000[[.x]]$par)), nrow = 2))

colMeans(alpha_est_100_1000) # Mean of alpha_0 and alpha_1 estimates, respectively.
var(alpha_est_100_1000[,1]) # Variance of alpha_0 estimates.
var(alpha_est_100_1000[,2]) # Variance of alpha_1 estimates.

mean((alpha_est_100_1000[,1]-0.8)^2) # MSQE of alpha_0 estimates.
mean((alpha_est_100_1000[,2]-0.7)^2) # MSQE of alpha_1 estimates.


### The End.