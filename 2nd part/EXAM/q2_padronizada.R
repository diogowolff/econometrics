set.seed(1337)

library(evd)
library(purrr)

######
# Q2 #
######

### (a)

# Check LaTeX pdf document.

### (b)

# Creating matrices containing 100 vectors of incomes, 
# with size 500 and 1000 each, respectively, generated from the 
# log-normal distribution with location 5 and scale sqrt(3).
X_mat_500 = matrix(rlnorm(100*500, log(25/sqrt(25+3)), log(1 + 3/25)), nrow = 500) 
X_mat_1000 = matrix(rlnorm(100*1000, log(25/sqrt(25+3)), log(1 + 3/25)), nrow = 1000)

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

## SIZE 500 DATASETS ESTIMATES, S = 30:

# Generating the "actual" moments, using the data generated in (b). 
data_y_mean_500 = colMeans(Y_mat_500) # The "actual" means.
data_cov_500 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_500), 
                      as.data.frame(X_mat_500)) # The "actual" covariances.

# Defining the SMM estimator function for datasets of size 500, with S = 30.
alpha_gmm_minimizer_30_500 = function(param, col_index) {
  sim_eps0 = rgumbel(30*500)
  sim_eps1 = rgumbel(30*500)
  
  estimated_choice = matrix(ifelse(param[1] + param[2]*X_mat_500[, col_index] 
                                   + sim_eps0 >= sim_eps1, 1, 0), nrow = 500)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_500[, col_index]))

  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_500[col_index] - mean_of_means, 
                      data_cov_500[col_index] - mean_of_covs)

  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-500 datasets
# generated in (a), for S = 30... 

# 1) With initial guess (0, 0), using BFGS: 

results_30_500_00 = map(1:100, ~ optim(c(0, 0), 
                                    alpha_gmm_minimizer_30_500, col_index = .x,
                                    method = 'BFGS',
                                    control = list('maxit' = '1000')),
                                   )

# Checking convergence status:
convergence_30_500_00 = purrr::map(1:100, ~ results_30_500_00[[.x]]$convergence) 

# Results: 
alpha_est_30_500_00 = t(matrix(unlist(purrr::map(1:100, ~ results_30_500_00[[.x]]$par)), nrow = 2))

means_30_500_00 <- colMeans(alpha_est_30_500_00) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_500_00 <- var(alpha_est_30_500_00[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_500_00 <- var(alpha_est_30_500_00[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_30_500_00 <- mean((alpha_est_30_500_00[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_500_00 <- mean((alpha_est_30_500_00[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_500_00
var_alpha0_30_500_00
var_alpha1_30_500_00
msqe_alpha0_30_500_00
msqe_alpha1_30_500_00

# 2) With initial guess (0.5, 0.5), using BFGS: 

results_30_500_55 = map(1:100, ~ optim(c(0.5, 0.5), 
                                              alpha_gmm_minimizer_30_500, col_index = .x,
                                              method = 'BFGS',
                                              control = list('maxit' = '1000')),
                               )

# Checking convergence status:
convergence_30_500_55 = purrr::map(1:100, ~ results_30_500_55[[.x]]$convergence) 

# Results: 
alpha_est_30_500_55 = t(matrix(unlist(purrr::map(1:100, ~ results_30_500_55[[.x]]$par)), nrow = 2))

means_30_500_55 <- colMeans(alpha_est_30_500_55) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_500_55 <- var(alpha_est_30_500_55[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_500_55 <- var(alpha_est_30_500_55[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_30_500_55 <- mean((alpha_est_30_500_55[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_500_55 <- mean((alpha_est_30_500_55[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_500_55
var_alpha0_30_500_55
var_alpha1_30_500_55
msqe_alpha0_30_500_55
msqe_alpha1_30_500_55

# 3) With initial guess (0.75, 0.75), using BFGS: 

results_30_500_7575 = map(1:100, ~ optim(c(0.75, 0.75), 
                                              alpha_gmm_minimizer_30_500, col_index = .x,
                                              method = 'BFGS',
                                              control = list('maxit' = '1000')),
                               )

# Checking convergence status:
convergence_30_500_7575 = purrr::map(1:100, ~ results_30_500_7575[[.x]]$convergence) 

# Results: 
alpha_est_30_500_7575 = t(matrix(unlist(purrr::map(1:100, ~ results_30_500_7575[[.x]]$par)), nrow = 2))

means_30_500_7575 <- colMeans(alpha_est_30_500_7575) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_500_7575 <- var(alpha_est_30_500_7575[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_500_7575 <- var(alpha_est_30_500_7575[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_30_500_7575 <- mean((alpha_est_30_500_7575[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_500_7575 <- mean((alpha_est_30_500_7575[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_500_7575
var_alpha0_30_500_7575
var_alpha1_30_500_7575
msqe_alpha0_30_500_7575
msqe_alpha1_30_500_7575

## SIZE 1000 DATASETS ESTIMATES, S = 30:

# Generating the "actual" moments, using the data generated in (b). 
data_y_mean_1000 = colMeans(Y_mat_1000)
data_cov_1000 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_1000), 
                       as.data.frame(X_mat_1000))

# Defining the SMM estimator function for datasets of size 1000, with S = 30.
alpha_gmm_minimizer_30_1000 = function(param, col_index) {
  sim_eps0 = rgumbel(30*1000)
  sim_eps1 = rgumbel(30*1000)
  
  estimated_choice = matrix(ifelse(param[1] + param[2]*X_mat_1000[, col_index] + sim_eps0 >= sim_eps1, 1, 0), nrow = 1000)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_1000[, col_index]))
  
  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_1000[col_index] - mean_of_means, 
                      data_cov_1000[col_index] - mean_of_covs)
  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-1000 datasets
# generated in (b), for S = 30... 

# 1) With initial guess (0, 0), using BFGS: 

results_30_1000_00 = map(1:100, ~ optim(c(0, 0), 
                                            alpha_gmm_minimizer_30_1000, col_index = .x,
                                            method = 'BFGS',
                                            control = list('maxit' = '1000')),
                                            )

# Checking convergence status: 
convergence_30_1000_00 = purrr::map(1:100, ~ results_30_1000_00[[.x]]$convergence)

# Results:
alpha_est_30_1000_00 = t(matrix(unlist(purrr::map(1:100, ~ results_30_1000_00[[.x]]$par)), nrow = 2))

means_30_1000_00 <- colMeans(alpha_est_30_1000_00) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_1000_00 <- var(alpha_est_30_1000_00[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_1000_00 <- var(alpha_est_30_1000_00[,2]) # Variance of alpha_1 estimates. 
msqe_alpha0_30_1000_00 <- mean((alpha_est_30_1000_00[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_1000_00 <- mean((alpha_est_30_1000_00[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_1000_00
var_alpha0_30_1000_00
var_alpha1_30_1000_00
msqe_alpha0_30_1000_00
msqe_alpha1_30_1000_00

# 2) With initial guess (0.5, 0.5), using BFGS: 

results_30_1000_55 = map(1:100, ~ optim(c(0.5, 0.5), 
                                               alpha_gmm_minimizer_30_1000, col_index = .x,
                                               method = 'BFGS',
                                               control = list('maxit' = '1000')),
                                )

# Checking convergence status: 
convergence_30_1000_55 = purrr::map(1:100, ~ results_30_1000_55[[.x]]$convergence)

# Results:
alpha_est_30_1000_55 = t(matrix(unlist(purrr::map(1:100, ~ results_30_1000_55[[.x]]$par)), nrow = 2))

means_30_1000_55 <- colMeans(alpha_est_30_1000_55) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_1000_55 <- var(alpha_est_30_1000_55[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_1000_55 <- var(alpha_est_30_1000_55[,2]) # Variance of alpha_1 estimates. 

msqe_alpha0_30_1000_55 <- mean((alpha_est_30_1000_55[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_1000_55 <- mean((alpha_est_30_1000_55[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_1000_55
var_alpha0_30_1000_55
var_alpha1_30_1000_55
msqe_alpha0_30_1000_55
msqe_alpha1_30_1000_55

# 3) With initial guess (0.75, 0.75), using BFGS: 

results_30_1000_7575 = map(1:100, ~ optim(c(0.75, 0.75), 
                                               alpha_gmm_minimizer_30_1000, col_index = .x,
                                               method = 'BFGS',
                                               control = list('maxit' = '1000')),
                                )

# Checking convergence status: 
convergence_30_1000_7575 = purrr::map(1:100, ~ results_30_1000_7575[[.x]]$convergence)

# Results:
alpha_est_30_1000_7575 = t(matrix(unlist(purrr::map(1:100, ~ results_30_1000_7575[[.x]]$par)), nrow = 2))

means_30_1000_7575 <- colMeans(alpha_est_30_1000_7575) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_30_1000_7575 <- var(alpha_est_30_1000_7575[,1]) # Variance of alpha_0 estimates.
var_alpha1_30_1000_7575 <- var(alpha_est_30_1000_7575[,2]) # Variance of alpha_1 estimates. 

msqe_alpha0_30_1000_7575 <- mean((alpha_est_30_1000_7575[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_30_1000_7575 <- mean((alpha_est_30_1000_7575[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_30_1000_7575
var_alpha0_30_1000_7575
var_alpha1_30_1000_7575
msqe_alpha0_30_1000_7575
msqe_alpha1_30_1000_7575

### (d)

### SIZE 500 DATASETS ESTIMATES, S = 100:

# Defining the SMM estimator function for datasets of size 500, with S = 100.
alpha_gmm_minimizer_100_500 = function(param, col_index) {
  sim_eps0 = rgumbel(100*500)
  sim_eps1 = rgumbel(100*500)
  
  estimated_choice = matrix(ifelse(param[1] + param[2]*X_mat_500[, col_index] + sim_eps0 >= sim_eps1, 1, 0), nrow = 500)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_500[, col_index]))

  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_500[col_index] - mean_of_means, 
                      data_cov_500[col_index] - mean_of_covs)
  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-500 datasets
# generated in (b), for S = 100... 

# 1) With initial guess (0, 0), using BFGS: 
results_100_500_00 = map(1:100, ~ optim(c(0, 0), 
                                          alpha_gmm_minimizer_100_500, col_index = .x,
                                          method = 'BFGS',
                                          control = list('maxit' = '1000')),
                                          )

# Checking convergence status:
convergence_100_500_00 = purrr::map(1:100, ~ results_100_500_00[[.x]]$convergence)

# Results: 
alpha_est_100_500_00 = t(matrix(unlist(purrr::map(1:100, ~ results_100_500_00[[.x]]$par)), nrow = 2))

means_100_500_00 <- colMeans(alpha_est_100_500_00) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_500_00 <- var(alpha_est_100_500_00[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_500_00 <- var(alpha_est_100_500_00[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_500_00 <- mean((alpha_est_100_500_00[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_500_00 <- mean((alpha_est_100_500_00[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_500_00
var_alpha0_100_500_00
var_alpha1_100_500_00
msqe_alpha0_100_500_00
msqe_alpha1_100_500_00

# 2) With initial guess (0.5, 0.5), using BFGS:
results_100_500_55 = map(1:100, ~ optim(c(0.5, 0.5), 
                                            alpha_gmm_minimizer_100_500, col_index = .x,
                                            method = 'BFGS',
                                            control = list('maxit' = '1000')),
                             )

# Checking convergence status:
convergence_100_500_55 = purrr::map(1:100, ~ results_100_500_55[[.x]]$convergence)

# Results: 
alpha_est_100_500_55 = t(matrix(unlist(purrr::map(1:100, ~ results_100_500_55[[.x]]$par)), nrow = 2))

means_100_500_55 <- colMeans(alpha_est_100_500_55) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_500_55 <- var(alpha_est_100_500_55[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_500_55 <- var(alpha_est_100_500_55[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_500_55 <- mean((alpha_est_100_500_55[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_500_55 <- mean((alpha_est_100_500_55[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_500_55
var_alpha0_100_500_55
var_alpha1_100_500_55
msqe_alpha0_100_500_55
msqe_alpha1_100_500_55

# 3) With initial guess (0.75, 0.75), using BFGS:
results_100_500_75 = map(1:100, ~ optim(c(0.75, 0.75), 
                                               alpha_gmm_minimizer_100_500, col_index = .x,
                                               method = 'BFGS',
                                               control = list('maxit' = '1000'))
                                )

# Checking convergence status:
convergence_100_500_75 = purrr::map(1:100, ~ results_100_500_75[[.x]]$convergence)

# Results: 
alpha_est_100_500_75 = t(matrix(unlist(purrr::map(1:100, ~ results_100_500_75[[.x]]$par)), nrow = 2))

means_100_500_75 <- colMeans(alpha_est_100_500_75) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_500_75 <- var(alpha_est_100_500_75[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_500_75 <- var(alpha_est_100_500_75[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_500_75 <- mean((alpha_est_100_500_75[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_500_75 <- mean((alpha_est_100_500_75[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_500_75
var_alpha0_100_500_75
var_alpha1_100_500_75
msqe_alpha0_100_500_75
msqe_alpha1_100_500_75

### SIZE 1000 DATASETS ESTIMATES, S = 100:

# Defining the SMM estimator function for datasets of size 1000, with S = 100.
alpha_gmm_minimizer_100_1000 = function(param, col_index) {
  sim_eps0 = rgumbel(100*1000)
  sim_eps1 = rgumbel(100*1000)
  
  estimated_choice = matrix(ifelse(param[1] + param[2]*X_mat_1000[, col_index] + sim_eps0 >= sim_eps1, 1, 0), nrow=1000)
  mean_of_sims = colMeans(estimated_choice)
  mean_of_means = mean(mean_of_sims)

  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_1000[, col_index]))
  mean_of_covs = mean(cov_of_sims_with_x)
  moment_diff_vec = c(data_y_mean_1000[col_index] - mean_of_means, 
                      data_cov_1000[col_index] - mean_of_covs)
  
  t(moment_diff_vec) %*% moment_diff_vec
}

# Optimizing with respect to (alpha_0, alpha_1) for each of the 100 size-1000 datasets
# generated in (b), for S = 100...

# 1) With initial guess (0, 0), using BFGS:

results_100_1000_00 = map(1:100, ~ optim(c(0, 0), 
                                          alpha_gmm_minimizer_100_1000, col_index = .x,
                                          method = 'BFGS',
                                          control = list('maxit' = '1000'))
                                          )

# Checking convergence status:
convergence_100_1000_00 = purrr::map(1:100, ~ results_100_1000_00[[.x]]$convergence)

# Results: 
alpha_est_100_1000_00 = t(matrix(unlist(purrr::map(1:100, ~ results_100_1000_00[[.x]]$par)), nrow = 2))

means_100_1000_00 <- colMeans(alpha_est_100_1000_00) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_1000_00 <- var(alpha_est_100_1000_00[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_1000_00 <- var(alpha_est_100_1000_00[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_1000_00 <- mean((alpha_est_100_1000_00[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_1000_00 <- mean((alpha_est_100_1000_00[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_1000_00
var_alpha0_100_1000_00
var_alpha1_100_1000_00
msqe_alpha0_100_1000_00
msqe_alpha1_100_1000_00

# 2) With initial guess (0.5, 0.5), using BFGS:

results_100_1000_55 = map(1:100, ~ optim(c(0.5, 0.5), 
                                                alpha_gmm_minimizer_100_1000, col_index = .x,
                                                method = 'BFGS',
                                                control = list('maxit' = '1000'))
                                 )

# Checking convergence status:
convergence_100_1000_55 = purrr::map(1:100, ~ results_100_1000_55[[.x]]$convergence)

# Results: 
alpha_est_100_1000_55 = t(matrix(unlist(purrr::map(1:100, ~ results_100_1000_55[[.x]]$par)), nrow = 2))

means_100_1000_55 <- colMeans(alpha_est_100_1000_55) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_1000_55 <- var(alpha_est_100_1000_55[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_1000_55 <- var(alpha_est_100_1000_55[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_1000_55 <- mean((alpha_est_100_1000_55[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_1000_55 <- mean((alpha_est_100_1000_55[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_1000_55
var_alpha0_100_1000_55
var_alpha1_100_1000_55
msqe_alpha0_100_1000_55
msqe_alpha1_100_1000_55

# 3) With initial guess (0.75, 0.75), using BFGS:

results_100_1000_7575 = map(1:100, ~ optim(c(0.75, 0.75), 
                                                alpha_gmm_minimizer_100_1000, col_index = .x,
                                                method = 'BFGS',
                                                control = list('maxit' = '1000')) 
                                 )

# Checking convergence status:
convergence_100_1000_7575 = purrr::map(1:100, ~ results_100_1000_7575[[.x]]$convergence)

# Results: 
alpha_est_100_1000_7575 = t(matrix(unlist(purrr::map(1:100, ~ results_100_1000_7575[[.x]]$par)), nrow = 2))

means_100_1000_7575 <- colMeans(alpha_est_100_1000_7575) # Mean of alpha_0 and alpha_1 estimates, respectively.
var_alpha0_100_1000_7575 <- var(alpha_est_100_1000_7575[,1]) # Variance of alpha_0 estimates.
var_alpha1_100_1000_7575 <- var(alpha_est_100_1000_7575[,2]) # Variance of alpha_1 estimates.
msqe_alpha0_100_1000_7575 <- mean((alpha_est_100_1000_7575[,1]-0.8)^2) # MSQE of alpha_0 estimates.
msqe_alpha1_100_1000_7575 <- mean((alpha_est_100_1000_7575[,2]-0.7)^2) # MSQE of alpha_1 estimates.

means_100_1000_7575
var_alpha0_100_1000_7575
var_alpha1_100_1000_7575
msqe_alpha0_100_1000_7575
msqe_alpha1_100_1000_7575





# Graphs of bootstrap estimations of SMM values given the last estimation
alpha_gmm_given_est_alpha1 = function(param, col_index) {
  sim_eps0 = rgumbel(100*1000)
  sim_eps1 = rgumbel(100*1000)
  
  estimated_choice = matrix(ifelse(param + colMeans(alpha_est_100_1000_7575)[2]*X_mat_1000[, col_index]
                                   + sim_eps0 >= sim_eps1, 1, 0), nrow = 1000)
  
  
  
  mean_of_sims = colMeans(estimated_choice)
  
  mean_of_means = mean(mean_of_sims)
  
  
  
  
  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_1000[, col_index]))
  
  mean_of_covs = mean(cov_of_sims_with_x)
  
  
  (data_y_mean_1000[col_index] - mean_of_means)^2 + (data_cov_1000[col_index] - mean_of_covs)^2
}


inner_map_alpha1 = function(param) {
  mean(map_dbl(1:100, ~ alpha_gmm_given_est_alpha1(param, .x)))
} 

graph_data_alpha1 = replicate(100, map_dbl(seq(0.1, 1.2, 0.05), inner_map_alpha1))

alpha_gmm_given_est_alpha0 = function(param, col_index) {
  sim_eps0 = rgumbel(100*1000)
  sim_eps1 = rgumbel(100*1000)
  
  estimated_choice = matrix(ifelse(colMeans(alpha_est_100_1000_7575)[1] + param*X_mat_1000[, col_index]
                                   + sim_eps0 >= sim_eps1, 1, 0), nrow = 1000)
  
  
  
  mean_of_sims = colMeans(estimated_choice)
  
  mean_of_means = mean(mean_of_sims)
  
  
  
  
  cov_of_sims_with_x = mapply(function(x, y){cov(x, y)}, as.data.frame(estimated_choice), 
                              as.data.frame(X_mat_1000[, col_index]))
  
  mean_of_covs = mean(cov_of_sims_with_x)
  
  
  (data_y_mean_1000[col_index] - mean_of_means)^2 + (data_cov_1000[col_index] - mean_of_covs)^2
}


inner_map_alpha0 = function(param) {
  mean(map_dbl(1:100, ~ alpha_gmm_given_est_alpha0(param, .x)))
} 

graph_data_alpha0 = replicate(100, map_dbl(seq(0.6, .9, 0.025), inner_map_alpha0))

library(ggplot2)
library(dplyr)


data.frame('Constant' = seq(0.1, 1.2, 0.05), 'upper' = apply(t(graph_data_alpha1), 2, max), 
           'lower' = apply(t(graph_data_alpha1), 2, min),
           'middle' = apply(t(graph_data_alpha1), 2, mean),
           check.names = FALSE) %>%
  ggplot(aes(x = Constant, y = middle)) + geom_line()+ 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)

data.frame('Constant' = seq(0.6, .9, 0.025), 'upper' = apply(t(graph_data_alpha0), 2, max), 
           'lower' = apply(t(graph_data_alpha0), 2, min),
           'middle' = apply(t(graph_data_alpha0), 2, mean),
           check.names = FALSE) %>%
  ggplot(aes(x = Constant, y = middle)) + geom_line()+ 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)


# Graph of sensitivity: plotting the estimated parameter for S=100 and n=1000 given the starting position
# in a grid

set.seed('2112')

estimator_with_100_sims = function(param) {
  results_100_1000 = map(1:100, ~ optim(param, alpha_gmm_minimizer_100_1000, col_index = .x,
                                              method = 'BFGS',
                                              control = list('maxit' = '1000')))
  
  convergence_100_1000 = purrr::map(1:100, ~ results_100_1000[[.x]]$convergence)
  
  alpha_est_100_1000 = t(matrix(unlist(purrr::map(1:100, ~ results_100_1000[[.x]]$par)), nrow = 2))
  
  return(c(colMeans(alpha_est_100_1000)))
}

end_locations = purrr::map2(expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02))[,1],
                            expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02))[,2],
                            ~ estimator_with_100_sims(c(.x, .y)))

df = cbind(expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02)), 
           t(matrix(unlist(end_locations), nrow = 3))) %>%
  rename(xend = `1`, yend = `2`)


ggplot(df) + geom_segment(aes(x=a0, y = a1, xend=xend, yend = yend),
                          arrow = arrow(length = unit(0.1, "cm"))) +
  xlim(0.5, 1) + ylim(0.5, 1) 