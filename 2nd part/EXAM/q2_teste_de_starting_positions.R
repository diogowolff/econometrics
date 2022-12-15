set.seed(1337)

library(evd)
library(furrr)
library(ggplot2)

X_mat_1000 = matrix(rlnorm(100*1000, log(25/sqrt(25+3)), log(1+3/25)), nrow = 1000)

eps0_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)
eps1_mat_1000 = matrix(rgumbel(100*1000), nrow = 1000)

Y_mat_1000 = ifelse(0.8 + 0.7*X_mat_1000 + eps1_mat_1000 >= eps0_mat_1000, 1, 0)

data_y_mean_1000 = colMeans(Y_mat_1000)
data_cov_1000 = mapply(function(x, y){cov(x, y)}, as.data.frame(Y_mat_1000), 
                       as.data.frame(X_mat_1000))

alpha_gmm_minimizer_20_1000 = function(param, col_index) {
  data_x = matrix(rep(X_mat_1000[, col_index], 20), nrow = 1000)
  sim_eps0 = matrix(rgumbel(20*1000), nrow = 1000)
  sim_eps1 = matrix(rgumbel(20*1000), nrow = 1000)
  
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


estimator_with_20_sims = function(param) {
  results_20_1000 = future_map(1:100, ~ optim(param, alpha_gmm_minimizer_20_1000, col_index = .x,
                                               method = 'BFGS',
                                               control = list('maxit' = '1000')), 
                                .options = furrr_options(seed = T))
  
  convergence_20_1000 = purrr::map(1:100, ~ results_20_1000[[.x]]$convergence)
  
  alpha_est_100_1000 = t(matrix(unlist(purrr::map(1:10, ~ results_20_1000[[.x]]$par)), nrow = 2))
  
  return(c(colMeans(alpha_est_100_1000), sum(unlist(convergence_20_1000))))
}

end_locations = purrr::map2(expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02))[,1],
                            expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02))[,2],
            ~ estimator_with_20_sims(c(.x, .y)))
            
df = cbind(expand.grid(a0 = seq(0.6, 1, .02), a1 = seq(.5, .9, .02)), 
           t(matrix(unlist(end_locations), nrow = 3))) %>%
  rename(xend = `1`, yend = `2`)


ggplot(df) + geom_segment(aes(x=a0, y = a1, xend=xend, yend = yend),
                          arrow = arrow(length = unit(0.1, "cm"))) +
  xlim(0.5, 1) + ylim(0.5, 1)


