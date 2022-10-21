library(tidyverse)
library(ggplot2)
library(purrr)

data = readRDS('PS2/adpw_2017_marathon.RDS')

############## q6

histogram_point_estimator = function(point, bandwidth = 5) {
  height = sum(abs(data$chiptime-point)/bandwidth < 1)
  n = nrow(data)
  
  return(height/(n*bandwidth*2))
}

x0.grid=seq(105,355,10)

q6_hist_point_estimate = map_dbl(x0.grid, histogram_point_estimator)
q6_aux_tibble = tibble(x0.grid, q6_hist_point_estimate)

q6_aux_tibble %>% ggplot(aes(x=x0.grid, y = q6_hist_point_estimate)) + geom_col(fill = 'white', color = 'gray') +
  geom_segment(aes(x=240, xend=240, y=0, yend=0.008)) +
  geom_segment(aes(x=270, xend=270, y=0, yend=0.008)) +
  geom_segment(aes(x=300, xend=300, y=0, yend=0.008)) + 
  xlab('Completion time') + ylab('')




############## q7

opt_h = (1/nrow(data))^(1/5)

quartic_function = function(data_point, central_point, bandwidth = opt_h) {
  z = (data_point - central_point)/bandwidth
  indicator =  abs(z) < 1
  
  return(15/16*(1-z^2)^2*indicator)
}

density_quartic_point_estimator = function(point, bandwidth = opt_h) {
  weigth = sum(map2_dbl(data$chiptime, point, quartic_function))
  
  return(1/(nrow(data)*bandwidth)*weigth)
}

x0.grid=seq(105,355,5)

q7_density_estimate = map_dbl(x0.grid, density_quartic_point_estimator)

q7_aux_tibble = tibble(x0.grid, q7_density_estimate)

q7_aux_tibble %>% ggplot(aes(x=x0.grid, y = q7_density_estimate)) + geom_col(fill = 'white', color = 'gray') +
  geom_segment(aes(x=240, xend=240, y=0, yend=0.008)) +
  geom_segment(aes(x=270, xend=270, y=0, yend=0.008)) +
  geom_segment(aes(x=300, xend=300, y=0, yend=0.008)) + 
  xlab('Completion time') + ylab('')







############## q8


q3_bandwidth = 2.778*(nrow(data))^(-1/5)*min(sd(data$chiptime), IQR(data$chiptime)/1.349)


q8_quartic_function = function(data_point, central_point, bandwidth = q3_bandwidth) {
  z = (data_point - central_point)/bandwidth
  indicator =  abs(z) < 1
  
  return(15/16*(1-z^2)^2*indicator)
}

q8_density_quartic_point_estimator = function(point, bandwidth = q3_bandwidth) {
  weigth = sum(map2_dbl(data$chiptime, point, q8_quartic_function))
  
  return(1/(nrow(data)*bandwidth)*weigth)
}

q8_density_estimate = map_dbl(x0.grid, q8_density_quartic_point_estimator)

q8_aux_tibble = tibble(x0.grid, q8_density_estimate)

q8_aux_tibble %>% ggplot(aes(x=x0.grid, y = q8_density_estimate)) + geom_col(fill = 'white', color = 'gray') +
  geom_segment(aes(x=240, xend=240, y=0, yend=0.008)) +
  geom_segment(aes(x=270, xend=270, y=0, yend=0.008)) +
  geom_segment(aes(x=300, xend=300, y=0, yend=0.008)) + 
  xlab('Completion time') + ylab('')


                                                




################### q9

q9_data = dplyr::filter(data, !is.na(age))


local_constant_estimator = function(point) {
  kernel_stat = map_dbl(q9_data$age, ~ quartic_function(.x, point, bandwidth = 3))
  
  return(sum(kernel_stat*q9_data$chiptime)/sum(kernel_stat))
}

q9.grid = seq(20, 70, 5)

q9_regression_estimate = map_dbl(q9.grid, local_constant_estimator)
q9_aux_tibble = tibble(q9.grid, q9_regression_estimate)

q9_aux_tibble %>% ggplot(aes(x=q9.grid, y = q9_regression_estimate)) + geom_line() +
  xlab('Age') + ylab('Completion time') + theme_bw()






################## q10

z_matrix_generator = function(point) {
  ones = rep(1, nrow(q9_data))
  xdiff = q9_data$age-point

  return(data.frame(ones, xdiff))
}

left_side_term_calculator = function(point, k_vec, z_mat) {
  
  aux_matrix = cbind(k_vec, z_mat)
  
  aux_matrix2 = subset(aux_matrix, aux_matrix$k_vec != 0)
  
  map_value = map(1:nrow(aux_matrix2), ~ aux_matrix2[.x, 1] * 
                    t(aux_matrix2[.x, 2:3]) %*% as.matrix(aux_matrix2[.x, 2:3])) 
  
  return(Reduce('+', map_value))
}

right_hand_term_calculator = function(point, k_vec, z_mat) {
  
  return(colSums(k_vec*z_mat*q9_data$chiptime))
}


local_linear_estimator = function(point) {
  k_vector = map_dbl(q9_data$age, ~ quartic_function(.x, point, bandwidth = 3))
  
  z_matrix = z_matrix_generator(point)
  
  lht = left_side_term_calculator(point, k_vector, z_matrix)
  
  rht = right_hand_term_calculator(point, k_vector, z_matrix)
  
  
  return(solve(lht) %*% rht)
}


system.time(beta_values <- map(q9.grid, local_linear_estimator))



#### visualization

q10_graph_tibble = tibble(as.data.frame(matrix(unlist(beta_values), ncol = 2, byrow = TRUE)), age = q9.grid) %>%
  mutate(lowery = V1-3*V2,
         lowerx = age - 3,
         uppery = V1+3*V2,
         upperx = age + 3)

q10_graph_tibble %>% ggplot(aes(x=age, y=V1)) + geom_line(alpha = 0.5, size = 1, linetype = 2) +
  geom_segment(aes(x=lowerx, xend=upperx, y=lowery, yend = uppery), size = 1.1) + theme_bw()
