three_eq_mme_condition_vector = function(param) {
  
  c_1 = mean(data_ps1_q7$Y - data_ps1_q7$K^param[1]*data_ps1_q7$L^param[2])
  c_2 = mean(param[1] * data_ps1_q7$p * data_ps1_q7$K^(param[1] - 1) * data_ps1_q7$L^(param[2]) - data_ps1_q7$r)
  c_3 = mean(param[2] * data_ps1_q7$p * data_ps1_q7$K^(param[1]) * data_ps1_q7$L^(param[2] - 1) - data_ps1_q7$w)
  
  return(c(c_1, c_2, c_3))  
}

three_eq_mme_statistic = function(param) {
  100*t(three_eq_mme_condition_vector(param)) %*% diag(3) %*% three_eq_mme_condition_vector(param)
}


three_eq_result_q7 = optim(c(0,0), three_eq_mme_statistic, method = "BFGS")$par

three_eq_result_q8 = optim(c(1,1), three_eq_mme_statistic, method = "BFGS")$par


three_eq_calculate_h = function(data) { #generates the fitted value
  
  h_1 = data$Y - data$K^three_eq_result_q7[1]*data$L^three_eq_result_q7[2]
  h_2 = three_eq_result_q7[1] * data$p * data$K^(three_eq_result_q7[1] - 1) * data$L^(three_eq_result_q7[2]) - data$r
  h_3 = three_eq_result_q7[2] * data$p * data$K^(three_eq_result_q7[1]) * data$L^(three_eq_result_q7[2] - 1) - data$w
  
  return(data.frame(h_1, h_2, h_3))
}

three_h_at_fit = three_eq_calculate_h(data_ps1_q7)

three_h_S = t(as.matrix(three_h_at_fit)) %*% as.matrix(three_h_at_fit)/100

three_g11 = mean(-three_eq_result_q7[1]*data_ps1_q7$K^(three_eq_result_q7[1]-1) * data_ps1_q7$L^(three_eq_result_q7[2]))
three_g12 = mean(-three_eq_result_q7[2]* data_ps1_q7$K^(three_eq_result_q7[1]) * data_ps1_q7$L^(three_eq_result_q7[2]-1))
three_g21 = mean(data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]-1) * data_ps1_q7$L^(three_eq_result_q7[2])*(1 + three_eq_result_q7[1] * log(data_ps1_q7$K)))
three_g22 = mean(three_eq_result_q7[1] * data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1] - 1) * data_ps1_q7$L^(three_eq_result_q7[2]) * log(data_ps1_q7$L))
three_g31 = mean(three_eq_result_q7[2] * data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]) * log(data_ps1_q7$K) * data_ps1_q7$L^(three_eq_result_q7[2]-1))
three_g32 = mean(data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]) * data_ps1_q7$L^(three_eq_result_q7[2]-1) * (1 + three_eq_result_q7[2]*log(data_ps1_q7$L)))
three_G = matrix(c(three_g11, three_g12, three_g21, three_g22, three_g31, three_g32), ncol = 2)

three_G
three_V = 1/100*solve(t(three_G) %*% three_G) %*% t(three_G) %*% three_h_S %*% three_G %*% solve(t(three_G) %*% three_G)
