#q6

#a)
# given ln y = a*ln(k) + (1-a)*ln(l) + e, we need only to swap a few terms to find g(param)

mme_synthetic_stat = function(param) {
  sum(data_ps1_q6$L*(log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))
}

#b)

fitted_curve =  purrr::map_dbl(seq(0, 1, 0.01), mme_synthetic_stat)

plot(fitted_curve)

#c)

index = match(min(abs(fitted_curve)), abs(fitted_curve))

seq(0, 1, 0.01)[31]

#labor share = 0.7, capital share = 0.3




#q7

mme_condition_vector = function(param) {
  c1 = sum(data_ps1_q7$L*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  c2 = sum(data_ps1_q7$K*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  return(c(c1, c2))
}

mme_condition_vector(c(1, 2))

mme_statistic = function(param) {
  100*norm(mme_condition_vector(param), type = "2")^2
}


result_q7 = optim(c(0,0), mme_statistic, method = "BFGS")$par



#q8

result_q8 = optim(c(1,1), mme_statistic, method = "BFGS")$par

# got them, works pretty well cuz globally concave



#q9

