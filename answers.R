
#q6

#a)
# we start from  y = k^a*l^(1-a); then take logs
# given ln y = a*ln(k) + (1-a)*ln(l) + e, we need only to swap a few terms to find g(param)
# note that we take the hypothesis E(e|x)=0 on the ln form, not the original one, as it wouldn't work
# otherwise.

mme_synthetic_stat = function(param) {             #function that calculates the sample value of g
  sum(data_ps1_q6$L*(log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))
}

#b)
# we use map to apply the function over every value of param that is requested, then plot the curve

fitted_curve =  purrr::map_dbl(seq(0, 1, 0.01), mme_synthetic_stat)

plot(fitted_curve)

#c)
#finding the value of param that results in the closest value of g to 0.

index = match(min(abs(fitted_curve)), abs(fitted_curve))

seq(0, 1, 0.01)[31]

#labor share = 0.7, capital share = 0.3




#q7
# now, there are 2 parameters, therefore we need two equations. Using the same hypothesis
# as in q6, note that both K and L are in x, therefore we may aswell generate two equations
# in a similar manner.

mme_condition_vector = function(param) {                          #analogous function for two param
  c1 = sum(data_ps1_q7$L*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  c2 = sum(data_ps1_q7$K*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  return(c(c1, c2))
}



mme_statistic = function(param) {                      #calculates the value of the stat for W = I
  100*norm(mme_condition_vector(param), type = "2")^2
}


result_q7 = optim(c(0,0), mme_statistic, method = "BFGS")$par



#q8

result_q8 = optim(c(1,1), mme_statistic, method = "BFGS")$par

# got them, works pretty well cuz globally concave



#q9
# we create the analytical jacobian by hand for the theoretical limit, and then apply the finite
# sample approximation, following the steps outlined on 6.3.4 of the book

#matrix entries by hand
g11 = mean(-data_ps1_q7$K*log(data_ps1_q7$K))
g12 = mean(-data_ps1_q7$K*log(data_ps1_q7$L))
g21 = mean(-data_ps1_q7$K*log(data_ps1_q7$L))
g22 = mean(-data_ps1_q7$L*log(data_ps1_q7$L))
G = matrix(c(g11, g12, g21, g22), ncol = 2)

calculate_h = function(data) {                                           #generates the fitted value
  h_1 = data[1]*(log(data[3]) - result_q8[1]*log(data[1]) - result_q8[2]*log(data[2]))
  h_2 = data[2]*(log(data[3]) - result_q8[1]*log(data[1]) - result_q8[2]*log(data[2]))
  
  return(c(h_1, h_2))
}

s11 = mean(h_at_fit$K^2)
s12 = mean(h_at_fit$K * h_at_fit$L)
s21 = s12
s22 = mean(h_at_fit$L^2)

S = matrix(c(s11, s12, s21, s22), ncol = 2)

#variance formula
V = 1/100*solve(t(G) %*% G) %*% t(G) %*% S %*% G %*% solve(t(G) %*% G)



#q10

# doing as outlined in 6.3.5

W = solve(S)

value_function = function(param) {
  h1 = mean(data_ps1_q7$K * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h2 = mean(data_ps1_q7$L * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h = c(h1, h2)
  
  Qn = t(h) %*% W %*% h
  
  return(Qn)
}


result_q10 = optim(c(0,0), value_function, method = "BFGS")$par
