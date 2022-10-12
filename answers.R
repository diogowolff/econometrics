###########################################
# Econometrics - Problem Set 1            #
# Students: Luan Borelli and Diogo Wolff  #
###########################################


# Importing useful libraries and data.

library(purrr)
library(ggplot2)
data_ps1_q6 <- readRDS('data_ps1_q6.Rds')
data_ps1_q7 <- readRDS('data_ps1_q7.Rds')

######
# Q6 #
######

# a) As discussed in the LaTeX document (check it!), we consider several different specifications for the estimation.
# All of them yield the same results.

## Without log-linearizing 

mme_function = function(param) {
  sum(data_ps1_q6$Y - data_ps1_q6$K^(param) * data_ps1_q6$L^(1-param))/100
}


## Log-linearizing 

mme_function_loglin = function(param) {             #function that calculates the sample value of g
  sum((log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
         (1-param)*log(data_ps1_q6$L)))/100
}

mme_function_loglin_K = function(param) {             #function that calculates the sample value of g
  sum(data_ps1_q6$L*(log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))/100
}

mme_function_loglin_L = function(param) {             #function that calculates the sample value of g
  sum(data_ps1_q6$L*(log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))/100
}



# b)
# We use map to apply the function over every value of param that is requested, then plot the curve.

x = seq(0, 1, 0.01)
fitted_curve =  purrr::map_dbl(x, mme_function)
fitted_curve_loglin =  purrr::map_dbl(x, mme_function_loglin)
fitted_curve_loglin_K = purrr::map_dbl(x, mme_function_loglin_K)
fitted_curve_loglin_L = purrr::map_dbl(x, mme_function_loglin_L)

# Plotting curves: 

# Non-linearized version, unconditional expectation hypothesis.  

ggplot(data = data.frame(x, fitted_curve), aes(x=x, y=fitted_curve)) +
  geom_line(size = 1.2) +
  labs(x = "Grid", y  = "MME moment function, w/o log-linearization")

# Linearized version, conditional expectation hypothesis. 
ggplot(data = data.frame(x, fitted_curve_loglin), aes(x=x, y=fitted_curve_loglin)) +
  geom_line(size = 1.2) +
  labs(x = "Grid", y  = "MME moment function, log-linearized")

ggplot(data = data.frame(x, fitted_curve_loglin_K), aes(x=x, y=fitted_curve_loglin_K)) +
  geom_line(size = 1.2) +
  labs(x = "Grid", y  = "MME moment function, log-linearized, K")

ggplot(data = data.frame(x, fitted_curve_loglin_L), aes(x=x, y=fitted_curve_loglin_L)) +
  geom_line(size = 1.2) +
  labs(x = "Grid", y  = "MME moment function, log-linearized, L")


# c)
# Finding the value of param that results in the closest value of g to 0 for each specification.

index = match(min(abs(fitted_curve)), abs(fitted_curve))
index_loglin = match(min(abs(fitted_curve_loglin)), abs(fitted_curve_loglin))
index_loglin_K = match(min(abs(fitted_curve_loglin_K)), abs(fitted_curve_loglin_K))
index_loglin_L = match(min(abs(fitted_curve_loglin_L)), abs(fitted_curve_loglin_L))

value = x[index]
value_loglin = x[index_loglin]
value_loglin_K = x[index_loglin_K]
value_loglin_L = x[index_loglin_L]

print(c('Non-linearized result:', c(value, 1-value)))
print(c('Linearized result, standard:', c(value_loglin, 1-value_loglin)))
print(c('Linearized result, condition on K:', c(value_loglin_K, 1-value_loglin_K)))
print(c('Linearized result, condition on L:', c(value_loglin_L, 1-value_loglin_L)))

# Capital share (alpha) = 0.3; Labor share (1-alpha) = 0.7
# for all cases. All the different specifications ended up being equivalent.

######
# Q7 #
######


## As discussed in the LaTeX document (check it!), we've studied three different approaches. 

# The first approach, which we call the "naive approach" requires only the assumption that the production 
# of the economy follows a Cobb-Douglas production function. The shortcoming of this approach is that it 
# discards all available price data (w, p and r) from the estimation. The second, which we call the 
# FOCs approach, requires some economic theory by assuming that firms behave as profit-maximizers. 
# In this way, price data are used in the estimation. However, we continue to discard some data: 
# output data is not considered in this approach. Finally, the third, which we call the "full approach", 
# gathers the two previous approaches and set an over-identified system of three moment
# equations -- although only two equations were required, as specified by the exercise. 
# This last approach is especially attractive, as it is the only one that incorporates the 
# entire available dataset in the estimation procedure. In addition, it also motivates the use of GMM 
# estimation, since this method is attractive precisely for the case of over-identified systems.


## NAIVE APPROACH

# We work here with the log-linearized version presented in Q6, but we 
# could consider the non-linearized version as well. The results would be equivalent.
 
# Now, there are 2 parameters, therefore we need two equations. Using the same hypothesis
# as in Q6, conditioning on K and L we are able to generate two moment equations.
# Check LaTeX document for more details.

mme_condition_vector = function(param) {                          #analogous function for two param
  c1 = sum(data_ps1_q7$L*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  c2 = sum(data_ps1_q7$K*(log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K) -
                            param[2]*log(data_ps1_q7$L)))
  
  return(c(c1, c2))
}


mme_statistic = function(param) { # This calculates the value of the stat for W = I.
  100*norm(mme_condition_vector(param), type = "2")^2
}


result_q7 = optim(c(0,0), mme_statistic, method = "BFGS")$par

print(c("(alpha, beta) = ", result_q7))
# Results: (alpha, beta) = (0.3532129, 2.8355161) 

## FOCS APPROACH: solving using FOCs from firm's profit maximization problem.
# Check LaTeX document for details. 

focs_mme_condition_vector = function(param) {
  
  c_1 = mean(param[1] * data_ps1_q7$p * data_ps1_q7$K^(param[1] - 1) * data_ps1_q7$L^(param[2]) - data_ps1_q7$r)
  c_2 = mean(param[2] * data_ps1_q7$p * data_ps1_q7$K^(param[1]) * data_ps1_q7$L^(param[2] - 1) - data_ps1_q7$w)

  return(c(c_1, c_2))  
}

focs_mme_statistic = function(param) {
  100*t(focs_mme_condition_vector(param)) %*% focs_mme_condition_vector(param)
} 

focs_result_q7 = optim(c(0,0), focs_mme_statistic, method = "BFGS")$par

print(c("(alpha, beta) =", focs_result_q7))

# Results: (alpha, beta) = (0.2, 0.6)

## FULL APPROACH: Using three moment conditions.

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

print(c("(alpha, beta) = ", three_eq_result_q7))
# Results: (alpha, beta) = (0.316859723130909, 0.859379092580883)

######
# Q8 #
######

result_q8 = optim(c(1,1), mme_statistic, method = "BFGS")$par
print(result_q8)
# Pretty consistent. The function is globally concave.

focs_result_q8 = optim(c(1,1), focs_mme_statistic, method = "BFGS")$par  
print(focs_result_q8)
# It changed dramatically! Multiple critical points. Depending on the initial guess, we hit
# the "wrong" (undesired) point: a local minimum, instead of a global minimum.

three_eq_result_q8 = optim(c(1,1), three_eq_mme_statistic, method = "BFGS")$par
print(three_eq_result_q8)
# Also pretty consistent.


######
# Q9 #
######

# We create the analytical jacobian by hand for the theoretical limit, and then apply the finite
# sample approximation.
# Check LaTeX document for further details. 


## NAIVE APPROACH:

# Jacobian entries (calculated analytically, by hand):
g11 = mean(-data_ps1_q7$K*log(data_ps1_q7$K))
g12 = mean(-data_ps1_q7$K*log(data_ps1_q7$L))
g21 = mean(-data_ps1_q7$L*log(data_ps1_q7$K))
g22 = mean(-data_ps1_q7$L*log(data_ps1_q7$L))

# Defining the Jacobian:
G = matrix(c(g11, g12, g21, g22), ncol = 2)

calculate_h = function(data) { # This generates the fitted values.
  h_1 = data[1]*(log(data[3]) - result_q8[1]*log(data[1]) - result_q8[2]*log(data[2]))
  h_2 = data[2]*(log(data[3]) - result_q8[1]*log(data[1]) - result_q8[2]*log(data[2]))
  
  return(c(h_1, h_2))
}

# Fitted values of the moment conditions. 
h_at_fit = data.frame(calculate_h(data_ps1_q7[, 4:6]))


# Calculating S. 
S = t(as.matrix(h_at_fit)) %*% as.matrix(h_at_fit)/100

# Calculating the variance-covariance matrix: 
V = 1/100*solve(t(G) %*% G) %*% t(G) %*% S %*% G %*% solve(t(G) %*% G)

print("Jacobian - Naive approach:")
print(G)
print("S hat matrix - Naive approach:")
print(S)
print("Variance-covariance matrix - Naive approach:")
print(V)

## FOCS APPROACH: 

# Jacobian entries (calculated analytically, by hand):
focs_g11 = mean(data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]-1) * data_ps1_q7$L^(focs_result_q7[2])*(1 + focs_result_q7[1] * log(data_ps1_q7$K)))
focs_g12 = mean(focs_result_q7[1] * data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1] - 1) * data_ps1_q7$L^(focs_result_q7[2]) * log(data_ps1_q7$L))
focs_g21 = mean(focs_result_q7[2] * data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]) * log(data_ps1_q7$K) * data_ps1_q7$L^(focs_result_q7[2]-1))
focs_g22 = mean(data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]) * data_ps1_q7$L^(focs_result_q7[2]-1) * (1 + focs_result_q7[2]*log(data_ps1_q7$L)))

# Defining the Jacobian:
focs_G = matrix(c(focs_g11, focs_g12, focs_g21, focs_g22), ncol = 2)

focs_calculate_h = function(data) { # Generates the fitted value.
  
  h_1 = focs_result_q7[1] * data[3] * data[4]^(focs_result_q7[1] - 1) * data[5]^(focs_result_q7[2]) - data[1]
  h_2 = focs_result_q7[2] * data[3] * data[4]^(focs_result_q7[1]) * data[5]^(focs_result_q7[2] - 1) - data[2]
  
  return(c(h_1, h_2))
}

# Fitted values of the moment conditions: 
focs_h_at_fit = data.frame(focs_calculate_h(data_ps1_q7))

# Calculating S: 
focs_S = t(as.matrix(focs_h_at_fit)) %*% as.matrix(focs_h_at_fit)/100

# Calculating the variance-covariance matrix:
focs_V = 1/100*solve(t(focs_G) %*% focs_G) %*% t(focs_G) %*% focs_S %*% focs_G %*% solve(t(focs_G) %*% focs_G)

print("Jacobian - FOCs approach:")
print(focs_G)
print("S hat matrix - FoCs approach:")
print(focs_S)
print("Variance-covariance matrix - FOCs approach:")
print(focs_V)

## FULL APPROACH:

# Jacobian entries (calculated analytically, by hand):
three_g11 = mean(-three_eq_result_q7[1]*data_ps1_q7$K^(three_eq_result_q7[1]-1) * data_ps1_q7$L^(three_eq_result_q7[2]))
three_g12 = mean(-three_eq_result_q7[2]* data_ps1_q7$K^(three_eq_result_q7[1]) * data_ps1_q7$L^(three_eq_result_q7[2]-1))
three_g21 = mean(data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]-1) * data_ps1_q7$L^(three_eq_result_q7[2])*(1 + three_eq_result_q7[1] * log(data_ps1_q7$K)))
three_g22 = mean(three_eq_result_q7[1] * data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1] - 1) * data_ps1_q7$L^(three_eq_result_q7[2]) * log(data_ps1_q7$L))
three_g31 = mean(three_eq_result_q7[2] * data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]) * log(data_ps1_q7$K) * data_ps1_q7$L^(three_eq_result_q7[2]-1))
three_g32 = mean(data_ps1_q7$p * data_ps1_q7$K^(three_eq_result_q7[1]) * data_ps1_q7$L^(three_eq_result_q7[2]-1) * (1 + three_eq_result_q7[2]*log(data_ps1_q7$L)))

# Defining the Jacobian:
three_G = matrix(c(three_g11, three_g12, three_g21, three_g22, three_g31, three_g32), ncol = 2)

three_eq_calculate_h = function(data) { # Generates the fitted values.
  
  h_1 = data$Y - data$K^three_eq_result_q7[1]*data$L^three_eq_result_q7[2]
  h_2 = three_eq_result_q7[1] * data$p * data$K^(three_eq_result_q7[1] - 1) * data$L^(three_eq_result_q7[2]) - data$r
  h_3 = three_eq_result_q7[2] * data$p * data$K^(three_eq_result_q7[1]) * data$L^(three_eq_result_q7[2] - 1) - data$w
  
  return(data.frame(h_1, h_2, h_3))
}

# Fitted values of the moment conditions: 
three_h_at_fit = three_eq_calculate_h(data_ps1_q7)

# Calculating S: 
three_h_S = t(as.matrix(three_h_at_fit)) %*% as.matrix(three_h_at_fit)/100

# Calculating the variance-covariance matrix: 
three_V = 1/100*solve(t(three_G) %*% three_G) %*% t(three_G) %*% three_h_S %*% three_G %*% solve(t(three_G) %*% three_G)

print("Jacobian - Full approach:")
print(three_G)
print("S hat matrix - Full approach:")
print(three_h_S)
print("Variance-covariance matrix - Full approach:")
print(three_V)

#######
# Q10 #
#######

## NAIVE APPROACH

# Calculating the optimal weighting matrix.
W = solve(S)

value_function = function(param) { # Setting up the objetive function using the optimal weighting matrix.
  h1 = mean(data_ps1_q7$K * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h2 = mean(data_ps1_q7$L * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h = c(h1, h2)
  
  Qn = t(h) %*% W %*% h
  
  return(Qn)
}

# Minimizing the objective function. 
result_q10 = optim(c(0,0), value_function, method = "BFGS")$par

print(c("(alpha_OGMM, beta_OGMM) = ", result_q10, "(Naive approach)"))


## FOCS approach 

# Calculating the optimal weighting matrix.
focs_W = solve(focs_S)

focs_value_function = function(param) { # Setting up the objetive function using the optimal weighting matrix.
  
  h1 = mean(param[1] * data_ps1_q7$p * data_ps1_q7$K^(param[1] - 1) * data_ps1_q7$L^(param[2]) - data_ps1_q7$r)
  h2 = mean(param[2] * data_ps1_q7$p * data_ps1_q7$K^(param[1]) * data_ps1_q7$L^(param[2] - 1) - data_ps1_q7$w)
  
  h = c(h1, h2)
  
  Qn = t(h) %*% focs_W %*% h
  
  return(Qn)
}

# Minimizing the objective function. 
focs_result_q10 = optim(c(0,0), focs_value_function, method = "BFGS")$par

print(c("(alpha_OGMM, beta_OGMM) = ", focs_result_q10, "(FOCs approach)"))


## FULL APPROACH

# Faltando.






######################
# Summary of results #
######################


# Q6 

print(c('Non-linearized result:', c(value, 1-value)))
print(c('Linearized result, standard:', c(value_loglin, 1-value_loglin)))
print(c('Linearized result, condition on K:', c(value_loglin_K, 1-value_loglin_K)))
print(c('Linearized result, condition on L:', c(value_loglin_L, 1-value_loglin_L)))


# Q7 

print(c("(alpha, beta) = ", result_q7, "(Naive approach)"))
print(c("(alpha, beta) =", focs_result_q7, "(FOCs approach)"))
print(c("(alpha, beta) =", three_eq_result_q7, "(Full approach)"))

# Q8 

print(c(result_q8, "(Naive approach)"))
print(c(focs_result_q8, "(FOCs approach)"))
print(c(three_eq_result_q8, "(Full approach)"))

# Q9 

print("Jacobian - Naive approach:")
print(G)
print("S hat matrix - Naive approach:")
print(S)
print("Variance-covariance matrix - Naive approach:")
print(V)

print("Jacobian - FOCs approach:")
print(focs_G)
print("S hat matrix - FOCs approach:")
print(focs_S)
print("Variance-covariance matrix - FOCs approach:")
print(focs_V)

print("Jacobian - Full approach:")
print(three_G)
print("S hat matrix - Full approach:")
print(three_h_S)
print("Variance-covariance matrix - Full approach:")
print(three_V)

# Q10 

print(c("(alpha_OGMM, beta_OGMM) = ", result_q10, "(Naive approach)"))
print(c("(alpha_OGMM, beta_OGMM) = ", focs_result_q10, "(FOCs approach)"))

# Faltando full approach.