library(purrr)
data_ps1_q6 <- readRDS('data_ps1_q6.Rds')
data_ps1_q7 <- readRDS('data_ps1_q7.Rds')

######
# Q6 #
######

### (!) We should consider redoing this without applying logs 
### and using E[u] instead of E[u | X] as assumption.

# a)
# We start from  y = k^a*l^(1-a); then take logs
# given ln y = a*ln(k) + (1-a)*ln(l) + e, we need only to swap a few terms to find g(param)
# note that we take the hypothesis E(e|x)=0 on the ln form, not the original one, as it wouldn't work
# otherwise.

mme_synthetic_stat = function(param) {             #function that calculates the sample value of g
  sum(data_ps1_q6$L*(log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))
}

mme_synthetic_stat_2 = function(param) {             #function that calculates the sample value of g
  sum((log(data_ps1_q6$Y) - param*log(data_ps1_q6$K) -
                       (1-param)*log(data_ps1_q6$L)))
}



# b)
# We use map to apply the function over every value of param that is requested, then plot the curve.

x = seq(0, 1, 0.01)
fitted_curve =  purrr::map_dbl(x, mme_synthetic_stat)
fitted_curve_2 =  purrr::map_dbl(x, mme_synthetic_stat_2)

plot(fitted_curve_2)



# c)
# Finding the value of param that results in the closest value of g to 0.

index = match(min(abs(fitted_curve)), abs(fitted_curve))
index_2 = match(min(abs(fitted_curve_2)), abs(fitted_curve_2))
seq(0, 1, 0.01)[31]

# Labor share = 0.7; Capital share = 0.3.


######
# Q7 #
######


# Solution 1: proceeding as in Q6. 

# Now, there are 2 parameters, therefore we need two equations. Using the same hypothesis
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


result_q7 = optim(c(0,0), mme_statistic, method = "BFGS")$par #BFGS dando erro?

# Solution 2: solving using FOCs from firm's profit maximization problem.

focs_mme_condition_vector = function(param) {
  
  c_1 = mean(param[1] * data_ps1_q7$p * data_ps1_q7$K^(param[1] - 1) * data_ps1_q7$L^(param[2]) - data_ps1_q7$r)
  c_2 = mean(param[2] * data_ps1_q7$p * data_ps1_q7$K^(param[1]) * data_ps1_q7$L^(param[2] - 1) - data_ps1_q7$w)

  return(c(c_1, c_2))  
}

focs_mme_statistic = function(param) {
  100*t(focs_mme_condition_vector(param)) %*% diag(2) %*% focs_mme_condition_vector(param)
} # Pro BFGS functionar precisei definir a fç objetivo usando g' I g. Usando norma dava erro. (???)
# Mas o resultado deu igual ao resultado usando normal sem usar method = "BFGS" (i.e., usando o método padrão do optim)

focs_result_q7 = optim(c(0,0), focs_mme_statistic, method = "BFGS")$par
focs_result_q7 # Results: (alpha, beta) = (0.2, 0.6)



######
# Q8 #
######

result_q8 = optim(c(1,1), mme_statistic, method = "BFGS")$par
# got them, works pretty well cuz globally concave

focs_result_q8 = optim(c(1,1), focs_mme_statistic, method = "BFGS")$par  
# It changed dramatically! Multiple critical points. Depending on the initial guess, we hit
# the "wrong" (undesired) point: a local minimum, instead of a global minimum.

print(result_q8)
print(focs_result_q8)

######
# Q9 #
######

# We create the analytical jacobian by hand for the theoretical limit, and then apply the finite
# sample approximation, following the steps outlined on 6.3.4 of the book.

# Obs.: Without "focs_": stands for variables defined using the estimation done without FOCs.
# "focs_" stands for variables defined using the estimation done considering FOCs moment conditions.

# Matrix entries (calculated by hand)

g11 = mean(-data_ps1_q7$K*log(data_ps1_q7$K))
g12 = mean(-data_ps1_q7$K*log(data_ps1_q7$L))
g21 = mean(-data_ps1_q7$K*log(data_ps1_q7$L))
g22 = mean(-data_ps1_q7$L*log(data_ps1_q7$L))
G = matrix(c(g11, g12, g21, g22), ncol = 2)

focs_g11 = mean(data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]-1) * data_ps1_q7$L^(focs_result_q7[2])*(1 + focs_result_q7[1] * log(data_ps1_q7$K)))
focs_g12 = mean(focs_result_q7[1] * data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1] - 1) * data_ps1_q7$L^(focs_result_q7[2]) * log(data_ps1_q7$L))
focs_g21 = mean(focs_result_q7[2] * data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]) * log(data_ps1_q7$K) * data_ps1_q7$L^(focs_result_q7[2]-1))
focs_g22 = mean(data_ps1_q7$p * data_ps1_q7$K^(focs_result_q7[1]) * data_ps1_q7$L^(focs_result_q7[2]-1) * (1 + focs_result_q7[2]*log(data_ps1_q7$L)))
focs_G = matrix(c(focs_g11, focs_g12, focs_g21, focs_g22), ncol = 2)

calculate_h = function(data) { #generates the fitted value
  h_1 = data$K*(log(data$Y) - result_q8[1]*log(data$K) - result_q8[2]*log(data$L))
  h_2 = data$L*(log(data$Y) - result_q8[1]*log(data$K) - result_q8[2]*log(data$L))
  
  return(c(h_1, h_2))
}

focs_calculate_h = function(data) { #generates the fitted value

  h_1 = focs_result_q7[1] * data$p * data$K^(focs_result_q7[1] - 1) * data$L^(focs_result_q7[2]) - data$r
  h_2 = focs_result_q7[2] * data$p * data$K^(focs_result_q7[1]) * data$L^(focs_result_q7[2] - 1) - data$w
  
  return(c(h_1, h_2))
}

h_at_fit = data.frame(calculate_h(data_ps1_q7))
focs_h_at_fit = data.frame(focs_calculate_h(data_ps1_q7))

s11 = mean(h_at_fit[1]^2) 
s12 = mean(h_at_fit[1] * h_at_fit[2])
s21 = s12
s22 = mean(h_at_fit[2]^2)

focs_s11 = sum(focs_h_at_fit[1]^2)/100 
focs_s12 = sum(focs_h_at_fit[1] * focs_h_at_fit[2])/100
focs_s21 = focs_s12
focs_s22 = sum(focs_h_at_fit[2]^2)/100

S = matrix(c(s11, s12, s21, s22), ncol = 2)
focs_S = matrix(c(focs_s11, focs_s12, focs_s21, focs_s22), ncol = 2)

# Variance formula
V = 1/100*solve(t(G) %*% G) %*% t(G) %*% S %*% G %*% solve(t(G) %*% G)
focs_V = 1/100*solve(t(focs_G) %*% focs_G) %*% t(focs_G) %*% focs_S %*% focs_G %*% solve(t(focs_G) %*% focs_G)

print(V)
print(focs_V)

#######
# Q10 #
#######

# Doing as outlined in 6.3.5.

W = solve(S)
focs_W = solve(focs_S)

value_function = function(param) {
  h1 = mean(data_ps1_q7$K * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h2 = mean(data_ps1_q7$L * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h = c(h1, h2)
  
  Qn = t(h) %*% W %*% h
  
  return(Qn)
}

focs_value_function = function(param) {
  h1 = mean(data_ps1_q7$K * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h2 = mean(data_ps1_q7$L * (log(data_ps1_q7$Y) - param[1]*log(data_ps1_q7$K)
                             - param[2]*log(data_ps1_q7$L)))
  
  h = c(h1, h2)
  
  Qn = t(h) %*% focs_W %*% h
  
  return(Qn)
}


result_q10 = optim(c(0,0), value_function, method = "BFGS")$par
focs_result_q10 = optim(c(0,0), focs_value_function, method = "BFGS")$par

print(result_q10)
print(focs_result_q10)
