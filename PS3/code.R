###########################################
# Econometrics - Problem Set 2            #
# Students: Luan Borelli and Diogo Wolff  #
###########################################

############# Q6

# Defining Nelder-mead parameters:
reflection <- 1
expansion <- 1.7
contraction <- 0.3 
shrinkage <- 0.6

# Defining tolerance variables:
max_run <- 10000 # Loop max run times.
run_count <- 0 # Loop counter.
tolerance <- 10^(-6) # Tolerance.
tol <- 999999999999 # An initial tolerance value to be updated.

# Defining the Himmelblau function, to be minimized:
f <- function(theta) { 
  (theta[1]^2 + theta[2] - 11)^2 + (theta[1]+theta[2]^2 - 7)^2
}

# Defining the algorithm function.
# This function receives a list of three 2-dimensional vectors (i.e., a list of vertices)
# E.g.: list(c(2, 3), c(7, 3), c(6,4))

nelder_mead <- function(verts) {
  vertices <- verts
  while(run_count <= max_run & tol >= tolerance) {
    
    f_at_vertices <- sapply(vertices,f)
    best_pos <- which.min(f_at_vertices)
    worst_pos <- which.max(f_at_vertices)
    new_vertices <- vertices[-worst_pos]
    f_at_new_vertices <- sapply(new_vertices, f)
    
    centroid <- 1/length(new_vertices) * (unlist(new_vertices[1]) + unlist(new_vertices[2]))
    reflection_point <- (1+reflection)*centroid - reflection*unlist(vertices[worst_pos])
    
    if(f(reflection_point) < f_at_vertices[best_pos]) {
      theta_e = (1 + reflection*expansion)*centroid - reflection*expansion*unlist(vertices[worst_pos])
      if(f(theta_e) < f(reflection_point)) {
        vertices <- c(new_vertices, list(theta_e))
      } else {
        vertices <- c(new_vertices, list(reflection_point))
      }
    } else if(f(reflection_point) >= f_at_vertices[best_pos] & f(reflection_point) < f(unlist(new_vertices[which.max(f_at_new_vertices)]))) {
      vertices <- c(new_vertices, list(reflection_point))  
    } else if(f(unlist(new_vertices[which.max(f_at_new_vertices)])) <= f(reflection_point & f(reflection_point) < f_at_vertices[worst_pos])) {
      theta_c <- (1 + reflection*contraction)*centroid - reflection*contraction*unlist(vertices[worst_pos])
      if(f(theta_c) < f(reflection_point)) {
        vertices <- c(new_vertices, list(theta_c))  
      } else {
        tilde_theta_1 <- unlist(vertices[best_pos])
        tilde_theta_2 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(new_vertices[which.max(f_at_new_vertices)])
        tilde_theta_3 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(vertices[worst_pos])
        vertices <- list(tilde_theta_1, tilde_theta_2, tilde_theta_3)
      } 
    } else if(f(reflection_point) >= f_at_vertices[worst_pos]) {
      theta_c2 <- (1 - contraction)*centroid + contraction*unlist(vertices[worst_pos])
      if(f(theta_c2)<f_at_vertices[worst_pos]) {
        vertices <- c(new_vertices, list(theta_c2))
      } else {
        tilde_theta_1 <- unlist(vertices[best_pos])
        tilde_theta_2 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(new_vertices[which.max(f_at_new_vertices)])
        tilde_theta_3 <- shrinkage*unlist(vertices[best_pos]) + (1-shrinkage)*unlist(vertices[worst_pos])
        vertices <- list(tilde_theta_1, tilde_theta_2, tilde_theta_3)
      }
    }
    run_count = run_count + 1
    tol = max(sapply(vertices,f)) - min(sapply(vertices,f))
    print(c("Round", run_count))
  }
  
  print('Pontos:')
  print(vertices)
  print('Valores:')
  print(format(sapply(vertices,f), scientific=F))
}

##############
# RESULTADOS #
##############

# Essa função tem quatro mínimos locais idênticos e um máximo local. 
# Vou tentar fazer o algoritmo pegar cada um desses pontos chutando diferentes vértices iniciais.

# Minima examples 
nelder_mead(list(c(-0.2, -1), c(-0.3, -0.9), c(0,-0.5))) # Pegou (3,2)
nelder_mead(list(c(2, 3), c(7, 3), c(6,4))) # Pegou (-2.8, 3.13)
nelder_mead(list(c(1, 7), c(2, 1), c(3,5))) # Pegou (3.58, -1.84)
nelder_mead(list(c(-5, -4), c(-3, -2), c(-1,-4))) # Pegou (-3.77. -3.2)

# Maxima (NÃO CONSIGO PEGAR ESSE MÁXIMO DE JEITO NENHUM)
nelder_mead(list(c(-0.27, -0.92), c(-0.26, -0.91), c(-0.28, -0.91))) # ... (-0.27, -0.92)

# What? Esse esquema de tolerância parece problemático, faz a função pegar uns pontos nonsense.

nelder_mead(list(c(1, 3), c(1, 5), c(1,1))) # ???
nelder_mead(list(c(2, 3), c(2, 3), c(1,1))) # ???



########## Q7

# loading used packages

library(tidyverse)
library(purrr)
library(tidymodels)

prisoner = load('PS3/prisoner.Rds')

#saving eulers constant

euler = -digamma(1)


# the functions below calculate the CCPs for a given set of parameters and prison times;
# the first one calculates the formulas presented on the PDF, and the second is a wrapper
# that applies the first to the dataset, and also calculates the conditional probability
# for t = 2

calculates_ccps_given_t = function(alpha, alpha1, alpha2, t1, t2) {
  U_deny_t1 = alpha*(alpha1*10 + alpha2*100)
  U_deny_t2 = alpha*(alpha1*10 + alpha2*100)
  U_confess_t1 = alpha1*t1 + alpha2*(t1)^2
  U_confess_t2 = alpha1*t2 + alpha2*(t2)^2
  v_deny = max(U_deny_t1 + euler, U_confess_t2 + euler)
  V_confess = euler
  
  ccp_confess_t1 = exp(V_confess + U_confess_t1)/(exp(V_confess + U_confess_t1) + 
                                                    exp(v_deny + U_deny_t1))
  
  ccp_confess_t2 = exp(U_confess_t2)/(exp(U_confess_t2) + exp(U_deny_t2))
  
  return(c(ccp_confess_t1, ccp_confess_t2))
}


calculates_ccps_of_data = function(alpha, alpha1, alpha2, data) {
  ccps = calculates_ccps_given_t(alpha, alpha1, alpha2, data$offer1, data$offer2)
  
  ccps[2] = ifelse(data$choice1 == "confess", 1, ccps[2])
  
  return(ccps)
}



# given the CCP functions, we calculate the log-likelihood of a given param vector
# by applying the funcionts to every data point, and then use the standard log-lik formula

log_lik_function = function(params, data) {
  confess_probabilities = do.call('rbind', map(1:1000, 
                                               ~ calculates_ccps_of_data(params[1], params[2], params[3], 
                                                                         data[.x,])))
  
  df_with_probabilities = cbind(data, confess_probabilities)
  
  -sum(log(ifelse(df_with_probabilities$choice1 == "deny", 
                  1-df_with_probabilities$`1`, 
                  df_with_probabilities$`1`)*
             ifelse(df_with_probabilities$choice2 == "deny", 
                    1-df_with_probabilities$`2`, 
                    df_with_probabilities$`2`)))
  
}


# to find the estimate, we simply use a standard optimizer

result = optim(c(1/2, -1, -0.2), log_lik_function, data = prisoner)




############# Q8



# to calculate the bootstraps, we use the bootstrap functions from tidymodels, which 
# saves us from having to code a big wrapper to map/ from using a loop

# the function below is a small wrapper designed to work for a given resample of the data

log_lik_bootstrap_estimator = function(data) {
  optim(c(1/2, -1, -0.2), log_lik_function, data = analysis(data))
}


# then we generate the resamples

boots = bootstraps(prisoner, times = 100)


# and apply the minimization over every resample (TAKES A WHILE TO RUN)

bootstrapestimate = boots %>%
  mutate(model = map(splits, log_lik_bootstrap_estimator))


# after that, we retrieve every coefficient estimate

bootstrap_coefficients = do.call("rbind", map(1:100, ~ bootstrapestimate$model[[.x]]$par))


# and generate the varcov matrix

varcov = t(apply(bootstrap_coefficients, 2, scale, scale=FALSE, center=TRUE))%*%
  apply(bootstrap_coefficients, 2, scale, scale=FALSE, center=TRUE)/99


# the estimates and variances are checked below

mean_est = colMeans(bootstrap_coefficients)
sd_est = sqrt(diag(varcov))


# as are each requested confidence intervals

alpha_conf_int = c(sort(bootstrap_coefficients[,1])[6], sort(bootstrap_coefficients[,1])[95])
alpha1_conf_int = c(sort(bootstrap_coefficients[,2])[6], sort(bootstrap_coefficients[,2])[95])
alpha2_conf_int = c(sort(bootstrap_coefficients[,3])[6], sort(bootstrap_coefficients[,3])[95])
