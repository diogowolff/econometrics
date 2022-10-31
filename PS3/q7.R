#PS3


########## q7

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
  V_confess = U_confess_t2 + euler
  
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




############# q8



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
