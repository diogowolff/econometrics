#PS3
library(tidyverse)
library(purrr)



euler = -digamma(1)

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

log_lik_estimator = function(params) {
  confess_probabilities = matrix(unlist(map(1:1000, 
                                            ~ calculates_ccps_of_data(params[1], params[2], params[3], 
                                                                      prisoner[.x,]))),
                                 ncol = 2, byrow = TRUE)
  
  df_with_probabilities = cbind(prisoner, confess_probabilities)
  
  -sum(log(ifelse(df_with_probabilities$choice1 == "deny", 
                  1-df_with_probabilities$`1`, 
                  df_with_probabilities$`1`)*
             ifelse(df_with_probabilities$choice2 == "deny", 
                    1-df_with_probabilities$`2`, 
                    df_with_probabilities$`2`)))
  
}

optim(c(1/2, -1, -0.2), log_lik_estimator)
