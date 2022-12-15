library(tidyverse)
library(ggplot2)
library(purrr)

setwd(this.path::here())

df = read_csv('progresaRDD_exam_q1.csv')

# xnorm already is xi-c


# a)

# its easy to implement using poly()

quad_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 2, raw = TRUE), data = df)
quad_fit = quad_reg$fitted.values
plot(df$xnorm, quad_fit)

cub_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 3, raw = TRUE), data = df)
cub_fit = cub_reg$fitted.values
plot(df$xnorm, cub_fit)

quar_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 4, raw = TRUE), data = df)
quar_fit = quar_reg$fitted.values
plot(df$xnorm, quar_fit)


# it's working, just have to make some pretty graphs in ggplot


# b)

data_over = df %>% filter(xnorm > 0)
data_under = df %>% filter(xnorm <= 0)

kernel_calculator = function(data, bandwidth) {
  (1-abs(data/bandwidth))*I(abs(data/bandwidth) <= 1)
}


theta_estimator = function(x, y, bandwidth) {
  K_vec = kernel_calculator(x, bandwidth)
  
  lh_mat = matrix(c(sum(K_vec), sum(K_vec*x), sum(K_vec*x), 
                    sum((K_vec*(x)^2)^2)), nrow = 2)
  
  rh_vec = c(sum(K_vec*y), sum(K_vec*y*x))
  
  solve(lh_mat) %*% rh_vec
}

alp_rdd = (theta_estimator(data_over$xnorm, data_over$y, 15) -
  theta_estimator(data_under$xnorm, data_under$y, 15))[1]

alpha_robustness = map_dbl(10:35, ~ (theta_estimator(data_over$xnorm, data_over$y, .x) -
                                   theta_estimator(data_under$xnorm, data_under$y, .x))[1])

plot(10:35, alpha_robustness)

