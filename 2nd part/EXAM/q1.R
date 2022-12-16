############################################
# Econometrics - 2nd part - Final Exam     #
# Students: Luan Borelli and Diogo Wolff   #
############################################

library(tidyverse)
library(ggplot2)
library(purrr)

setwd(this.path::here())

######
# Q1 #
######

df = read_csv('progresaRDD_exam_q1.csv')

# Notice that the xnorm variable already is "xi-c".

### (a)

# It is easy to implement these regressions using the poly() function.
# We are considering here full-specified (complete) polynomials; 
# i.e., polynomials considering all the degrees.

## Quadratic specification, full polynomial 

quad_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 2, raw = TRUE), data = df) # Regressing
quad_fit = quad_reg$fitted.values # Fitted values 
quad_cons = quad_reg$coefficients[1] # Constant (for plotting purposes)
quad_alpha_rdd = quad_reg$coefficients[2] # alpha_RDD estimate.

quad_alpha_rdd

# Plotting (using xnorm in the x-axis): 
ggplot(data.frame(df$xnorm, quad_fit), aes(x = df.xnorm, y = quad_fit)) + 
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', size=1, color="red") + 
  geom_hline(yintercept = quad_cons, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  geom_hline(yintercept = quad_cons + quad_alpha_rdd, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quadratic specification with full polynomial")

# Including a confidence interval: 
ggplot(data.frame(df$xnorm, quad_fit), aes(x = df.xnorm, y = quad_fit)) + 
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', size=.1, color="black") +
  geom_segment(aes(x = 0, xend = 0, y = quad_cons + confint(quad_reg)[2,1], 
                   yend = quad_cons + confint(quad_reg)[2,2])) +
  geom_segment(aes(x=-1, xend=1, y=quad_cons + confint(quad_reg)[2,1],
                   yend = quad_cons + confint(quad_reg)[2,1])) +
  geom_segment(aes(x=-1, xend=1, y=quad_cons + confint(quad_reg)[2,2],
                   yend = quad_cons + confint(quad_reg)[2,2])) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quadratic specification with full polynomial")

# Plotting (using xbin in the x-axis): 
ggplot(data.frame(df$xbin, quad_fit), aes(x = df.xbin, y = quad_fit)) + 
  geom_point(alpha = 0.2) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (bin'd xnorm values)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quadratic specification with full polynomial (using bin'd xnom values)")

## Cubic specification, full polynomial 

cub_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 3, raw = TRUE), data = df) # Regressing
cub_fit = cub_reg$fitted.values # Fitted values
cub_cons = cub_reg$coefficients[1] # Constant (for plotting purposes)
cub_alpha_rdd = cub_reg$coefficients[2] # alpha_RDD estimate.

cub_alpha_rdd

# Plotting (using xnorm in the x-axis): 
ggplot(data.frame(df$xnorm, cub_fit), aes(x = df.xnorm, y = cub_fit)) + 
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', size=1, color="red") + 
  geom_hline(yintercept = cub_cons, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  geom_hline(yintercept = cub_cons + cub_alpha_rdd, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, cubic specification with full polynomial")

# Including a confidence interval:  
ggplot(data.frame(df$xnorm, cub_fit), aes(x = df.xnorm, y = cub_fit)) + 
  geom_point(alpha = 0.2)  +
  geom_vline(xintercept = 0, linetype = 'dashed', size=.1, color="black") +
  geom_segment(aes(x = 0, xend = 0, y = cub_cons + confint(cub_reg)[2,1], 
                   yend = cub_cons + confint(cub_reg)[2,2])) +
  geom_segment(aes(x=-1, xend=1, y=cub_cons + confint(cub_reg)[2,1],
                   yend = cub_cons + confint(cub_reg)[2,1])) +
  geom_segment(aes(x=-1, xend=1, y=cub_cons + confint(cub_reg)[2,2],
                   yend = cub_cons + confint(cub_reg)[2,2])) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, cubic specification with full polynomial")

# Plotting (using xbin in the x-axis): 
ggplot(data.frame(df$xbin, cub_fit), aes(x = df.xbin, y = cub_fit)) + 
  geom_point(alpha = 0.2) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (bin'd xnorm values)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, cubic specification with full polynomial (using bin'd xnom values)")


## Quartic specification, full polynomial 

quar_reg = lm(y ~ I(xnorm <= 0)*poly(xnorm, 4, raw = TRUE), data = df) # Regressing
quar_fit = quar_reg$fitted.values # Fitted values
quar_cons = quar_reg$coefficients[1] # Constant (for plotting purposes)
quar_alpha_rdd = quar_reg$coefficients[2] # alpha_RDD estimate.

quar_alpha_rdd

# Plotting (using xnorm in the x-axis): 
ggplot(data.frame(df$xnorm, quar_fit), aes(x = df.xnorm, y = quar_fit)) + 
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', size=1, color="red") + 
  geom_hline(yintercept = quar_cons, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  geom_hline(yintercept = quar_cons + quar_alpha_rdd, linetype = 'dashed', size=0.5, color="black", alpha=0.3) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quartic specification with full polynomial")

# Including a confidence interval:  
ggplot(data.frame(df$xnorm, quar_fit), aes(x = df.xnorm, y = quar_fit)) + 
  geom_point(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', size=.1, color="black") +
  geom_segment(aes(x = 0, xend = 0, y = quar_cons + confint(quar_reg)[2,1], 
                                            yend = quar_cons + confint(quar_reg)[2,2])) +
  geom_segment(aes(x=-1, xend=1, y=quar_cons + confint(quar_reg)[2,1],
                   yend = quar_cons + confint(quar_reg)[2,1])) +
  geom_segment(aes(x=-1, xend=1, y=quar_cons + confint(quar_reg)[2,2],
                   yend = quar_cons + confint(quar_reg)[2,2])) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (xnorm, elegible if <= 0)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quartic specification with full polynomial")

# Plotting (using xbin in the x-axis): 
ggplot(data.frame(df$xbin, quar_fit), aes(x = df.xbin, y = quar_fit)) + 
  geom_point(alpha = 0.2) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Elegibility (bin'd xnorm values)") +
  ylab("Scholl enrollment status (Y)") +
  ggtitle("RDD fitted values, quartic specification with full polynomial (using bin'd xnom values)")


# b)

# Separating data conditional on elegibility.
data_over = df %>% filter(xnorm > 0) 
data_under = df %>% filter(xnorm <= 0)

# Defining the triangular kernel function.
kernel_calculator = function(data, bandwidth) {
  (1-abs(data/bandwidth))*I(abs(data/bandwidth) <= 1)
}

# Defining the theta estimator function.
theta_estimator = function(x, y, bandwidth) {
  K_vec = kernel_calculator(x, bandwidth)
  
  lh_mat = matrix(c(sum(K_vec), sum(K_vec*x), sum(K_vec*x), 
                    sum((K_vec*(x)^2)^2)), nrow = 2)
  
  rh_vec = c(sum(K_vec*y), sum(K_vec*y*x))
  
  solve(lh_mat) %*% rh_vec
}

# Estimating alpha:
alp_rdd_15 = (theta_estimator(data_under$xnorm, data_under$y, 15) -
  theta_estimator(data_over$xnorm, data_over$y, 15))[1]

alp_rdd_30 = (theta_estimator(data_under$xnorm, data_under$y, 30) -
                theta_estimator(data_over$xnorm, data_over$y, 30))[1]

alp_rdd # h = 15
alp_rdd_30 # h = 30


# Computing alpha estimators for bandwiths (h) running from 15 to 30
alpha_robustness = map_dbl(15:30, ~ (theta_estimator(data_under$xnorm, data_under$y, .x) -
                                   theta_estimator(data_over$xnorm, data_over$y, .x))[1])

alpha_robustness

ggplot(data.frame("bandwith" = 15:30, "alpha_rdd" = alpha_robustness), aes(x = bandwith, y = alpha_rdd)) + 
  geom_point(size=3) +
  geom_line(size=1, alpha=0.3) +
  theme_bw() + scale_color_brewer(palette="Paired") + 
  xlab("Bandwidth (h)") +
  ylab("Kernel RDD estimate (alpha_RDD)") +
  ggtitle("Triangular Kernel RDD estimates as bandwidth (h) varies from 15 to 30")


