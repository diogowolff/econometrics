###########################################
# Econometrics - Problem Set 2            #
# Students: Luan Borelli and Diogo Wolff  #
###########################################

# Importing useful packages. 

library(purrr)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(knitr)
library(kableExtra)

# Importing data. 

df <- readRDS('PS2/adpw_2017_marathon.Rds')


######
# Q6 #
######

f_hist = function(data, x0, h) { # Defining the function f_hist for the histogram.
  n = length(data)
  ind = (data > x0 - h) & (data < x0 + h)
  f_histogram = (1/n)*(sum(ind)/(2*h))
  return(f_histogram)
}

x0_grid = seq(105,355,10) # Defining the grid.
f_hist_values = map_dbl(x0_grid, f_hist, data = df$chiptime, h = 2.5) # Calculating f_hist at each point of the grid.


# Now we finally plot the histogram: 

hist_plot <- ggplot(data = data.frame(x0_grid, f_hist_values), aes(x=x0_grid, y=f_hist_values)) +
  geom_col(size = 1) +
  labs(x = "Finishing time (minutes)", y  = "Percentile") +
  geom_vline(xintercept = 4*60, size=1) + 
  geom_vline(xintercept = 4.5*60, size=1) +
  geom_vline(xintercept = 5*60, size=1) +
  annotate(geom="text", x=244, y=0.008, label="4H", angle=90) +
  annotate(geom="text", x=4.5*60+4, y=0.008, label="4.5H", angle=90) +
  annotate(geom="text", x=5*60+4, y=0.008, label="5H", angle=90)+
  ggtitle('Histogram of marathon finishing times') + 
  theme_bw()

hist_plot

######
# Q7 #
######

# The Quartic (or biweight) Kernel function is given by
# K(z) = (15/16)*(1-z^2)^2 * 1(|z| < 1)
# where 1(.) is the indicator function.

K_quartic <- function(z) { # Defining the Quartic Kernel Function.
  ind = abs(z) < 1
  value = ((15/16)*(1-z^2)^2)*ind
  return(value)
}

f_quartic = function(data_q, x0, h_q) { # Defining the analog "histogram function" using the Quartic Kernel.
  n = length(data_q)
  K_quartic_values =  purrr::map_dbl((data_q - x0)/h_q, K_quartic)
  f_quartic_value = (1/(n*h_q))*(sum(K_quartic_values))
  return(f_quartic_value)
}

# The "heuristic" optimal bandwidth that depends only on the 
# sample size is given by h* = (1/n)^0.2.

h_star = (1/length(df$chiptime))^0.2 # Defining the "heuristic" optimal bandwidth.
x0_grid_q7 = seq(105, 355, 5) # Defining the new grid.

# Obtaining the f_quartic values for each point in the grid.
f_quartic_values =  purrr::map_dbl(x0_grid_q7, f_quartic, data = df$chiptime, h = h_star)

# Plotting
quartic_plot <- ggplot() +
  geom_line(data = data.frame(x0_grid_q7, f_quartic_values), aes(x=x0_grid_q7, y=f_quartic_values), size=1) + 
  labs(x = "Finishing time (minutes)", y  = "Percentile") +
  geom_vline(xintercept = 4*60, size=1) + 
  geom_vline(xintercept = 4.5*60, size=1) +
  geom_vline(xintercept = 5*60, size=1) +
  annotate(geom="text", x=244, y=0.008, label="4H", angle=90) +
  annotate(geom="text", x=4.5*60+4, y=0.008, label="4.5H", angle=90) +
  annotate(geom="text", x=5*60+4, y=0.008, label="5H", angle=90)+
  ggtitle('Density estimates for marathon finishing times (Quartic + Heuristic Bandwith)') + 
  theme_bw()

quartic_plot


######
# Q8 #
######

# As calculated in Q3: h^* = 2.778*n^{-0.2} * min(s, IQR/1.349).
# Defining the optimal bandwith using the Quartic Kernel:
h_optimal = 2.778*(length(df$chiptime))^(-1/5)*min(sd(df$chiptime), IQR(df$chiptime)/1.349)
# Calculating the f_quartic optimal values for each point in the grid:
f_quartic_opt_values =  purrr::map_dbl(x0_grid_q7, f_quartic, data = df$chiptime, h = h_optimal) 

# Plotting the density estimates:
quartic_opt_plot <- ggplot() +
  geom_line(data = data.frame(x0_grid_q7, f_quartic_opt_values), aes(x=x0_grid_q7, y=f_quartic_opt_values), size=1) + 
  geom_line(size = 1) +
  labs(x = "Finishing time (minutes)", y  = "Percentile") +
  geom_vline(xintercept = 4*60, size=1) + 
  geom_vline(xintercept = 4.5*60, size=1) +
  geom_vline(xintercept = 5*60, size=1) +
  annotate(geom="text", x=244, y=0.008, label="4H", angle=90) +
  annotate(geom="text", x=4.5*60+4, y=0.008, label="4.5H", angle=90) +
  annotate(geom="text", x=5*60+4, y=0.008, label="5H", angle=90)+
  ggtitle('Density estimates for marathon finishing times (Quartic + Optimal Bandwith)') + 
  theme_bw()

quartic_opt_plot 

######
# Q9 #
######

# Cleaning the dataset, as requested:
clean_df <- dplyr::filter(df, !is.na(age))

# Defining the local-constant estimator function:
local_constant <- function(y, x, x0, h_lc) {
  n = length(y)
  K_quartic_values =  purrr::map_dbl((x - x0)/h_lc, K_quartic)
  sum_elements = K_quartic_values * y 
  m = sum(sum_elements)/sum(K_quartic_values)
  return(m)
}

# Defining a grid for the ages. 
# We decided to define  the same grid that will be used in Question 10.
x0_grid_q9 = seq(20, 70, 5) 

# Obtaining the values of the local-constant estimator for each point in the grid:
local_constant_values <-  purrr::map_dbl(x0_grid_q9, local_constant, y = clean_df$chiptime, x = clean_df$age, h_lc = 3)

# Plotting:
lc_plot <- ggplot() +
  geom_line(data = data.frame(x0_grid_q9, local_constant_values), aes(x=x0_grid_q9, y=local_constant_values), size=1) + 
  geom_line(size = 1) +
  labs(x = "Age", y  = "Finishing time (minutes)") +
  ggtitle("Local-constant estimator") + 
  theme_bw()

lc_plot

#######
# Q10 #
#######

### BE CAREFUL (or patient): on our computer this part of the code took some time to run.
# Outline of the solution:
# Calculate each matrix of the left and right sides, then sum over all of them 
# and solve the matrix multiplication.

x0_grid_q10 <- seq(20, 70, 5) # Defining the grid.


left_side_for_xi = function(xi, x0) {
  K = K_quartic((xi-x0)/3)
  
  return(matrix(c(K, K*(xi-x0), K*(xi-x0), K*(xi-x0)^2), nrow=2))
}

left_side_over_all_xi = function(x0) {
  Reduce('+', map(clean_df$age, ~ left_side_for_xi(.x, x0 = x0)))
}

right_side_for_xi = function(xi, yi, x0) {
  K = K_quartic((xi-x0)/3)
  
  return(matrix(c(K*yi, K*yi*(xi-x0))))
}

right_side_over_all_xi = function(x0) {
  Reduce('+', map2(clean_df$age, clean_df$chiptime, ~ right_side_for_xi(.x, .y, x0 = x0)))
}

somewhat_fast_beta_estimator = function(x0) {
  solve(left_side_over_all_xi(x0)) %*% right_side_over_all_xi(x0)
}


betas = lapply(x0_grid_q10, somewhat_fast_beta_estimator)

print("Coefficient vector results:")
print(betas)


# Plotting:
q10_graph_tibble = tibble(as.data.frame(matrix(unlist(betas), ncol = 2, byrow = TRUE)), age = x0_grid_q10) %>%
  mutate(lowery = V1-3*V2,
         lowerx = age - 3,
         uppery = V1+3*V2,
         upperx = age + 3)

ll_plot <- q10_graph_tibble %>% ggplot(aes(x=age, y=V1)) + 
  geom_line(alpha = 0.5, size = 1, linetype = 2) +
  geom_segment(aes(x=lowerx, xend=upperx, y=lowery, yend = uppery), size = 1.1) + 
  ggtitle('Local-linear estimator') + 
  labs(x='Age', y='Finishing time (minutes)') +
  theme_bw()

ll_plot

# Comparison plot of local-constant and local-linear estimators:

lc_ll_plot <- tibble(q10_graph_tibble, local_constant_values) %>%
  select(age, local_constant_values, V1) %>%
  pivot_longer(-age, names_to = 'values') %>%
  ggplot() + 
  geom_line(aes(x=age, y=value, linetype = values, color = values)) +
  scale_linetype_discrete(labels = c('LC Estimator', 'LL Estimator')) +
  scale_colour_discrete(labels = c('LC Estimator', 'LL Estimator')) +
  labs(linetype = 'Estimators evaluated at grid points',
       color = 'Estimators evaluated at grid points') +
  geom_segment(data = q10_graph_tibble,
               aes(x=lowerx, xend = upperx, y=lowery, yend=uppery), color="cyan") +
  xlab('Age') + ylab('Finishing time (minutes)') +
  theme_bw()

lc_ll_plot

# Generating the table of coefficients in LaTeX:

as.data.frame(matrix(unlist(betas), nrow=11, byrow = TRUE)) %>% 
  kbl(caption = 'Marginal effect of age on finishing time',
      format = 'latex')