library(tidyverse)
library(ggplot2)
library(RANN)

setwd(this.path::here())


#################### q1


df = read_csv("jtpa_ps2.csv")


# a)

did_estimator = function(timeframe) {
  y_of_treated_after_treatment = colMeans(df[df$d == 1 & df$qtr > 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_treated_before_treatment = colMeans(df[df$d == 1 & df$qtr < 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_control_after_treatment = colMeans(df[df$d == 0 & df$qtr > 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_control_before_treatment = colMeans(df[df$d == 0 & df$qtr < 0 & abs(df$qtr) <= timeframe, 'earn'])
  
  alpha = (y_of_treated_after_treatment - y_of_treated_before_treatment) -
    (y_of_control_after_treatment - y_of_control_before_treatment)

  return(alpha)  
}

did_estimator(4)
did_estimator(5)
did_estimator(6)

# it seems that the longer we collect data on, the smaller the effect becomes.


# b)

df_itemb = df %>%
  mutate(post = ifelse(qtr>0, 1, 0))

reg_itemb_calculator = function(timeframe) {
  lm(earn ~ d + post + post*d, data = df_itemb[abs(df_itemb$qtr) <= timeframe, ])
}

summary(reg_itemb_calculator(4))
summary(reg_itemb_calculator(5))
summary(reg_itemb_calculator(6))

# the intercept of the full regression is precisely the average income of the control group before treatment, while
# the intercept plus d is the average of the treated group before treatment. intercept plus post is the 
# average of the control group after treatment, and intercept + d + post + d:post is the treated after treatment.


# c)

df_itemc = df %>% 
  group_by(d, qtr) %>%
  summarise(mean = mean(earn)) %>%
  mutate(post = ifelse(qtr>0, 1, 0))
  
ggplot(df_itemc, aes(x = qtr, y = mean, color = as.factor(d) )) + geom_line(linetype = 'dashed') +
  geom_line(data = df_itemc[df_itemc$qtr <= -1, ], aes(x = qtr, y = mean, color =
                                                                             as.factor(d)), size = 1.1) +
  geom_line(data = df_itemc[df_itemc$qtr >= 1, ], aes(x = qtr, y = mean, color =
                                                         as.factor(d)), size = 1.1) +
  labs(color = 'Treatment Status')


# as can be seen, the trends seem to be somewhat parallel after the treatment, but before treatment, the two groups
# have clearly different behaviours over time. This is clearly seen as the control group's average earnings are
# increasing throughout the sample, while the treatment group's average earnings are decreasing before treatment.


# d)

# 2-NN estimator

score_control = df[df$d == 0, ]$p
score_treatment =  df[df$d == 1, ]$p

nn_indexes = nn2(score_control, score_treatment, k = 2)

nn_did_estimator = function(timeframe) {
  #mapeia sobre os tratados: pro tratado i, pega a media timeframe p frente - timeframe pra tras
  # - soma de [i, 1]*media timeframe do controle p frente - ... p tras com [i, 2] do mesmo
}


# e)

# achei que era pra fazer common support ja na d), entao acabei fazendo a e) antes da hora

df %>% ggplot(aes(x = p, group = d)) + 
  geom_density() + facet_wrap(~d)

range_control = c(min(score_control), max(score_control))

range_treatment = c(min(score_treatment), max(score_treatment))

lower_bound_score = max(min(score_control), min(score_treatment))
upper_bound_score = min(max(score_control), max(score_treatment))

