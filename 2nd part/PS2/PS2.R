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



# first observation: in order to subtract individual by individual, i need to select only 
# those in the sample that have both data before and after the intervention


# IMPORTANT!!!! THIS MATCHES PEOPLE WHO MIGHT GET CUT WHEN LOOKING AT T=4 OR 5
# WE MUST CHANGE THE DF GENERATION TO BE DEPENDENT ON TIMEFRAME
# DOING THAT WE WILL NOT NEED ALL THE OTHER FILTERS I THINK


df_itemd = df %>% group_by(id) %>% summarise(max = max(qtr),
                                                        min = min(qtr)) %>%
  mutate(bool = ifelse(max*min < 0, 1, 0)) %>%   # this checks if there's a positive and a negative
  right_join(df) %>%                             # quarter in the data, for every id
  filter(bool == 1) %>%
  select(-c(bool, max, min))



# 2-NN estimator

score_control = df_itemd[df_itemd$d == 0, ]     # i generate which is which to get the
score_treatment =  df_itemd[df_itemd$d == 1, ]  # nearest neighbors using RANN

nn_indexes = nn2(score_control$p, score_treatment$p, k = 2)
df_nn = as.data.frame(nn_indexes$nn.idx)


nn_did_estimator = function(timeframe) {
  
  # given a timeframe, we (for example) take the treatment data that is after treat,
  # subset it to data that is in the timeframe, then calculate the mean for did
  
  yt_t1 = score_treatment[score_treatment$qtr > 0,] %>% filter(abs(qtr) <= timeframe) %>%
    group_by(id) %>% summarise(mean = mean(earn))
  yt_t0 = score_treatment[score_treatment$qtr < 0,] %>% filter(abs(qtr) <= timeframe) %>%
    group_by(id) %>% summarise(mean = mean(earn))
  
  yc_t1 = score_control[score_control$qtr > 0,] %>% filter(abs(qtr) <= timeframe) %>%
    group_by(id) %>% summarise(mean = mean(earn))
  yc_t0 = score_control[score_control$qtr < 0,] %>% filter(abs(qtr) <= timeframe) %>%
    group_by(id) %>% summarise(mean = mean(earn))
  
  # we then join all this data in the original dataframe to get the dataset with 
  # differences for every individual
  # this dataset contains the mean earn for ind. at t1 and t0, whether he's C or T
  y = df_itemd %>% left_join(yt_t1, by = 'id') %>% 
    left_join(yt_t0, by = 'id') %>% left_join(yc_t1, by = 'id') %>% 
    left_join(yc_t0, by = 'id') %>% filter(rowMeans(is.na(.)) < 0.5) %>%
    rename('yt_t1' = mean.x, 'yt_t0' = mean.y, 'yc_t1' = mean.x.x, 'yc_t0' = mean.y.y) %>%
    mutate(yt_diff = yt_t1 - yt_t0,
           yc_diff = yc_t1 - yc_t0)
  
  
  # given the NN we found, we rearrange the control dataset accordingly to the match
  # in treatment data, then join with all the data and get the yc_diff in the correct order
  
  nn1 = score_control[df_nn$V1,] %>% left_join(y) %>% select(yc_diff)
  nn2 = score_control[df_nn$V2,] %>% left_join(y) %>% select(yc_diff)

  
  # we now have the vector of yt diffs and the vector of yc diffs all in the right order
  # so just apply the formula with weight 1/2 for each neighbour and 0 otherwise (just a mean really)
  
  alpha = sum(y[y$d == 1, 'yt_diff'] - (nn1+nn2)/2, na.rm = T)/nrow(y)
  
  return(alpha)
}

nn_did_estimator(6)



# e)

# achei que era pra fazer common support ja na d), entao acabei fazendo a e) antes da hora

df %>% ggplot(aes(x = p, group = d)) + 
  geom_density() + facet_wrap(~d)

range_control = c(min(score_control), max(score_control))

range_treatment = c(min(score_treatment), max(score_treatment))

lower_bound_score = max(min(score_control), min(score_treatment))
upper_bound_score = min(max(score_control), max(score_treatment))

