library(tidyverse)
library(ggplot2)
library(RANN)
library(dotwhisker)

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


match_did_nn_estimator = function(timeframe, data) {
  post_vec = 1:timeframe
  pre_vec = -post_vec
  
  valid_indexes_post = data %>% filter(qtr %in% post_vec) %>%
    select(id) %>% unique()
  
  valid_indexes_pre = data %>% filter(qtr %in% pre_vec) %>%
    select(id) %>% unique()
  
  valid_indexes = intersect(valid_indexes_post, valid_indexes_pre) %>% unlist()
  
  df_itemd = data %>% filter(id %in% valid_indexes) %>%
    filter(qtr>= -timeframe) %>% filter(qtr <= timeframe) %>%
    mutate(post = ifelse(qtr > 0, "post", "pre")) %>%
    group_by(post, id) %>%
    mutate(mean = mean(earn)) %>% slice(1) %>%
    ungroup() %>%
    select(id, d, mean, p, post) %>%
    pivot_wider(names_from = post, values_from = mean) %>%
    mutate(diff_in_means = post-pre) %>%
    select(-c(post, pre))
  
  
  score_control = df_itemd[df_itemd$d == 0, ]   # i generate which is which to get the
  score_treatment =  df_itemd[df_itemd$d == 1, ]  # nearest neighbors using RANN
  
  nn_indexes = nn2(score_control$p, score_treatment$p, k = 2)
  df_nn = as.data.frame(nn_indexes$nn.idx)
  
  control_diff = (score_control[df_nn$V1, 'diff_in_means'] + score_control[df_nn$V2, 'diff_in_means'])/2
  
  alpha = colMeans(score_treatment$diff_in_means - control_diff)
  
  return(alpha)
}

match_did_nn_estimator(4, df)
match_did_nn_estimator(5, df)
match_did_nn_estimator(6, df)




# kernel matching

epanechnikov = function(z) {0.75*(1-z^2)*as.numeric(abs(z)<1)}

match_did_kernel_estimator = function(timeframe, data) {
  post_vec = 1:timeframe
  pre_vec = -post_vec
  
  valid_indexes_post = data %>% filter(qtr %in% post_vec) %>%
    select(id) %>% unique()
  
  valid_indexes_pre = data %>% filter(qtr %in% pre_vec) %>%
    select(id) %>% unique()
  
  valid_indexes = intersect(valid_indexes_post, valid_indexes_pre) %>% unlist()
  
  df_itemd = data %>% filter(id %in% valid_indexes) %>%
    filter(qtr>= -timeframe) %>% filter(qtr <= timeframe) %>%
    mutate(post = ifelse(qtr > 0, "post", "pre")) %>%
    group_by(post, id) %>%
    mutate(mean = mean(earn)) %>% slice(1) %>%
    ungroup() %>%
    select(id, d, mean, p, post) %>%
    pivot_wider(names_from = post, values_from = mean) %>%
    mutate(diff_in_means = post-pre) %>%
    select(-c(post, pre))
  
  
  score_control = df_itemd[df_itemd$d == 0, ]   # i generate which is which to get the
  score_treatment =  df_itemd[df_itemd$d == 1, ]  #
  
  
  combinations = expand.grid(score_treatment$p, score_control$p)
  h = 2.345*nrow(df_itemd)^(-1/5)*sd(df_itemd$p)
  
  Kern_mat = matrix(epanechnikov((combinations[, 2] - combinations[, 1])/h), nrow = nrow(score_treatment))
  W_mat = apply(Kern_mat, 2, function(i) i/sum(i))
  
  alpha = colMeans(score_treatment$diff_in_means - W_mat %*% score_control$diff_in_means)
  
  return(alpha)
}

match_did_kernel_estimator(4, df)
match_did_kernel_estimator(5, df)
match_did_kernel_estimator(6, df)





# e)

# achei que era pra fazer common support ja na d), entao acabei fazendo a e) antes da hora

df %>% ggplot(aes(x = p, group = d)) + 
  geom_density() + facet_wrap(~d)

control_p = df %>% filter(d == 0) %>% select(p)
control_t = df %>% filter(d == 1) %>% select(p)

range_control = c(min(control_p), max(control_p))
range_treatment = c(min(control_t), max(control_t))

lower_bound_score = max(min(control_p), min(control_t))
upper_bound_score = min(max(control_p), max(control_t))




# f)

common_support_df = df %>% filter(between(p, lower_bound_score, upper_bound_score))

match_did_nn_estimator(4, common_support_df)
match_did_nn_estimator(5, common_support_df)
match_did_nn_estimator(6, common_support_df)

match_did_kernel_estimator(4, common_support_df)
match_did_kernel_estimator(5, common_support_df)
match_did_kernel_estimator(6, common_support_df)







################### q2


df2 = read_csv('enoe_q219-q122_married_female.csv')

glimpse(df2)


# a)
table(df2$time)
df2_aug = df2 %>% mutate(event = time - 4,
                         Dm3 = ifelse(event == -3, 1, 0),
                         Dm2 = ifelse(event == -2, 1, 0),
                         Dm1 = ifelse(event == -1, 1, 0),
                         D0 = ifelse(event == 0, 1, 0),
                         D1 = ifelse(event == 1, 1, 0),
                         D2 = ifelse(event == 2, 1, 0),
                         D3 = ifelse(event == 3, 1, 0),
                         D4 = ifelse(event == 4, 1, 0),
                         D5 = ifelse(event == 5, 1, 0),
                         D6 = ifelse(event == 6, 1, 0),
                         D7 = ifelse(event == 7, 1, 0),
                         edusq = edu^2)

glimpse(df2_aug)


# b)

data_reg = df2_aug %>% select(-c(quarter, ent, Dm1, event,
                                 time, newid, dmarr))
           
data_unemp = data_reg %>% select(-c(formal_new, informal_new, inact))
  
reg_unemp = lm(unemp ~ ., data_unemp)

broom::tidy(reg_unemp) %>% filter(str_detect(term, 'D')) %>% dwplot()


# unemployment shoots up right after covid and then goes down


data_inac = data_reg %>% select(-c(formal_new, informal_new, unemp))

reg_inac = lm(inact ~ ., data_inac)

broom::tidy(reg_inac) %>% filter(str_detect(term, 'D')) %>% dwplot()


# inactivity shoots up after covid and takes a while to come back down


data_formal_emp = data_reg %>% select(-c(inact, informal_new, unemp))

reg_form = lm(formal_new ~., data_formal_emp)

broom::tidy(reg_form) %>% filter(str_detect(term, 'D')) %>% dwplot()


# formal employment is hit and stays down in the sample



data_inform = data_reg %>% select(-c(inact, formal_new, unemp))

reg_infor = lm(informal_new ~ ., data_inform)

broom::tidy(reg_infor) %>% filter(str_detect(term, 'D')) %>% dwplot()


# informal employment is hit but quickly recovers; it seems that the recovery
# in unemployment and activity level is driven by informal employment only



# c)

# the event variable is already a quarter-specific fixed effect, so there's no point
# adding another one.