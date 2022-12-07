############################################
# Econometrics - 2nd part - Problem Set 2  #
# Students: Luan Borelli and Diogo Wolff   #
############################################

## (!) Check LaTeX PDF document for further details on these questions.

# Importing useful packages: 
library(tidyverse)
library(ggplot2)
library(RANN)
library(dotwhisker)
library(stargazer)

######
# Q1 #
######

# Importing data:
setwd(this.path::here())
df = read_csv("jtpa_ps2.csv")

### (a)

# Outcome means for each case: 

y_of_treated_after_treatment_4 = colMeans(df[df$d == 1 & df$qtr > 0 & abs(df$qtr) <= 4, 'earn'])
y_of_treated_before_treatment_4 = colMeans(df[df$d == 1 & df$qtr < 0 & abs(df$qtr) <= 4, 'earn'])
y_of_control_after_treatment_4 = colMeans(df[df$d == 0 & df$qtr > 0 & abs(df$qtr) <= 4, 'earn'])
y_of_control_before_treatment_4 = colMeans(df[df$d == 0 & df$qtr < 0 & abs(df$qtr) <= 4, 'earn'])

y_of_treated_after_treatment_5 = colMeans(df[df$d == 1 & df$qtr > 0 & abs(df$qtr) <= 5, 'earn'])
y_of_treated_before_treatment_5 = colMeans(df[df$d == 1 & df$qtr < 0 & abs(df$qtr) <= 5, 'earn'])
y_of_control_after_treatment_5 = colMeans(df[df$d == 0 & df$qtr > 0 & abs(df$qtr) <= 5, 'earn'])
y_of_control_before_treatment_5 = colMeans(df[df$d == 0 & df$qtr < 0 & abs(df$qtr) <= 5, 'earn'])

y_of_treated_after_treatment_6 = colMeans(df[df$d == 1 & df$qtr > 0 & abs(df$qtr) <= 6, 'earn'])
y_of_treated_before_treatment_6 = colMeans(df[df$d == 1 & df$qtr < 0 & abs(df$qtr) <= 6, 'earn'])
y_of_control_after_treatment_6 = colMeans(df[df$d == 0 & df$qtr > 0 & abs(df$qtr) <= 6, 'earn'])
y_of_control_before_treatment_6 = colMeans(df[df$d == 0 & df$qtr < 0 & abs(df$qtr) <= 6, 'earn'])


# Defining a function for the difference-in-differences estimator:

did_estimator = function(timeframe) {
  y_of_treated_after_treatment = colMeans(df[df$d == 1 & df$qtr > 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_treated_before_treatment = colMeans(df[df$d == 1 & df$qtr < 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_control_after_treatment = colMeans(df[df$d == 0 & df$qtr > 0 & abs(df$qtr) <= timeframe, 'earn'])
  y_of_control_before_treatment = colMeans(df[df$d == 0 & df$qtr < 0 & abs(df$qtr) <= timeframe, 'earn'])
  
  alpha = (y_of_treated_after_treatment - y_of_treated_before_treatment) -
    (y_of_control_after_treatment - y_of_control_before_treatment)

  return(alpha)  
}

did_4 <- did_estimator(4)
did_5 <- did_estimator(5)
did_6 <- did_estimator(6)

did_4
did_5
did_6

# It seems that the longer we collect data on, the smaller the effect becomes.


### (b) 


df_itemb = df %>%
  mutate(post = ifelse(qtr>0, 1, 0))

reg_itemb_calculator = function(timeframe) {
  lm(earn ~ d + post + post*d, data = df_itemb[abs(df_itemb$qtr) <= timeframe, ])
}

summary(reg_itemb_calculator(4))
summary(reg_itemb_calculator(5))
summary(reg_itemb_calculator(6))

# Generating LaTeX table with all the results.
stargazer(reg_itemb_calculator(4), reg_itemb_calculator(5), reg_itemb_calculator(6), no.space=TRUE, digits=2)

# The intercept of the full regression is precisely the average income of the control group before treatment, while
# the intercept plus d is the average of the treated group before treatment. Intercept plus post is the 
# average of the control group after treatment, and intercept + d + post + d:post is the treated after treatment.

### (c)

df_itemc = df %>% 
  group_by(d, qtr) %>%
  summarise(mean = mean(earn)) %>%
  mutate(post = ifelse(qtr>0, 1, 0))
  
ggplot(df_itemc, aes(x = qtr, y = mean, color = as.factor(d) )) + geom_line(linetype = 'dashed') +
  geom_line(data = df_itemc[df_itemc$qtr <= -1, ], aes(x = qtr, y = mean, color =
                                                                             as.factor(d)), size = 1.1) +
  geom_line(data = df_itemc[df_itemc$qtr >= 1, ], aes(x = qtr, y = mean, color =
                                                         as.factor(d)), size = 1.1) +
  labs(color = 'Treatment Status') +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha=0.5) + 
  annotate(geom = "text", x = -0.3, y = 750, label = "TREATMENT ASSIGNMENT", color = "black", alpha=0.5,
           angle = 90) + 
  theme_bw() + scale_color_brewer(palette="Paired") + 
  ggtitle("Mean earnings per quarter for treated and non-treated individuals")


# As can be seen, the trends seem to be somewhat parallel after the treatment, but before treatment, the two groups
# have clearly different behaviours over time. This is clearly seen as the control group's average earnings are
# increasing throughout the sample, while the treatment group's average earnings are decreasing before treatment.

### (d)

# 2-NN estimator


match_did_nn_estimator = function(timeframe, data) {
  post_vec = 1:timeframe        # Generate the array of necessary observations.
  pre_vec = -post_vec
  
  valid_indexes_post = data %>% filter(qtr %in% post_vec) %>%  # Check if this ID has a valid post 
    select(id) %>% unique()                                    # observation.
  
  valid_indexes_pre = data %>% filter(qtr %in% pre_vec) %>%   # The same, but for before treatment.
    select(id) %>% unique()
  
  valid_indexes = intersect(valid_indexes_post, valid_indexes_pre) %>% unlist() # Get IDs with both.
  
  df_itemd = data %>% filter(id %in% valid_indexes) %>%        # Get only correct IDs.
    filter(qtr>= -timeframe) %>% filter(qtr <= timeframe) %>%  # Get only relevant obs from them.
    mutate(post = ifelse(qtr > 0, "post", "pre")) %>%          # Split data into pre and post treat.
    group_by(post, id) %>%
    mutate(mean = mean(earn)) %>% slice(1) %>%     # Calculate average per id per pre/post treat,
    ungroup() %>%                                    # then convert the df into a table with only
    select(id, d, mean, p, post) %>%                 # id, treatment status, and avg before/after treat.
    pivot_wider(names_from = post, values_from = mean) %>%
    mutate(diff_in_means = post-pre) %>%            # Calculate difference in means for each ind.
    select(-c(post, pre))                           # Drop unnecessary variables.
  
  score_control = df_itemd[df_itemd$d == 0, ]   # i generate which is which to get the
  score_treatment =  df_itemd[df_itemd$d == 1, ]  # nearest neighbors using RANN
  
  nn_indexes = nn2(score_control$p, score_treatment$p, k = 2)
  df_nn = as.data.frame(nn_indexes$nn.idx)
  
  control_diff = (score_control[df_nn$V1, 'diff_in_means'] + score_control[df_nn$V2, 'diff_in_means'])/2
  
  alpha = colMeans(score_treatment$diff_in_means - control_diff)
  
  return(alpha)
}

did_nn_4 <- match_did_nn_estimator(4, df)
did_nn_5 <- match_did_nn_estimator(5, df)
did_nn_6 <- match_did_nn_estimator(6, df)

did_nn_4 
did_nn_5 
did_nn_6




match_did_nn_graph = function(timeframe, data) {
  post_vec = 1:timeframe        # Generate the array of necessary observations.
  pre_vec = -post_vec
  
  valid_indexes_post = data %>% filter(qtr %in% post_vec) %>%  # Check if this ID has a valid post 
    select(id) %>% unique()                                    # observation.
  
  valid_indexes_pre = data %>% filter(qtr %in% pre_vec) %>%   # The same, but for before treatment.
    select(id) %>% unique()
  
  valid_indexes = intersect(valid_indexes_post, valid_indexes_pre) %>% unlist() # Get IDs with both.
  
  df_itemd = data %>% filter(id %in% valid_indexes) %>%        # Get only correct IDs.
    filter(qtr>= -timeframe) %>% filter(qtr <= timeframe) %>%  # Get only relevant obs from them.
    mutate(post = ifelse(qtr > 0, "post", "pre"))
  
  score_control = df_itemd[df_itemd$d == 0, ] # i generate which is which to get the
  score_treatment =  df_itemd[df_itemd$d == 1, ]  # nearest neighbors using RANN
  
  nn_indexes = nn2(score_control$p, score_treatment$p, k = 2)
  df_nn = as.data.frame(nn_indexes$nn.idx)
  
  df_with_match = tibble(score_treatment, id1 = score_control[df_nn$V1, ]$id, 
                         id2 = score_control[df_nn$V2, ]$id) %>%
    left_join(select(score_control, c('id', 'earn', 'qtr')), by = c("id1" = "id", "qtr")) %>%
    left_join(select(score_control, c('id', 'earn', 'qtr')), by = c("id2" = "id", "qtr")) %>%
    mutate(earn_control = earn.y + earn)
  
  match_control = df_with_match %>% select(id, qtr, earn_control) %>%
    mutate(d=0) %>% rename('earn' = earn_control)
    
  graph_df = score_treatment %>% select(id, qtr, d, earn) %>% rbind(match_control) %>%
    group_by(qtr, d) %>% summarise(mean = mean(earn, na.rm = T)) %>%
    mutate(post = ifelse(qtr>0, 1, 0))
  
  
  ggplot(graph_df, aes(x = qtr, y = mean, color = as.factor(d) )) + geom_line(linetype = 'dashed') +
    geom_line(data = graph_df[graph_df$qtr <= -1, ], aes(x = qtr, y = mean, color =
                                                           as.factor(d)), size = 1.1) +
    geom_line(data = graph_df[graph_df$qtr >= 1, ], aes(x = qtr, y = mean, color =
                                                          as.factor(d)), size = 1.1) +
    labs(color = 'Treatment Status') +
    geom_vline(xintercept = 0, linetype = 'dashed', alpha=0.5) + 
    annotate(geom = "text", x = -0.3, y = 750, label = "TREATMENT ASSIGNMENT", color = "black", alpha=0.5,
             angle = 90) + 
    theme_bw() + scale_color_brewer(palette="Paired") + 
    ggtitle("Mean earnings per quarter for treated and non-treated individuals")
  
}



match_did_nn_graph(6, df)











# Kernel matching

epanechnikov = function(z) {0.75*(1-z^2)*as.numeric(abs(z)<1)}
timeframe = 4
data = df
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
  h = .9*nrow(score_control)^(-1/5)*sd(score_control$p)
  
  Kern_mat = matrix(epanechnikov((combinations[, 2] - combinations[, 1])/h), nrow = nrow(score_treatment))
  W_mat = t(apply(Kern_mat, 1, function(i) i/sum(i)))
  
  index_of_treat_with_near_controls = !is.na(W_mat[, 1])
  
  alpha = mean(score_treatment$diff_in_means[index_of_treat_with_near_controls] - 
                 (W_mat %*% score_control$diff_in_means)[index_of_treat_with_near_controls])
  
  return(alpha)
}

did_kernel_4 <- match_did_kernel_estimator(4, df)
did_kernel_5 <- match_did_kernel_estimator(5, df)
did_kernel_6 <- match_did_kernel_estimator(6, df)

did_kernel_4
did_kernel_5
did_kernel_6

### (e) 

df %>% ggplot(aes(x = p, group = d)) + 
  geom_density() + facet_wrap(~d) + 
  theme_bw() + scale_color_brewer(palette="Paired") + 
  ggtitle("Distribution of p for treated (d = 1) and untreated (d = 0) individuals")

control_p = df %>% filter(d == 0) %>% select(p)
control_t = df %>% filter(d == 1) %>% select(p)

range_control = c(min(control_p), max(control_p))
range_treatment = c(min(control_t), max(control_t))

range_control
range_treatment

lower_bound_score = max(min(control_p), min(control_t))
upper_bound_score = min(max(control_p), max(control_t))

lower_bound_score
upper_bound_score

### (f)

common_support_df = df %>% filter(between(p, lower_bound_score, upper_bound_score))

did_nn_cs_4 <- match_did_nn_estimator(4, common_support_df)
did_nn_cs_5 <- match_did_nn_estimator(5, common_support_df)
did_nn_cs_6 <- match_did_nn_estimator(6, common_support_df)

did_nn_cs_4
did_nn_cs_5
did_nn_cs_6

did_kernel_cs_4 <- match_did_kernel_estimator(4, common_support_df)
did_kernel_cs_5 <- match_did_kernel_estimator(5, common_support_df)
did_kernel_cs_6 <- match_did_kernel_estimator(6, common_support_df)

did_kernel_cs_4
did_kernel_cs_5
did_kernel_cs_6

######
# Q2 #
######

df2 = read_csv('enoe_q219-q122_married_female.csv')
glimpse(df2)


### (a)

table(df2$time)
df2_aug = df2 %>% mutate(event = time - 4, # Creating the event-time variable.
                         Dm3 = ifelse(event == -3, 1, 0), # Generating the set of dummies
                         Dm2 = ifelse(event == -2, 1, 0), # specific to each event time,
                         Dm1 = ifelse(event == -1, 1, 0), # from -3 to +7.
                         D0 = ifelse(event == 0, 1, 0),
                         D1 = ifelse(event == 1, 1, 0),
                         D2 = ifelse(event == 2, 1, 0),
                         D3 = ifelse(event == 3, 1, 0),
                         D4 = ifelse(event == 4, 1, 0),
                         D5 = ifelse(event == 5, 1, 0),
                         D6 = ifelse(event == 6, 1, 0),
                         D7 = ifelse(event == 7, 1, 0),
                         edusq = edu^2) # Creating the edusq (education squared) variable.

glimpse(df2_aug)

# Computing the mean of variable "event" by "time", as requested: 
check <- df2_aug %>% group_by(time) %>% summarize(event_mean = mean(event, na.rm=TRUE))
check

### (b)

data_reg = df2_aug %>% select(-c(quarter, ent, Dm1, event,
                                 time, newid, dmarr))
           
data_unemp = data_reg %>% select(-c(formal_new, informal_new, inact))
  
reg_unemp = lm(unemp ~ ., data_unemp)

broom::tidy(reg_unemp) %>% filter(str_detect(term, 'D')) %>% dwplot(dot_args= list(color = "black", size=2), whisker_args = list(color="black", size=1), size=2) + 
            theme_bw() + 
            geom_hline(yintercept='D0', linetype='dashed', alpha=0.5) +
            geom_vline(xintercept=reg_unemp$coefficients[38],linetype='dashed', alpha=0.5) + 
            geom_vline(xintercept=0,linetype='dashed', color="red", alpha=0.5) +
            ggtitle("Event-study (dot-and-whisker) plot: Unemployment") +
            xlab('Unemployment') + 
            ylab('Event-time dummy (from -3 to +7)')


# unemployment shoots up right after covid and then goes down


data_inac = data_reg %>% select(-c(formal_new, informal_new, unemp))

reg_inac = lm(inact ~ ., data_inac)

broom::tidy(reg_inac) %>% filter(str_detect(term, 'D')) %>% dwplot(dot_args= list(color = "black", size=2), whisker_args = list(color="black", size=1), size=2) + 
  theme_bw() + 
  geom_hline(yintercept='D0', linetype='dashed', alpha=0.5) +
  geom_vline(xintercept=reg_inac$coefficients[38],linetype='dashed', alpha=0.5) +
  geom_vline(xintercept=0,linetype='dashed', color="red", alpha=0.5) +
  ggtitle("Event-study (dot-and-whisker) plot: Inactivity") +
  xlab('Inactivity') + 
  ylab('Event-time dummy (from -3 to +7)')


# inactivity shoots up after covid and takes a while to come back down


data_formal_emp = data_reg %>% select(-c(inact, informal_new, unemp))

reg_form = lm(formal_new ~., data_formal_emp)

broom::tidy(reg_form) %>% filter(str_detect(term, 'D')) %>% dwplot(dot_args= list(color = "black", size=2), whisker_args = list(color="black", size=1), size=2) + 
  theme_bw() + 
  geom_hline(yintercept='D0', linetype='dashed', alpha=0.5) + 
  geom_vline(xintercept=reg_form$coefficients[38],linetype='dashed', alpha=0.5) +
  geom_vline(xintercept=0,linetype='dashed', color="red", alpha=0.5) +
  ggtitle("Event-study (dot-and-whisker) plot: Formal employment") +
  xlab('Formal employment') + 
  ylab('Event-time dummy (from -3 to +7)')


# formal employment is hit and stays down in the sample

data_inform = data_reg %>% select(-c(inact, formal_new, unemp))

reg_infor = lm(informal_new ~ ., data_inform)

broom::tidy(reg_infor) %>% filter(str_detect(term, 'D')) %>% dwplot(dot_args= list(color = "black", size=2), whisker_args = list(color="black", size=1), size=2) + 
  theme_bw() + 
  geom_hline(yintercept='D0', linetype='dashed', alpha=0.5) + 
  geom_vline(xintercept=reg_infor$coefficients[38],linetype='dashed', alpha=0.5) +
  geom_vline(xintercept=0,linetype='dashed', color="red", alpha=0.5) +
  ggtitle("Event-study (dot-and-whisker) plot: Informal employment") +
  xlab('Informal employment') + 
  ylab('Event-time dummy (from -3 to +7)')

# informal employment is hit but quickly recovers; it seems that the recovery
# in unemployment and activity level is driven by informal employment only


### (c)

# the event variable is already a quarter-specific fixed effect, so there's no point
# adding another one. this is specific to this example.
