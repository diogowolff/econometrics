library(Matching)

df_pacote = df %>% mutate(post = ifelse(qtr > 0, "post", "pre")) %>%
  group_by(post, id) %>%
  mutate(mean = mean(earn)) %>%
  slice(1) %>% 
  dplyr::select(d, p, post, mean) %>% 
  pivot_wider(names_from = post, values_from = mean) %>%
  mutate(y = post - pre) %>% drop_na()

teste = Match(Y = df_pacote$y, Tr = df_pacote$d, X = df_pacote$p, CommonSupport = TRUE)
summary(teste)
