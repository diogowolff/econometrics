data = readr::read_csv('WAGES_PSID1976-1982_IDs.csv')

library(tidyverse)

data_demean = data %>% mutate(EXP2 = EXP^2) %>%
  group_by(ID) %>%
  mutate_all(funs(.-mean(.)))

reg = lm(LWAGE ~ -1 + EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED, data = data_demean)

summary(reg)




# b

mean_data = data %>% mutate(EXP2 = EXP^2) %>%
  group_by(ID) %>%
  summarise(across(everything(), mean))

between_reg = lm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED, data = mean_data)

sigma_est = sqrt(1/(4165-595-10)*sum(reg$residuals^2))

sigma_bet_est = sqrt(1/595*sum(between_reg$residuals^2))

sigma_u_est = max(0, sigma_bet_est^2 - sigma_est^2/7)

rho = sigma_est/sqrt(sigma_est^2 + 7*sigma_u_est)

gls_data = data %>% mutate(EXP2 = EXP^2) %>%
  group_by(ID) %>%
  mutate_all(funs(. - (1-rho)*mean(.)))

gls_est = lm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED, data = gls_data)

summary(gls_est)



# c

dif_coefs = reg$coefficients[1:9]-gls_est$coefficients[2:10]
t(dif_coefs) %*% solve(vcov(reg)[1:9, 1:9] - vcov(gls_est)[2:10, 2:10]) %*% dif_coefs


# d

data_time_mean = data %>% mutate(EXP2 = EXP^2) %>%
  group_by(YEAR) %>%
  summarise(across(everything(), mean))

data_time_mean_stack = purrr::map_dfr(seq_len(595), ~data_time_mean)

mean_data_stack = mean_data %>% slice(rep(1:n(), each = 7))

data_twoways = data_demean - data_time_mean_stack + mean_data_stack

twols = lm(LWAGE ~ -1 + EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED,
           data = data_twoways)

summary(twols)

