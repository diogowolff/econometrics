# Importing useful packages
library(dplyr)
library(plm)

setwd('C:/Users/Luan Borelli/Desktop/EPGE/Metrics/2nd part/PS1/PS1')
raw_df <- read.csv("WAGES_PSID1976-1982_IDs.csv")

# We are interested in the following specifiction
# log(w_it) = beta_1 * X_it + beta_2 Z_it + gamma ED + alpha_i + epsilon_it

# X includes variables that are considered endognous a priori: 
# i) working experience (EXP); 
# ii) experience squared (EXP2); 
# iii) a dummy for union status (UNION);
# iv) number of weeks worked in the year (WKS); and 
# v) marital status dumy (MS).

# Z includes variables that are considered exogenous a priori: 
# i) blue-collar occupation dummy (OCC);
# ii) southern region indicator (SOUTH); 
# iii) Standard Metropolitan Statistical Area indicator (SMSA);
# iv) industry dummy (IND);
# v) black race dummy variable (BLK); and 
# vi) female dummy variable (FEM).


# Creating EXP2 variable, adding it to the dataframe and selecting variables of interest:
df <- raw_df %>% mutate(EXP2 = EXP^2) %>% select(ID, EXP, EXP2, UNION, WKS, MS, OCC, SOUTH, SMSA, IND, BLK, FEM, LWAGE, ED)


#####
# a #
#####

# We want to estimate the model implementing a fixed effects approach
# by using the within estimator for beta_1, beta_2 and gamma.

# For this we need to take the individual specific means
# with respect to time for each variable of the data set.

# Notice that the panel is balanced, with 
# T = 7 (1976 - 1982) for each individual i.

# Computing individual specific means (ISMs):

df_ISM <- df %>% group_by(ID) %>% summarise(across(everything(), list(mean))) 
df_ISM_repeated_rows <- df_ISM %>% slice(rep(1:n(), each = 7))
df_within_transformed <- df - df_ISM_repeated_rows

within_estimation <- lm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED - 1, data = df_within_transformed)

# Checking results using the 'plm' package: 
df_check <- raw_df %>% mutate(EXP2 = EXP^2) %>% select(ID, YEAR, EXP, EXP2, UNION, WKS, MS, OCC, SOUTH, SMSA, IND, BLK, FEM, LWAGE, ED)
within_estimation_check <- plm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED - 1, data = df_check, model = "within")

within_estimation_check

# Notice that dummy variable coefficients cannot be identified, since they are time invariant.

#####
# b #
#####


# Directly using the 'plm' package: 

random_effects_estimation <- plm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM + ED, data = df_check, model = "random")

random_effects_estimation


# From scratch: 

## (...)

#####
# c #
#####

# Now we want to implement the Hausman test.
# The Hausman test statistic is given by
# [(bhat_FE - bhat_RE)^2]/(var_betaFE - var_betaRE).
# It follows a chi squared distribution of k 
# degrees of freedom, where k is the dimension of 
# the coefficient vectors (bhat's).

within_coefficients <- within_estimation$coefficients[-c(10,11,12)]
re_coefficients <- c(random_effects_estimation$coefficients)[-1][-c(10,11,12)]
diff_coefficients <- within_coefficients - re_coefficients

var_betaFE <- vcov(within_estimation_check)
var_betaRE <- vcov(random_effects_estimation)[-1,][,-1][-10,][-10,][-10,][,-10][,-10][,-10]
diff_var <- var_betaFE - var_betaRE

h <- t(as.matrix(diff_coefficients)) %*% solve(diff_var) %*% as.matrix(diff_coefficients)

h

#####
# d # 
#####

# Now we want to extend the fixed effects (within) estimator
# of item (a) to include time-specific fixed effects. 

# For this, we use the TWFE estimator. 

# In order to compute the TWFE estimator we need to
# "double demean" (with respecto to individuals AND time specific means)
# each variable to be included in the regression. Then, we run OLS.

# We already have the individual-specific means.
# We need to compute the time-specific means and the "full means".


# Computing time specific means (TSMs):
df_TWFE <- raw_df %>% mutate(EXP2 = EXP^2) %>% select(ID, YEAR, EXP, EXP2, UNION, WKS, MS, OCC, SOUTH, SMSA, IND, BLK, FEM, LWAGE)
df_TSM <- df_TWFE %>% group_by(YEAR) %>% summarise(across(everything(), list(mean))) 
df_TSM_repeated_rows <- do.call("rbind", replicate(595, df_TSM, simplify = FALSE)) %>% select(ID_1, EXP_1, EXP2_1, UNION_1, WKS_1, MS_1, OCC_1, SOUTH_1, SMSA_1, IND_1, BLK_1, FEM_1, LWAGE_1)

# Computing full means:
df_full_means <- data.frame(do.call("rbind", replicate(4165, t(data.frame("Means" = colMeans(df_TWFE))), simplify = FALSE))) %>% select(ID, EXP, EXP2, UNION, WKS, MS, OCC, SOUTH, SMSA, IND, BLK, FEM, LWAGE)

# Computing the transformed variables for TWFE estimation:
df_TWFE_transformed <- df - df_ISM_repeated_rows - df_TSM_repeated_rows + df_full_means

# Estimating by least squares:
TWFE_estimation <- lm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM - 1, data = df_TWFE_transformed)

TWFE_estimation

# Notice that now the EXP coefficient is also not identified.

# Checking using the 'plm' package: 

TWFE_estimation_check <- plm(LWAGE ~ EXP + EXP2 + UNION + WKS + MS + OCC + SOUTH + SMSA + IND + BLK + FEM - 1, data = df_check, model = "within", effect = "twoways")

TWFE_estimation_check


