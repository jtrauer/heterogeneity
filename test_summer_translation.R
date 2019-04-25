
# setwd("C:/Users/jtrauer/Desktop/summer")
setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/summer")
source("summer_model.R")

setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/heterogeneity")
# setwd("C://Users/jtrauer/Desktop/heterogeneity/heterogeneity")
source("function_tool_kit.R")
# library(lhs)
local_repo_directory <- getwd()
setwd(local_repo_directory)

# initially setting parameter values to what were the upper limits in sensitivity.R
params <- list(
            P_epsilon = 0.128,
            P_nu = 0.077, 
            mu = 1 / 75,
            Time_L1 = 1,
            Time_L2 = 20,
            Time_I = 3,
            P_mui = 0.787,
            beta = 30,
            cdr_b = 0.8,
            treatment_success = 0.8,
            r = 0.21,
            prop_I0 = 0.58,
            prop_I1 = 0.31,
            prop_I2 = 0.11,
            P_j = 0.058,
            P_h = 0.058,
            alpha_1 = 0.22,
            alpha_0 = 0,
            P_mui0 = 0.096,
            P_mui1 = 0.096,
            P_mui2 = 0.787)

# uncertainty_param_ranges <- c(P_epsilon = list(min = 0.074, max = 0.128),
#                               P_nu = list(min = 0.018, max = 0.077))
# 
# 
# n_runs <- 2
# lhs <- maximinLHS(n_runs, length(uncertainty_param_ranges))


# parameter processing
params$epsilon <- params$P_epsilon / params$Time_L1
params$kappa <- (1 - params$P_epsilon) / params$Time_L1
params$nu <- params$P_nu / params$Time_L2
params$j <- params$P_j / params$Time_I
params$h <- params$P_h / params$Time_I
params$universal_death_rate <- params$mu
params$mui0 <- params$P_mui0 / params$Time_I
params$mui1 <- params$P_mui1 / params$Time_I
params$mui2 <- params$P_mui2 / params$Time_I
params$gamma <- (1 - params$P_mui) / params$Time_I
params$gamma0 <- 1 / params$Time_I - (params$mu + params$h + params$mui0)
params$gamma1 <- 1 / params$Time_I - (params$mu + params$j + params$mui1)
params$gamma2 <- 1 / params$Time_I - (params$mui2 + params$mu)
params$beta_reinfection <- params$beta * params$r
params$delta <- find_delta_from_cdr(params$cdr_b, params$gamma + params$mu, 1)

# model intial conditions and integration time specification
S_init = 1
I_init = 1e-6
initial_values = c(S = S_init - I_init, L1 = 0, L2 = 0, I = I_init)
initial_values_yaye_version <- c(S = S_init - I_init, L1 = 0, L2 = 0, I0 = I_init / 3, I1 = I_init / 3, I2 = I_init / 3)
times <- seq(0, 1e2)

# create summer model and define yaye model
source("model_functions.R")

# yaye version new
yaye_version <- as.data.frame(lsoda(initial_values_yaye_version, times, YayeModel, params))

# run summer version
summer_version$run_model()

# print comparison of outputs
print("yaye version new")
print(yaye_version$I0)
print("summer version")
print(summer_version$outputs$IXinfect_0)



