
source("model_functions.R")
local_repo_directory <- getwd()
setwd("C:/Users/jtrauer/Desktop/summer")
source("summer_model.R")
setwd(local_repo_directory)

# initially setting parameter values to what were the upper limits
params <- list(N = 1,
            P_epsilon = 0.128,
            P_nu = 0.077, 
            alpha = 0.22,
            mu = 1 / 75,
            Time_L1 = 1.0,
            Time_L2 = 20,
            Time_I = 3,
            P_mui0 = 0.096,
            P_mui1 = 0.096,
            P_mui2 =0.787,
            r = 0.21,
            beta = 60,
            cdr = 0.8,
            treatment_success = 0.8,
            p1 = 0,
            p2 = 0,
            crude_birth_rate = 1 / 75,
            universal_death_rate = 1 / 75,
            place_holder = .1)

params$epsilon <- params$P_epsilon / params$Time_L1
params$kappa <- (1 - params$P_epsilon) / params$Time_L1
params$nu <- params$P_nu / params$Time_L2


S_init = 1
L1_init = 0
L2_init = 0
I_init = 1e-6

initial_values = c(S = S_init - I_init, L1 = 0, L2 = 0, I = I_init)
initial_model_run_duration <- 1e2
times <- seq(0, initial_model_run_duration)

# yaye version
yaye_version <- as.data.frame(lsoda(initial_values, times, Abbreviated_Model, params))
print("yaye version")
print(yaye_version$I)


# summer version
summer_version <- EpiModel$new(times, names(initial_values), as.list(initial_values), params,
                              list(c("infection_frequency", "beta", "S", "L1"),
                                   c("standard_flows", "kappa", "L1", "L2"),
                                   c("standard_flows", "epsilon", "L1", "I"),
                                   c("standard_flows", "nu", "L2", "I")),
                              infectious_compartment="I", initial_conditions_sum_to_total = FALSE, report_progress = FALSE, 
                              birth_approach = "replace_deaths", entry_compartment = "S")

summer_version$run_model()
print("summer version")
print(summer_version$outputs$I)




