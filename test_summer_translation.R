
# setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/heterogeneity")
setwd("C://Users/jtrauer/Desktop/heterogeneity/heterogeneity")
source("model_functions.R")
source("function_tool_kit.R")
local_repo_directory <- getwd()
setwd("C:/Users/jtrauer/Desktop/summer")
# setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/summer")
source("summer_model.R")
setwd(local_repo_directory)

# initially setting parameter values to what were the upper limits
params <- list(
            P_epsilon = 0.128,
            P_nu = 0.077, 
            mu = 1 / 75,
            Time_L1 = 1.0,
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
            alpha_0 = 0)

# parameter processing
params$epsilon <- params$P_epsilon / params$Time_L1

# params$epsilon <- 0

params$kappa <- (1 - params$P_epsilon) / params$Time_L1
params$nu <- params$P_nu / params$Time_L2

params$j <- params$P_j / params$Time_I
params$h <- params$P_h / params$Time_I

# params$nu <- 0

params$universal_death_rate <- params$mu
params$mui <- params$P_mui / params$Time_I
params$gamma <- (1 - params$P_mui) / params$Time_I
params$beta_reinfection <- params$beta * params$r

params$delta <- find_delta_from_cdr(params$cdr_b, params$gamma + params$mu, 1)

# model intial conditions and integration time specification
S_init = 1
L1_init = 0
L2_init = 0
I_init = 1e-6
initial_values = c(S = S_init - I_init, L1 = 0, L2 = 0, I = I_init)
initial_values_yaye_version <- c(S = S_init - I_init, L1 = 0, L2 = 0, I0 = I_init / 3, I1 = I_init / 3, I2 = I_init / 3)
initial_model_run_duration <- 1e2
times <- seq(0, initial_model_run_duration)

# yaye version old
# yaye_version <- as.data.frame(lsoda(initial_values, times, Abbreviated_Model_old, params))
# print("yaye version")
# print(yaye_version$I)

# yaye version new
yaye_version <- as.data.frame(lsoda(initial_values_yaye_version, times, Abbreviated_Model, params))
print("yaye version new")
print(yaye_version$I0)

# summer version
summer_version <- EpiModel$new(times, names(initial_values), as.list(initial_values), params,
                              list(c("infection_frequency", "beta", "S", "L1"),
                                   c("infection_frequency", "beta_reinfection", "L2", "L1"),
                                   c("standard_flows", "kappa", "L1", "L2"),
                                   c("standard_flows", "epsilon", "L1", "I"),
                                   c("standard_flows", "nu", "L2", "I"),
                                   c("standard_flows", "gamma", "I", "L2"),
                                   c("standard_flows", "delta", "I", "S"),
                                   c("compartment_death", "mui", "I")),
                              infectious_compartment="I", initial_conditions_sum_to_total = FALSE, report_progress = FALSE, reporting_sigfigs = 6,
                              birth_approach = "replace_deaths", entry_compartment = "S")

summer_version$stratify("infect", seq(0, 2), c("I"), 
                        list(epsilon=list(adjustments=list("0"=params$prop_I0,
                                                           "1"=params$prop_I1,
                                                           "2"=params$prop_I2)),
                             nu=list(adjustments=list("0"=params$prop_I0,
                                                      "1"=params$prop_I1,
                                                      "2"=params$prop_I2))),
                        infectiousness_adjustments = c("0"=params$alpha_0,
                                                       "1"=params$alpha_1,
                                                       "2"=1),
                        report = FALSE)


summer_version$add_transition_flow(c("standard_flows", "h", "IXinfect_0", "IXinfect_1"))
summer_version$add_transition_flow(c("standard_flows", "j", "IXinfect_1", "IXinfect_2"))


# print(summer_version$flows)

summer_version$run_model()
print("summer version")

# print(summer_version$outputs$I)


print(summer_version$outputs$IXinfect_0)




