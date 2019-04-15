
source("model_functions.R")
source("summer_model.R")

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
            beta0 = 0,
            beta2 = 60,
            cdr_b = 0.8,
            treatment_success = 0.8,
            P_h = 0.058,
            P_j = 0.058,
            p1 = 0,
            p2 = 0,
            crude_birth_rate = 1 / 75,
            place_holder = .1)

S_init = 1
L1_init = 0
L2_init = 0
infectious_seed = 1e-6
prop_I0 = 0.58
I0_init = infectious_seed * prop_I0
I1_init = infectious_seed * prop_I1
I2_init = infectious_seed * prop_I2
initial_values = c(S = S_init - I0_init - I1_init - I2_init, L1 = 0, L2 = 0, I0 = I0_init, I1 = 0, I2 = 0)
initial_model_run_duration <- 1e2
times <- seq(0, initial_model_run_duration)
baseline_output <- as.data.frame(lsoda(initial_values, times, Abbreviated_Model, params))

print(baseline_output$I0)

summer_output <- EpiModel$new(times, names(initial_values), as.list(initial_values), params,
                              list(c("infection_frequency", "beta2", "S", "I0")),
                              infectious_compartment="I0", initial_conditions_sum_to_total = FALSE, report_progress = FALSE, 
                              birth_approach = "add_crude_birth_rate", entry_compartment = "S")

summer_output$run_model()

print(summer_output$outputs$I0)




