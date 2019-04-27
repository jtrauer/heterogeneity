
# setwd("C:/Users/jtrauer/Desktop/summer")
setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/summer")
source("summer_model.R")

setwd("//ad.monash.edu/home/User096/jtrauer/Documents/GitHub/heterogeneity")
# setwd("C://Users/jtrauer/Desktop/heterogeneity/heterogeneity")
source("function_tool_kit.R")
source("model_functions.R")
library(lhs)
local_repo_directory <- getwd()
setwd(local_repo_directory)

# process the parameters from the epidemiological parameters being varied into model flow rates
process_parameters <- function(parameters) {
  
  # latency progression parameters
  parameters$epsilon <- parameters$P_epsilon / parameters$Time_L1
  parameters$kappa <- (1 - parameters$P_epsilon) / parameters$Time_L1
  parameters$nu <- parameters$P_nu / parameters$Time_L2

  # between-disease compartment progression parameters
  parameters$j <- parameters$P_j / parameters$Time_I
  parameters$h <- parameters$P_h / parameters$Time_I
  
  # death and spontaneous recovery-related parameters
  parameters$universal_death_rate <- parameters$mu
  parameters$mui0 <- parameters$P_mui0 / parameters$Time_I
  parameters$mui1 <- parameters$P_mui1 / parameters$Time_I
  parameters$mui2 <- parameters$P_mui2 / parameters$Time_I
  parameters$gamma <- (1 - parameters$P_mui) / parameters$Time_I
  parameters$gamma0 <- 1 / parameters$Time_I - (parameters$mui0 + parameters$mu + parameters$h)
  parameters$gamma1 <- 1 / parameters$Time_I - (parameters$mui1 + parameters$mu + parameters$j)
  parameters$gamma2 <- 1 / parameters$Time_I - (parameters$mui2 + parameters$mu)

  # infection-related parameters
  parameters$beta_reinfection <- parameters$beta * parameters$r
  parameters$prop_I0 <- 1 - parameters$prop_I1 - parameters$prop_I2
  
  # case detection parameter
  parameters$delta <- find_delta_from_cdr(parameters$cdr_b, parameters$gamma + parameters$mu, 1)
  parameters
}

# first set the fixed parameter values
params <- list(Time_L1 = 1,
               Time_L2 = 20,
               Time_I = 3,
               P_mui = 0.787,
               treatment_success = 0.8,
               r = 0.21,
               prop_I1 = 0.31,
               prop_I2 = 0.11,
               alpha_1 = 0.22,
               alpha_0 = 0)

# now set the ones to vary through the LHS process
uncertainty_params <- list(P_epsilon = c(min = 0.074, max = 0.128),
                           P_nu = c(min = 0.018, max = 0.077),
                           Time_L1 = list(min = 0.167, max = 1.0),
                           mu = c(min = 1 / 75, max = 1 / 55),
                           P_mui0 = list(min = 0.05, max = 0.096),
                           P_mui1 = list(min = 0.05, max = 0.096),
                           P_mui2 = list(min = 0.327, max = 0.787),
                           beta = list(min = 40, max = 60),
                           cdr_b = list(min = 0.5, max = 0.8),
                           P_j = list(min = 0, max = 0.058),
                           P_h = list(min = 0, max = 0.058))

# user to request number of runs, set model intial conditions and specify integration time
n_runs <- 5
S_init = 1
I_init = 1e-6
initial_values = c(S = S_init - I_init, L1 = 0, L2 = 0, I = I_init)  # summer will split the I compartment
initial_values_yaye_version <- 
  c(S = S_init - I_init, L1 = 0, L2 = 0, I0 = I_init / 3, I1 = I_init / 3, I2 = I_init / 3, cumulative_incidence = 0)
times <- seq(0, 2e2)

# do the LHS sampling for all parameters and all runs
lhs_samples <- as.data.frame(maximinLHS(n_runs, length(uncertainty_params)))
colnames(lhs_samples) <- names(uncertainty_params)

# set stopping condition, based on largest absolute rate of change in compartment sizes being less than a specified value
tolerance <- 1e-5
stopping_condition <- function(t, state, parms) {
  max(abs(unlist(YayeModel(t, state, parms))[1:6])) - tolerance
}

# loop over the requested number of runs
for (run in seq(n_runs)) {
  writeLines(paste("\n____________\nrun number", run))

  # use the LHS sample values to determine parameter values and adapt parameters from epidemiological values to model flows
  for (param in names(uncertainty_params)) {
    params[[param]] <- adjust_lhs_to_range(lhs_samples[[param]][run], param, uncertainty_params)
  }
  params <- process_parameters(params)

  # create summer model and define yaye model
  source("summer_version_creation.R")

  # yaye version new
  yaye_version <- as.data.frame(lsodar(initial_values_yaye_version, times, YayeModel, params, rootfunc = stopping_condition))

  # find incidence
  yaye_version$time_diff <- c(0, diff(yaye_version$time))
  yaye_version$incidence <- c(0, diff(yaye_version$cumulative_incidence) * 1e5) / yaye_version$time_diff
  
  # cut out last row because adjusting the integration time mucks up the calculation of the last incidence value
  yaye_version <- yaye_version[1: nrow(yaye_version) - 1,]
    
  # run summer version
  summer_version$run_model()
  
  # print and report comparison of outputs
  plot(yaye_version$time, yaye_version$incidence, 
       xlab = "Time in years", ylab = "Overall incidence per 100,000 per year")
  lines(summer_version$incidence$times, summer_version$incidence$incidence * 1e5, "l", col = "blue")
  lines(summer_version$outputs$time, summer_version$outputs$IXinfect_2 * 1e5, "p")
  lines(yaye_version$time, yaye_version$I2 * 1e5, "l", col = "red", lwd = 2)
  writeLines(paste("direct ode-based version, prevalence of I0:", tail(yaye_version$I0, 1) * 1e5, 
                   "per 100,000"))
  writeLines(paste("summer interpretation, prevalence of I0:   ", tail(summer_version$outputs$IXinfect_0, 1) * 1e5, 
                   "per 100,000"))
  writeLines(paste("time to equilibrium:", tail(yaye_version$time, 1), "years"))
}

