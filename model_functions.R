
# Planning to move all the functions that represent the acgtual model to be run to this file - to avoid duplication

YayeModel <- function(current_timepoint, state_values, parameters)
  {
  # create state variables (local variables)
  S=state_values[1] # fully susceptible 
  L1=state_values[2] # early latency
  L2=state_values[3] # late latency
  I0=state_values[4] # infectious non-spreaders
  I1=state_values[5] # infectious non-spreaders
  I2=state_values[6] # infectious non-spreaders
    
  N <- sum(state_values)

    with(
      as.list(parameters), #variable names within parameters can be used
      {
        
        #compute derivative
        dS = mu * N + 
          mui0 * I0 + mui1 * I1 + mui2 * I2 +
          delta * (I0 + I1 + I2) -
          beta * S * (I0 * alpha_0 + I1 * alpha_1 + I2) / N - 
          mu * S
        
        dL1 = (beta * S + beta * r * L2) * (I0 * alpha_0 + I1 * alpha_1 + I2) / N - 
          (epsilon + kappa + mu) * L1
        
        dL2 = kappa * L1 + 
          gamma0 * I0 + gamma1 * I1 + gamma2 * I2 - 
          (beta * r * (I0 * alpha_0 + I1 * alpha_1 + I2) + nu + mu) * L2
        
        dI0 = epsilon * L1 * prop_I0 + 
          nu * L2 * prop_I0 - 
          (gamma0 + delta + h + mui0 + mu) * I0
        
        dI1 = epsilon * L1 * prop_I1 + 
          nu * L2 * prop_I1 +
          h * I0 - 
          (gamma1 + delta + j + mui1 + mu) * I1
        
        dI2 = epsilon * L1 * prop_I2 + 
          nu * L2 * prop_I2 + 
          j * I1 - 
          (gamma2 + delta + mui2 + mu) * I2

        #combine results
        results = c(dS, dL1, dL2, dI0, dI1, dI2)
        list(results)
      }
    )
  }
  
# summer version
summer_version <- EpiModel$new(times, names(initial_values), as.list(initial_values), params,
                               list(c("infection_frequency", "beta", "S", "L1"),
                                    c("infection_frequency", "beta_reinfection", "L2", "L1"),
                                    c("standard_flows", "kappa", "L1", "L2"),
                                    c("standard_flows", "epsilon", "L1", "I"),
                                    c("standard_flows", "nu", "L2", "I"),
                                    c("standard_flows", "gamma", "I", "L2"),
                                    c("standard_flows", "delta", "I", "S"),
                                    c("compartment_death", "mui0", "I")),
                               infectious_compartment="I", initial_conditions_sum_to_total = FALSE, report_progress = FALSE, reporting_sigfigs = 6,
                               birth_approach = "replace_deaths", entry_compartment = "S")

# stratify by infectiousness and modify relevant parameters
summer_version$stratify("infect", seq(0, 2), c("I"),
                        list(epsilon=list(adjustments=list("0"=params$prop_I0,
                                                           "1"=params$prop_I1,
                                                           "2"=params$prop_I2)),
                             nu=list(adjustments=list("0"=params$prop_I0,
                                                      "1"=params$prop_I1,
                                                      "2"=params$prop_I2)),
                             gamma=list(adjustments=list("0"=params$gamma0,
                                                         "1"=params$gamma1,
                                                         "2"=params$gamma2),
                                        overwrite=seq(0, 2)),
                             mui=list(adjustments=list("0"=params$mui0,
                                                       "1"=params$mui1,
                                                       "2"=params$mui2),
                                      overwrite=seq(0, 2))),
                        infectiousness_adjustments = c("0"=params$alpha_0,
                                                       "1"=params$alpha_1,
                                                       "2"=1),
                        report = FALSE)

# add the between-infectious compartment flows
summer_version$add_transition_flow(c("standard_flows", "h", "IXinfect_0", "IXinfect_1"))
summer_version$add_transition_flow(c("standard_flows", "j", "IXinfect_1", "IXinfect_2"))


