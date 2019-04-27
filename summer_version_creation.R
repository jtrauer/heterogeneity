
# specification of summer version of model
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
                               birth_approach = "replace_deaths", entry_compartment = "S", equilibrium_stopping_tolerance = tolerance,
                               output_connections = list(incidence = c(from = "", to = "I")))

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