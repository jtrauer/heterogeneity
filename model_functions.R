
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
  

