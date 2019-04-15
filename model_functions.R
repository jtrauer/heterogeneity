
# Planning to move all the functions that represent the acgtual model to be run to this file - to avoid duplication

Abbreviated_Model <- function(current_timepoint, state_values, parameters)
  {
  # create state variables (local variables)
  S=state_values[1] #fully susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I0=state_values[4] #infectious non-spreaders
  I1=state_values[5] #infectious spreaders
  I2=state_values[6] #infectious Super-spreaders

    with(
      as.list(parameters), #variable names within parameters can be used
      {
        
        N <- sum(state_values)
        
        #compute derivative
        dS = mu * N - beta2 * S * I0 / N
        dL1 = 0
        dL2 = 0
        dI0 = beta2 * S * I0 / N
        dI1 = 0
        dI2 = 0

        #combine results
        results = c(dS, dL1, dL2, dI0, dI1, dI2)
        list(results)
      }
    )
  }
  


# Model Func
Baseline_model <- function(current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  S=state_values[1] #fully susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I0=state_values[4] #infectious non-spreaders
  I1=state_values[5] #infectious spreaders
  I2=state_values[6] #infectious Super-spreaders
  
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      dS = mu * N + mui0 * I0 + mui1 * I1 + mui2 * I2 +
        delta0_b * I0 + delta1_b * I1 + delta2_b * I2 +
        r * (p1 * beta1 * I1 + p2 * beta2 * I2) / N * L2 - 
        (beta0 * I0 + (1 - p1) * beta1 * I1 + (1 - p2) * beta2 * I2) / N * S -
        mu*S
      dL1 = (beta0 * I0 + (1 - p1) * beta1 * I1 + (1 - p2) * beta2 * I2) / N * S +
        r * (beta0 * I0 + (1 - p1) * beta1 * I1 + (1 - p2) * beta2 * I2) / N * L2 -
        (epsilon0 + epsilon1 + epsilon2 + mu + kappa) * L1
      dL2 = kappa * L1 + 
        gamma0 * I0 + gamma1 * I1 + gamma2 * I2 -
        r * (beta0 * I0 + beta1 * I1 + beta2 * I2) / N * L2 -
        (nu0 + nu1 + nu2 + mu) * L2
      dI0 = epsilon0 * L1 + nu0 * L2 -
        (mui0 + mu + gamma0 + delta0_b + h) * I0
      dI1 = epsilon1 * L1 + nu1 * L2 + h * I0 -
        (mui1 + mu + gamma1 + delta1_b + j) * I1
      dI2 = epsilon2 * L1 + nu2 * L2 + j * I1 -
        (mui2 + mu + gamma2 + delta2_b) * I2
      
      dinc = epsilon0 * L1 + nu0 * L2 + 
        epsilon1 * L1 + nu1 * L2 + epsilon2 * L1 + nu2 * L2
      
      #combine results
      results = c(dS, dL1, dL2, dI0, dI1, dI2, dinc)
      list(results)
    }
  )
}