
library(deSolve)
library(reshape2)
# library(ggplot2)
# library(ggpubr)
library(lhs)
library(sensitivity)

source("function_tool_kit.R")
source("model_functions.R")

# SENSITIVITY ANALYSIS

# beginning to shift this code over to being based on lists, to avoid repeated calculations using LHS to sample from prop_I0 particular window
param_value_limits <- list(N = list(min = 1, max = 1),
                           prop_I0 = list(min = 0.58, max = 0.58),
                           prop_I1 = list(min = 0.31, max = 0.31),
                           prop_I2 = list(min = 0.11, max = 0.11),
                           P_epsilon = list(min = 0.074, max = 0.128),
                           P_nu = list(min = 0.018, max = 0.077), 
                           alpha = list(min = 0.22, max = 0.22),
                           mu = list(min = 1 / 55, max = 1 / 75),
                           Time_L1 = list(min = 0.167, max = 1.0),
                           Time_L2 = list(min = 20, max = 20),
                           Time_I = list(min = 3, max = 3),
                           P_mui0 = list(min = 0.05, max = 0.096),
                           P_mui1 = list(min = 0.05, max = 0.096),
                           P_mui2 = list(min = 0.327, max = 0.787),
                           r = list(min = 0.21, max = 0.21),
                           beta0 = list(min = 0, max = 0),
                           beta2 = list(min = 40, max = 60),
                           cdr_b = list(min = 0.5, max = 0.8), # baseline CDR
                           treatment_success = list(min = 0.8, max = 0.8),
                           P_h = list(min = 0.0, max = 0.058),
                           P_j = list(min = 0.0, max = 0.058),
                           p1 = list(min = 0, max = 0),
                           p2 = list(min = 0, max = 0),
                           alpha_1 = list(min = 0.22, max = 0.22),
                           alpha_0 = list(min = 0, max= 0))

# Latin hypercube sampling
n_runs <- 2 # choose number of points to simulate
set.seed(6242015) # random number generator
# To map these points in the unit cube to our parameters, we need minimum and maximum values for each.

lhs <- maximinLHS(n_runs, length(param_value_limits)) # simulate h= number of simulations

# Initialise empty data frame with correct number of rows to store parameter values
params_matrix <- data.frame(matrix(NA, nrow = n_runs, ncol = 0))

# Populate the parameter data frame sequentially by column from lhs sampling
parameter_names <- names(param_value_limits)
for (parameter in seq(length(parameter_names))) {
  param_name <- parameter_names[parameter]
  params_matrix[param_name] <- 
    adjust_lhs_to_range(lhs[, parameter], param_name, param_value_limits)
}


# Add colums for parameters and outputs
output_matrix_equi <- params_matrix
derived_params <- c("beta1", "epsilon0", "epsilon1", "epsilon2", "kappa", "nu0", "nu1", "nu2", "mui0", "mui1", "mui2", "gamma0",
                    "gamma1", "gamma2", "h", "j", "delta0_b", "delta1_b", "delta2_b", "Equi_incidence")
for (i in derived_params) {
  output_matrix_equi[[i]] <- NA
}

# compute beta1
output_matrix_equi$beta1 = output_matrix_equi$beta2 * output_matrix_equi$alpha

# compute parameters that are derived from proportions and sojourn times
output_matrix_equi$epsilon0 <- find_rate_from_proportion_params("P_epsilon", "Time_L1", output_matrix_equi, "prop_I0")
output_matrix_equi$epsilon1 <- find_rate_from_proportion_params("P_epsilon", "Time_L1", output_matrix_equi, "prop_I1")
output_matrix_equi$epsilon2 <- find_rate_from_proportion_params("P_epsilon", "Time_L1", output_matrix_equi, "prop_I2")
output_matrix_equi$epsilon <- find_rate_from_proportion_params("P_epsilon", "Time_L1", output_matrix_equi)
output_matrix_equi$nu0 <- find_rate_from_proportion_params("P_nu", "Time_L2", output_matrix_equi, "prop_I0")
output_matrix_equi$nu1 <- find_rate_from_proportion_params("P_nu", "Time_L2", output_matrix_equi, "prop_I1")
output_matrix_equi$nu2 <- find_rate_from_proportion_params("P_nu", "Time_L2", output_matrix_equi, "prop_I2")
output_matrix_equi$nu <- find_rate_from_proportion_params("P_nu", "Time_L2", output_matrix_equi)
output_matrix_equi$mui0 <- find_rate_from_proportion_params("P_mui0", "Time_I", output_matrix_equi)
output_matrix_equi$mui1 <- find_rate_from_proportion_params("P_mui1", "Time_I", output_matrix_equi)
output_matrix_equi$mui2 <- find_rate_from_proportion_params("P_mui2", "Time_I", output_matrix_equi)
output_matrix_equi$h <- find_rate_from_proportion_params("P_h", "Time_I", output_matrix_equi)
output_matrix_equi$j <- find_rate_from_proportion_params("P_j", "Time_I", output_matrix_equi)

# other parameter calculations that have now changed
output_matrix_equi$kappa <- (1 - output_matrix_equi$P_epsilon) / output_matrix_equi$Time_L1
output_matrix_equi$gamma0 <- 1 / output_matrix_equi$Time_I - (output_matrix_equi$mu + output_matrix_equi$h + output_matrix_equi$mui0)
output_matrix_equi$gamma1 <- 1 / output_matrix_equi$Time_I - (output_matrix_equi$mu + output_matrix_equi$j + output_matrix_equi$mui1)
output_matrix_equi$gamma2 <- 1 / output_matrix_equi$Time_I - (output_matrix_equi$mui2 + output_matrix_equi$mu)

# compute baseline deltas
output_matrix_equi$delta0_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b, 
                      output_matrix_equi$gamma0 + output_matrix_equi$mui0 +
                        output_matrix_equi$mu + output_matrix_equi$h, 
                      output_matrix_equi$treatment_success)
output_matrix_equi$delta1_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b,
                      output_matrix_equi$gamma1 + output_matrix_equi$mui1 + 
                        output_matrix_equi$mu + output_matrix_equi$j, 
                      output_matrix_equi$treatment_success)
output_matrix_equi$delta2_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b,
                      output_matrix_equi$gamma2 + output_matrix_equi$mui2 + 
                        output_matrix_equi$mu,
                      output_matrix_equi$treatment_success)

#
output_matrix_equi$delta <- output_matrix_equi$delta0_b
output_matrix_equi$beta <- output_matrix_equi$beta2

# Initial conditions for each compartment: 
S_init = 1
L1_init = 0
L2_init = 0
infectious_seed = 1e-6
I0_init = infectious_seed * output_matrix_equi$prop_I0[1]
I1_init = infectious_seed * output_matrix_equi$prop_I1[1]
I2_init = infectious_seed * output_matrix_equi$prop_I2[1]

initial_model_run_duration <- 2e3
additional_model_run_durations <- 5e2
equilibrium_prevalence_threshold <- 1e-6

# Loop up to equilibrium
for (run in seq(n_runs)) {
  
  # Run baseline
  initial_values = c(S = S_init - I0_init - I1_init - I2_init,
                     L1 = L1_init, L2 = L2_init, I0 = I0_init, I1 = I1_init, I2 = I2_init)
  params <- as.list(c(output_matrix_equi[run, ]))
  times <- seq(0, initial_model_run_duration)
  
  baseline_output <- as.data.frame(lsoda(initial_values, times, YayeModel, params))
  
  population_size <- rowSums(baseline_output[, 2:7])

  baseline_output$prevalence <- baseline_output$I0 + baseline_output$I1 + baseline_output$I2
  prevalence_baseline_output <- diff(baseline_output$prevalence) / population_size[-1] * 1e5
  prevalence_change <- abs(prevalence_baseline_output[length(prevalence_baseline_output) - 1] - 
                             prevalence_baseline_output[length(prevalence_baseline_output) - 2])
  
  # Loop over serial periods of time while checking to see if equilibrium has been reached
  while(prevalence_change > equilibrium_prevalence_threshold){
    initial_values <- c(S = min(baseline_output$S),
                        L1 = max(baseline_output$L1), L2 = max(baseline_output$L2), I0 = max(baseline_output$I0), I1 = max(baseline_output$I1), I2 = max(baseline_output$I2),
                        inc = 0)
    times=seq(0, additional_model_run_durations)
    baseline_output <- as.data.frame(lsoda(initial_values, times, YayeModel, params))
    population_size <- rowSums(baseline_output[, 2:7])
    incidence_baseline_output <- diff(baseline_output$inc) / population_size[-1] * 1e5
    incidence_change <- abs(incidence_baseline_output[length(incidence_baseline_output) - 1] - 
                     incidence_baseline_output[length(incidence_baseline_output) - 2])
  }
  end_baseline_prevalence <- tail(prevalence_baseline_output, n=1)
  output_matrix_equi$Equi_prevalence[run] <- end_baseline_prevalence
}

#now we have incidence and each parameters in one matrix

#save output as csv

write.csv(x=output_matrix_equi,file='..//heterogeneity/output/output_matrix_equi.csv')

#test significance of correlation, select parameters from csv manually 

Test_corelation<-read.csv(file.choose())#open csv with selected params and incidence
View(Test_corelation)
bonferroni.alpha <- 0.05/11 #11 parametrs

prcc <- pcc(Test_corelation[,1:11], Test_corelation[,12], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')
load('prcc.Rdata')
summary <- print(prcc)
Corl=data.frame(summary) 
View(Corl)
par(mar=c(9,4,4,2)+0.1)
plot(Corl$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',pch=19,
     axes=FALSE)
axis(2)
axis(1, at=seq(1:14), labels=row.names(Corl), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:11) lines(c(i,i),c(Corl[i,4], Corl[i,5]))
abline(h=0)

#ggplot for parameter space Vs equilibrium incidence
#significan colour= coral4, non-significant chocolate1 
#CDR, significant
e_cdr_b=data.frame(cbind(output_matrix_equi$cdr_b, output_matrix_equi$Equi_incidence))

e_cdr_b_g=melt(e_cdr_b,id='X1')
e_cdr_b_gp<-ggplot(e_cdr_b_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  
  xlab("CDR")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_cdr_b_gp

#nu
e_nu=data.frame(cbind(output_matrix_equi$P_nu, output_matrix_equi$Equi_incidence))
#View(e_nu)
e_nu_g=melt(e_nu,id='X1')
e_nu_gp<-ggplot(e_nu_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  
  xlab("nu")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_nu_gp

#beta2

e_beta2=data.frame(cbind(output_matrix_equi$beta2,output_matrix_equi$Equi_incidence))
e_beta2_g=melt(e_beta2,id='X1')
e_beta2_gp<-ggplot(e_beta2_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("beta")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_beta2_gp
#epsiolon
e_epsilon=data.frame(cbind(output_matrix_equi$P_epsilon,output_matrix_equi$Equi_incidence))
e_epsilon_g=melt(e_epsilon,id='X1')
e_epsilon_gp<-ggplot(e_epsilon_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("epsilon")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_epsilon_gp
#Time L1
e_T_L1=data.frame(cbind(output_matrix_equi$Time_L1,output_matrix_equi$Equi_incidence))
e_TL1_g=melt(e_T_L1,id='X1')
e_TL1_gp<-ggplot(e_TL1_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  xlab("Time L1")+
  # theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_TL1_gp
#j
e_j=data.frame(cbind(output_matrix_equi$P_j,output_matrix_equi$Equi_incidence))
e_j_g=melt(e_j,id='X1')
e_j_gp<-ggplot(e_j_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  xlab("j")+
  # theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_j_gp
#mu
e_mu=data.frame(cbind(output_matrix_equi$mu,output_matrix_equi$Equi_incidence))

e_mu_g=melt(e_mu,id='X1')
e_mu_gp<-ggplot(e_mu_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("mu")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mu_gp

#mui2
e_mui2=data.frame(cbind(output_matrix_equi$P_mui2,output_matrix_equi$Equi_incidence))
e_mui2_g=melt(e_mui2,id='X1')
e_mui2_gp<-ggplot(e_mui2_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("mui2")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mui2_gp

#mui1
e_mui1=data.frame(cbind(output_matrix_equi$P_mui1,output_matrix_equi$Equi_incidence))

e_mui1_g=melt(e_mui1,id='X1')
e_mui1_gp<-ggplot(e_mui1_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("mui1")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mui1_gp
#mui0
e_mui0=data.frame(cbind(output_matrix_equi$P_mui0,output_matrix_equi$Equi_incidence))
e_mui0_g=melt(e_mui0,id='X1')
e_mui0_gp<-ggplot(e_mui0_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("mui0")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mui0_gp
#h
e_h=data.frame(cbind(output_matrix_equi$P_h,output_matrix_equi$Equi_incidence))
e_h_g=melt(e_h,id='X1')
e_h_gp<-ggplot(e_h_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab("h")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_h_gp



GG<-ggarrange(e_cdr_b_gp,e_beta2_gp,
              e_epsilon_gp,e_nu_gp,e_TL1_gp,e_j_gp,e_mui2_gp,
              e_mui1_gp,e_mui0_gp,e_mu_gp,e_h_gp,
              ncol=3, nrow=4, common.legend = TRUE,
              legend="bottom")

annotate_figure(GG,left = text_grob("Equilibrium incidence per 100,000popn",rot = 90))
#



 









