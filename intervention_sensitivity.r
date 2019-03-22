#not yet done
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
require(lhs)
library(sensitivity)
source("function_tool_kit.R")



Baseline_model=function(current_timepoint, state_values, parameters)
{
  #creat state variables (local variables)
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
      dS=(mu*N+mui0*I0+mui1*I1+mui2*I2)+(delta0_b*I0+delta1_b*I1+delta2_b*I2)+((r*(p1*beta1*I1+p2*beta2*I2)/N)*L2)-((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S-mu*S
      dL1=((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S+(r*(beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*L2-(epsilon0+epsilon1+epsilon2+mu+kappa)*L1
      dL2=kappa*L1+(gamma0*I0+gamma1*I1+gamma2*I2)-(r*(beta0*I0+beta1*I1+beta2*I2)/N)*L2-(nu0+nu1+nu2+mu)*L2
      dI0=(epsilon0*L1+nu0*L2)-(mui0+mu+gamma0+delta0_b+h)*I0
      dI1=(epsilon1*L1+nu1*L2+h*I0)-(mui1+mu+gamma1+delta1_b+j)*I1
      dI2=(epsilon2*L1+nu2*L2+j*I1)-(mui2+mu+gamma2+delta2_b)*I2
      
      dinc=(epsilon0*L1+nu0*L2)+(epsilon1*L1+nu1*L2)+(epsilon2*L1+nu2*L2)
      
      #combine results
      results=c(dS, dL1, dL2, dI0,dI1,dI2, dinc)
      list(results)
    }
  )
}

#Intervention model

intervention_model=function(current_timepoint, state_values, parameters)
{
  
  S=state_values[1] #fully susceptible 
  L1=state_values[2] #early latency
  L2=state_values[3] #late latency
  I0=state_values[4] #infectious non-spreaders
  I1=state_values[5] #infectious spreaders
  I2=state_values[6] #infectious Super-spreaders
  
  #creat state variables (local variables)
  
  
  with(
    as.list(parameters), #variable names within parameters can be used
    {
      #compute derivative
      dS=(mu*N+mui0*I0+mui1*I1+mui2*I2)+(delta0_i*I0+delta1_i*I1+delta2_i*I2)+((r*(p1*beta1*I1+p2*beta2*I2)/N)*L2)-((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S-mu*S
      dL1=((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S+(r*(beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*L2-(epsilon0+epsilon1+epsilon2+mu+kappa)*L1
      dL2=kappa*L1+(gamma0*I0+gamma1*I1+gamma2*I2)-(r*(beta0*I0+beta1*I1+beta2*I2)/N)*L2-(nu0+nu1+nu2+mu)*L2
      dI0=(epsilon0*L1+nu0*L2)-(mui0+mu+gamma0+delta0_i+h)*I0
      dI1=(epsilon1*L1+nu1*L2+h*I0)-(mui1+mu+gamma1+delta1_i+j)*I1
      dI2=(epsilon2*L1+nu2*L2+j*I1)-(mui2+mu+gamma2+delta2_i)*I2
      
      dinc=(epsilon0*L1+nu0*L2)+(epsilon1*L1+nu1*L2)+(epsilon2*L1+nu2*L2)
      
      #combine results
      results=c(dS, dL1, dL2, dI0,dI1,dI2, dinc)
      list(results)
    }
  )
}

# Time in each compartment
T_L1=1/4#time in L1 in years
T_L2=20#time in L2 in years
T_I=3#Time in I in years


#Latin hypercube sampling
z <- 5 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(z,26) #simulate h= number of simulations, 35=number of parameters
#To map these points in the unit cube to our parameters, we need minimum and maximum values for each.

# beginning to shift this code over to being based on lists, to avoid repeated calculations using LHS to sample from a particular window
param_value_limits <- list(N= list(min = 1, max = 1),
                           a = list(min = 0.58, max = 0.58),
                           b = list(min = 0.31, max = 0.31),
                           c = list(min = 0.11, max = 0.11),
                           alpha = list(min = 0.22, max = 0.22),
                           mu = list(min = 0.0133, max = 0.0182),
                           P_mui0 = list(min = 0.049, max = 0.091),
                           P_mui1 = list(min = 0.049, max = 0.091),
                           P_mui2 = list(min = 0.279, max = 0.544),
                           r = list(min = 0.21, max = 0.21),
                           beta0 = list(min = 0, max = 0),
                           beta2 = list(min = 30, max = 60),
                           P_epsilon = list(min = 0.074, max = 0.128),
                           P_kappa = list(min = 0.54, max = 0.721),
                           P_gamma0 = list(min = 0.548, max = 0.595),
                           P_gamma1 = list(min = 0.548, max = 0.595),
                           P_gamma2 = list(min = 0.155, max = 0.466),
                           P_nu = list(min = 0.018, max = 0.077),
                           cdr_b = list(min = 0.5, max = 0.8),#baseline CDR
                           s = list(min = 0.8, max = 0.8),#Rx success
                           P_h = list(min = 0.0, max = 0.058),
                           P_j = list(min = 0.0, max = 0.058),
                           p1 = list(min = 0, max = 0),
                           p2 = list(min = 0, max = 0),
                           q = list(min = 0, max = 0.75),
                           d = list(min = 0, max = 0.75))


#Now we can generate a ?parameter set? by rescaling our simulated latin hypercube sample
params.set_Sint <- cbind(
  a = adjust_lhs_to_range(lhs[, 1], "a", param_value_limits),
  b = adjust_lhs_to_range(lhs[, 2], "b", param_value_limits),
  c = adjust_lhs_to_range(lhs[, 3], "c", param_value_limits),
  alpha = adjust_lhs_to_range(lhs[, 4], "alpha", param_value_limits),
  mu = adjust_lhs_to_range(lhs[, 5], "mu", param_value_limits),
  P_mui0 = adjust_lhs_to_range(lhs[, 6], "P_mui0", param_value_limits),
  P_mui1 = adjust_lhs_to_range(lhs[, 7], "P_mui1", param_value_limits),
  P_mui2 = adjust_lhs_to_range(lhs[, 8], "P_mui2", param_value_limits),
  r = adjust_lhs_to_range(lhs[, 9], "r", param_value_limits),
  beta0 = adjust_lhs_to_range(lhs[, 10], "beta0", param_value_limits),
  beta2 = adjust_lhs_to_range(lhs[, 11], "beta2", param_value_limits),
  P_epsilon = adjust_lhs_to_range(lhs[, 12], "P_epsilon", param_value_limits),
  P_kappa = adjust_lhs_to_range(lhs[, 13], "P_kappa", param_value_limits),
  P_gamma0 = adjust_lhs_to_range(lhs[, 14], "P_gamma0", param_value_limits),
  P_gamma1 = adjust_lhs_to_range(lhs[, 15], "P_gamma1", param_value_limits),
  P_gamma2 = adjust_lhs_to_range(lhs[, 16], "P_gamma2", param_value_limits),
  P_nu = adjust_lhs_to_range(lhs[, 17], "P_nu", param_value_limits),
  cdr_b = adjust_lhs_to_range(lhs[, 18], "cdr_b", param_value_limits),
  s = adjust_lhs_to_range(lhs[, 19], "s", param_value_limits),
  P_h = adjust_lhs_to_range(lhs[, 20], "P_h", param_value_limits),
  P_j = adjust_lhs_to_range(lhs[, 21], "P_j", param_value_limits),
  p1 = adjust_lhs_to_range(lhs[, 22], "p1", param_value_limits),
  p2 = adjust_lhs_to_range(lhs[, 23], "p2", param_value_limits),
  N= adjust_lhs_to_range(lhs[, 24], "N", param_value_limits),
  q= adjust_lhs_to_range(lhs[, 25], "q", param_value_limits),
  d= adjust_lhs_to_range(lhs[, 26], "d", param_value_limits))

View(params.set_Sint)


#creat matrix to save whole info
params_matrix_Sint = data.frame(params.set_Sint)
#add colums for parameters 
beta1 = data.frame('beta1'=rep(NA,z))

epsilon0 = data.frame('epsilon0'=rep(NA,z))
epsilon1 = data.frame('epsilon1'=rep(NA,z))
epsilon2 = data.frame('epsilon2'=rep(NA,z))

kappa = data.frame('kappa'=rep(NA,z))

nu0 = data.frame('nu0'=rep(NA,z))
nu1 = data.frame('nu1'=rep(NA,z))
nu2 = data.frame('nu2'=rep(NA,z))

mui0 = data.frame('mui0'=rep(NA,z))
mui1 = data.frame('mui1'=rep(NA,z))
mui2 = data.frame('mui2'=rep(NA,z))

gamma0 = data.frame('gamma0'=rep(NA,z))
gamma1 = data.frame('gamma1'=rep(NA,z))
gamma2 = data.frame('gamma2'=rep(NA,z))

h = data.frame('h'=rep(NA,z))
j = data.frame('j'=rep(NA,z))

CDR_inter=data.frame('CDR_inter'=rep(NA,z))#intervention level CDR
B_incidence = data.frame('B_incidence'=rep(NA,z))#baseline equilibrium incidence
Int_incidence = data.frame('Int_incidence'=rep(NA,z))#10years intervention incidence
Change_incidence= data.frame('Change_incidence'=rep(NA,z))#absolute incidence difference
Relative_change=data.frame('Relative_change'=rep(NA,z))#relative incidence difference
#add colums for delta base line


params_matrix_Sint=cbind(params_matrix_Sint,beta1,epsilon0,epsilon1,epsilon2,
                         mui0,mui1,mui2,kappa, nu0,nu1,nu2,gamma0,gamma1,gamma2,h,j,
                         CDR_inter,B_incidence,Int_incidence,Change_incidence,
                         Relative_change)
View(params_matrix_Sint)



#View(output_matrix_equi)
#compute beta1
params_matrix_Sint$beta1=(params_matrix_Sint$beta2)*(params_matrix_Sint$alpha)
#compute epsilon
params_matrix_Sint$epsilon0=log(1-(params_matrix_Sint$P_epsilon*params_matrix_Sint$a))/(-T_L1)
params_matrix_Sint$epsilon1=log(1-(params_matrix_Sint$P_epsilon*params_matrix_Sint$b))/(-T_L1)
params_matrix_Sint$epsilon2=log(1-(params_matrix_Sint$P_epsilon*params_matrix_Sint$c))/(-T_L1)
#compute kappa
params_matrix_Sint$kappa=log(1-(params_matrix_Sint$P_kappa))/(-T_L1)
#compute nu
params_matrix_Sint$nu0=log(1-(params_matrix_Sint$P_nu*params_matrix_Sint$a))/(-T_L2)
params_matrix_Sint$nu1=log(1-(params_matrix_Sint$P_nu*params_matrix_Sint$b))/(-T_L2)
params_matrix_Sint$nu2=log(1-(params_matrix_Sint$P_nu*params_matrix_Sint$c))/(-T_L2)
#compute gamma
params_matrix_Sint$gamma0=log(1-(params_matrix_Sint$P_gamma0))/(-T_I)
params_matrix_Sint$gamma1=log(1-(params_matrix_Sint$P_gamma1))/(-T_I)
params_matrix_Sint$gamma2=log(1-(params_matrix_Sint$P_gamma2))/(-T_I)
#compute mui
params_matrix_Sint$mui0=log(1-(params_matrix_Sint$P_mui0))/(-T_I)
params_matrix_Sint$mui1=log(1-(params_matrix_Sint$P_mui1))/(-T_I)
params_matrix_Sint$mui2=log(1-(params_matrix_Sint$P_mui2))/(-T_I)
#compute h and j
params_matrix_Sint$h=log(1-(params_matrix_Sint$P_h))/(-T_I)
params_matrix_Sint$j=log(1-(params_matrix_Sint$P_j))/(-T_I)



View(params_matrix_Sint)
##Initial values for sub population: 
a=0.58 #proportion for I0
b=0.31#proportion for I1
c=0.11#proportion for I2
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6*a    #active TB hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # active TB Normal spreaders
F=1e-6*c # active TB Super-spreader 10% of TB patients are super-spreaders
G=0 # diagnosed 



for(i in 1:z){
  #run baseline 
  params <- as.list(c(params_matrix_Sint[i,]))
  cdr0_b=params$cdr_b#baseline CDR
  cdr1_b=params$cdr_b
  cdr2_b=params$cdr_b
  #HERE I Am
  CDR_i_posibl=params$cdr_b+(1-params$cdr_b)*(params$q)#possible intervention size
  CDR_max=0.9#maximum possible CDR that the health system acan do
  CDR_i=min(CDR_i_posibl,CDR_max)#intervention CDR that can be used in model
  params_matrix_Sint$CDR_inter[i]=CDR_i#fill the parameter matrix
  
  cdr2_imax_possible=(CDR_i-(params$a+params$b)*params$CDR_b)/params$c#super-spreader CDR intervention 
  cdr2_imax=min(cdr2_imax_possible,CDR_max)#the limit that SS CDR can take
  
  cdr2_i=CDR_i+((cdr2_imax-CDR_i)*params$d)# Super-Spreaders CDR to be used in the model
  cdr1_i=(CDR_i-(params$c*cdr2_i))/(params$a+params$b)#low-spreaders CDR to be used in the model
  cdr0_i=cdr1_i#non-spreaders CDR
  #compute deltas for the Baseline_model
  delta0_b=(cdr0_b*(params$gamma0+params$mui0+params$mu+params$h)/(1-cdr0_b))
  delta1_b=(cdr1_b*(params$gamma1+params$mui1+params$mu+params$j)/(1-cdr1_b))
  delta2_b=(cdr2_b*(params$gamma2+params$mui2+params$mu)/(1-cdr2_b))
  #compute deltas for the BIntervention_model
  delta0_i=(cdr0_i*(params$gamma0+params$mui0+params$mu+params$h)/(1-cdr0_i))
  delta1_i=(cdr1_i*(params$gamma1+params$mui1+params$mu+params$j)/(1-cdr1_i))
  delta2_i=(cdr2_i*(params$gamma2+params$mui2+params$mu)/(1-cdr2_i))
  
  #RUN Baseline Model
  initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
  
  times=seq(0, 5000, by = 1)
  
  B_out <- as.data.frame(lsoda(initial_values, times, Baseline_model, params))
  
  #Record Baseline equilibrium incidence
  Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
  
  incidence_B_out=(diff(B_out$inc)/Nq)*100000
  
  B_Incidence_time_n = incidence_B_out[4999] #the model reaches equilibrium at time around 500
  output_matrix_j$B_incidence[i]=B_Incidence_time_n
  #record compartment values at equilibrium
  Equi_initial_values=c(S=min(B_out$S), L1=max(B_out$L1), L2=max(B_out$L2), I0=max(B_out$I0),I1=max(B_out$I1),I2=max(B_out$I2),inc=max(B_out$I0)+max(B_out$I1)+max(B_out$I2))
  
  #RUN Intervention model, 
  Int_out <- as.data.frame(lsoda(Equi_initial_values, times, intervention_model, params))
  
  Nq=Int_out$S+Int_out$L1+Int_out$L2+Int_out$I0+Int_out$I1+Int_out$I2
  Int_incidence_out=(diff(Int_out$inc)/Nq)*100000
  
  #considering a 10 years impact of intervention
  Int_Incidence_time_n = Int_incidence_out[10]
  output_matrix_j$Int_incidence[i]=Int_Incidence_time_n
  #calculate the difference between baseline incidence and intervention incidence
  Change_inc=B_Incidence_time_n-Int_Incidence_time_n 
  output_matrix_j$Change_incidence[i]=Change_inc
  Relative_inc=((B_Incidence_time_n-Int_Incidence_time_n)/B_Incidence_time_n)*100
  output_matrix_j$Relative_change[i]=Relative_inc
}

View(output_matrix_j)
write.csv(x=output_matrix_j,file = '..//Enchik.com/intervention/output_matrix_j.csv')


#View(output_matrix_equi) #now we have incidence and each parameters in one matrix

#save output as csv

write.csv(x=output_matrix_equi,file='..//heterogeneity/output/output_matrix_equi.csv')




