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

# beginning to shift this code over to being based on lists, to avoid repeated calculations using LHS to sample from a particular window
param_value_limits <- list(N= list(min = 1, max = 1),
                           a = list(min = 0.58, max = 0.58),
                           b = list(min = 0.31, max = 0.31),
                           c = list(min = 0.11, max = 0.11),
                           P_epsilon = list(min = 0.074, max = 0.128),
                           P_nu = list(min = 0.018, max = 0.077), 
                           alpha = list(min = 0.22, max = 0.22),
                           mu = list(min = 1/55, max = 1/75),
                           Time_L1=list(min = 0.167, max = 0.9),
                           Time_L2=list(min = 20, max = 20),
                           Time_I=list(min = 3, max = 3),
                           P_mui0 = list(min = 0.05, max = 0.096),
                           P_mui1 = list(min = 0.05, max = 0.096),
                           P_mui2 = list(min = 0.327, max = 0.787),
                           r = list(min = 0.21, max = 0.21),
                           beta0 = list(min = 0, max = 0),
                           beta2 = list(min = 40, max = 60),
                           cdr_b = list(min = 0.5, max = 0.8),#baseline CDR
                           s = list(min = 0.8, max = 0.8),#Rx success
                           P_h = list(min = 0.0, max = 0.058),
                           P_j = list(min = 0.0, max = 0.058),
                           q = list(min = 0, max = 0.75),
                           d = list(min = 0, max = 0.75),
                           p1 = list(min = 0, max = 0),
                           p2 = list(min = 0, max = 0)) 

# Latin hypercube sampling
z <- 10 # choose number of points to simulate
set.seed(6242015) # random number generator
# To map these points in the unit cube to our parameters, we need minimum and maximum values for each.

lhs <- maximinLHS(z, length(param_value_limits)) # simulate h= number of simulations

# Initialise empty data frame with correct number of rows to store parameter values
params_matrix_i <- data.frame(matrix(NA, nrow = z, ncol = 0))

# Populate the parameter data frame sequentially by column from lhs sampling
parameter_names <- names(param_value_limits)
for (parameter in seq(length(parameter_names))) {
  param_name <- parameter_names[parameter]
  params_matrix_i[param_name] <- 
    adjust_lhs_to_range(lhs[, parameter], param_name, param_value_limits)
}
View(params_matrix_i)


#add colums for parameters to be calculated from lhs sample
beta1 = data.frame('beta1'=rep(NA,z))
Total_L1=data.frame('Total_L1'=rep(NA,z))
Total_L2=data.frame('Total_L2'=rep(NA,z))
Total_I=data.frame('Total_I'=rep(NA,z))
epsilon = data.frame('epsilon'=rep(NA,z))
epsilon0 = data.frame('epsilon0'=rep(NA,z))
epsilon1 = data.frame('epsilon1'=rep(NA,z))
epsilon2 = data.frame('epsilon2'=rep(NA,z))

kappa = data.frame('kappa'=rep(NA,z))
nu = data.frame('nu'=rep(NA,z))
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



params_matrix_Sint=cbind(params_matrix_i,
                         beta1,Total_L1,Total_L2,Total_I,
                         epsilon,epsilon0,epsilon1,epsilon2,
                         mui0,mui1,mui2,kappa, 
                         nu,nu0,nu1,nu2,
                         gamma0,gamma1,gamma2,h,j,
                         CDR_inter,B_incidence,Int_incidence,
                         Change_incidence,
                         Relative_change)
View(params_matrix_Sint)
#View(output_matrix_equi)
#compute beta1
params_matrix_Sint$beta1=(params_matrix_Sint$beta2)*(params_matrix_Sint$alpha)
#compute total outflows
params_matrix_Sint$Total_L1=params_matrix_Sint$Time_L1^(-1)
params_matrix_Sint$Total_L2=params_matrix_Sint$Time_L2^(-1)
params_matrix_Sint$Total_I=params_matrix_Sint$Time_I^(-1)

params_matrix_Sint$epsilon=params_matrix_Sint$P_epsilon*params_matrix_Sint$Total_L1
params_matrix_Sint$epsilon0=params_matrix_Sint$epsilon*params_matrix_Sint$a
params_matrix_Sint$epsilon1=params_matrix_Sint$epsilon*params_matrix_Sint$b
params_matrix_Sint$epsilon2=params_matrix_Sint$epsilon*params_matrix_Sint$c
#compute kappa
params_matrix_Sint$kappa=params_matrix_Sint$Total_L1-params_matrix_Sint$epsilon+params_matrix_Sint$mu  
#compute nu
params_matrix_Sint$nu=params_matrix_Sint$P_nu*params_matrix_Sint$Total_L2
params_matrix_Sint$nu0=params_matrix_Sint$nu*params_matrix_Sint$a
params_matrix_Sint$nu1=params_matrix_Sint$nu*params_matrix_Sint$b
params_matrix_Sint$nu2=params_matrix_Sint$nu*params_matrix_Sint$c
#compute mui
params_matrix_Sint$mui0=params_matrix_Sint$P_mui0*params_matrix_Sint$Total_I
params_matrix_Sint$mui1=params_matrix_Sint$P_mui1*params_matrix_Sint$Total_I
params_matrix_Sint$mui2=params_matrix_Sint$P_mui2*params_matrix_Sint$Total_I
#compute h and j
params_matrix_Sint$h=params_matrix_Sint$P_h*params_matrix_Sint$Total_I
params_matrix_Sint$j=params_matrix_Sint$P_j*params_matrix_Sint$Total_I

#compute gamma
params_matrix_Sint$gamma0=params_matrix_Sint$Total_I-
  (params_matrix_Sint$mu+params_matrix_Sint$h+params_matrix_Sint$mui0)
params_matrix_Sint$gamma1=params_matrix_Sint$Total_I-
  (params_matrix_Sint$mu+params_matrix_Sint$j+params_matrix_Sint$mui1)
params_matrix_Sint$gamma2=params_matrix_Sint$Total_I-
  (params_matrix_Sint$mui2+params_matrix_Sint$mu)

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
  CDR_i_posibl=params$cdr_b+(1-params$cdr_b)*(params$q)#possible intervention size
  CDR_max=0.9#maximum possible CDR that the health system acan do
  CDR_i=min(CDR_i_posibl,CDR_max)#intervention CDR that can be used in model
  params_matrix_Sint$CDR_inter[i]=CDR_i#fill the parameter matrix
  
  cdr2_imax_possible=(CDR_i-(params$a+params$b)*params$cdr_b)/params$c#super-spreader CDR intervention 
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
  times=seq(0, 1000, by = 1)
  B_out <- as.data.frame(lsoda(initial_values, times, Baseline_model, params))
 Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
 incidence_B_out=(diff(B_out$inc)/Nq)*100000
  
 ChangeI=abs(incidence_B_out[length(incidence_B_out)-1]-
   incidence_B_out[length(incidence_B_out)-2])
 
 while(ChangeI>1.0e-6){
   initial_values=c(S=min(B_out$S),L1=max(B_out$L1),L2=max(B_out$L2),I0=max(B_out$I0),I1=max(B_out$I1),
                    I2=max(B_out$I2),inc=max(B_out$I0)+max(B_out$I1)+max(B_out$I2))
   times=times=seq(0, 1000, by = 1)
   B_out <- as.data.frame(lsoda(initial_values, times, Baseline_model, params))
   Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
  incidence_B_out=(diff(B_out$inc)/Nq)*100000
   ChangeI=abs(incidence_B_out[length(incidence_B_out)-1]-
     incidence_B_out[length(incidence_B_out)-2])
 }
 B_Incidence_time_n = incidence_B_out[length(incidence_B_out)-1] #the model reaches equilibrium at time around 500
 params_matrix_Sint$B_incidence[i]=B_Incidence_time_n
  
  #record compartment values at equilibrium
  Equi_initial_values=c(S=min(B_out$S), L1=max(B_out$L1), L2=max(B_out$L2), I0=max(B_out$I0),I1=max(B_out$I1),I2=max(B_out$I2),inc=max(B_out$I0)+max(B_out$I1)+max(B_out$I2))
  
  #RUN Intervention model, 
  Int_out <- as.data.frame(lsoda(Equi_initial_values, times, intervention_model, params))
  
  Nq=Int_out$S+Int_out$L1+Int_out$L2+Int_out$I0+Int_out$I1+Int_out$I2
  Int_incidence_out=(diff(Int_out$inc)/Nq)*100000
  
  #considering a 20 years impact of intervention
  Int_Incidence_time_n = Int_incidence_out[20]
  params_matrix_Sint$Int_incidence[i]=Int_Incidence_time_n
  #calculate the difference between baseline incidence and intervention incidence
  Change_inc=B_Incidence_time_n-Int_Incidence_time_n 
  params_matrix_Sint$Change_incidence[i]=Change_inc
  Relative_inc=((B_Incidence_time_n-Int_Incidence_time_n)/B_Incidence_time_n)*100
  params_matrix_Sint$Relative_change[i]=Relative_inc
}

View(params_matrix_Sint)
write.csv(x=params_matrix_Sint,file = '..//heterogeneity/params_matrix_Sint.csv')

#Plot Relative incidence % vs parameter values 

#CDR baseline

i_cdr_b=data.frame(cbind(params_matrix_Sint$cdr_b, params_matrix_Sint$Relative_change))
View(i_cdr_b)
i_cdr_b_g=melt(i_cdr_b,id='X1')
i_cdr_b_gp<-ggplot(i_cdr_b_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  
  xlab("baselien CDR")+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")+
ylim(0,100)
i_cdr_b_gp
#level of intervention (q)

i_q=data.frame(cbind(params_matrix_Sint$q, params_matrix_Sint$Relative_change))

i_q_g=melt(i_q,id='X1')
i_q_gp<-ggplot(i_q_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("ACF level (q)")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_q_gp
#level of intervention (d)

i_d=data.frame(cbind(params_matrix_Sint$d, params_matrix_Sint$Relative_change))

i_d_g=melt(i_d,id='X1')
i_d_gp<-ggplot(i_d_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("targeting level (d)")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_d_gp
#beta2

i_beta2=data.frame(cbind(params_matrix_Sint$beta2, params_matrix_Sint$Relative_change))

i_beta2_g=melt(i_beta2,id='X1')
i_beta2_gp<-ggplot(i_beta2_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("beta2")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_beta2_gp


#Epsilon
i_epsilon=data.frame(cbind(params_matrix_Sint$P_epsilon, params_matrix_Sint$Relative_change))

i_epsilon_g=melt(i_epsilon,id='X1')
i_epsilon_gp<-ggplot(i_epsilon_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("epsilon")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_epsilon_gp
#Nu
i_nu=data.frame(cbind(params_matrix_Sint$P_nu, params_matrix_Sint$Relative_change))

i_nu_g=melt(i_nu,id='X1')
i_nu_gp<-ggplot(i_nu_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("nu")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_nu_gp
#
#kappa
i_kappa=data.frame(cbind(params_matrix_Sint$P_kappa, params_matrix_Sint$Relative_change))

i_kappa_g=melt(i_kappa,id='X1')
i_kappa_gp<-ggplot(i_kappa_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("kappa")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_kappa_gp

#
#gamma2
i_gamma2=data.frame(cbind(params_matrix_Sint$P_gamma2, params_matrix_Sint$Relative_change))

i_gamma2_g=melt(i_gamma2,id='X1')
i_gamma2_gp<-ggplot(i_gamma2_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("i_gamma2")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_gamma2_gp
#
#gamma1
i_gamma1=data.frame(cbind(params_matrix_Sint$P_gamma1, params_matrix_Sint$Relative_change))

i_gamma1_g=melt(i_gamma1,id='X1')
i_gamma1_gp<-ggplot(i_gamma1_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("gamma1")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_gamma1_gp
#
#gamma0
i_gamma0=data.frame(cbind(params_matrix_Sint$P_gamma0, params_matrix_Sint$Relative_change))
i_gamma0_g=melt(i_gamma0,id='X1')
i_gamma0_gp<-ggplot(i_gamma0_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("gamma0")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_gamma0_gp
#
#mui2
i_mui2=data.frame(cbind(params_matrix_Sint$P_mui2, params_matrix_Sint$Relative_change))
i_mui2_g=melt(i_mui2,id='X1')
i_mui2_gp<-ggplot(i_mui2_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("mui2")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_mui2_gp
#
#mui1
i_mui1=data.frame(cbind(params_matrix_Sint$P_mui1, params_matrix_Sint$Relative_change))
i_mui1_g=melt(i_mui1,id='X1')
i_mui1_gp<-ggplot(i_mui1_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("mui1")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_mui1_gp
#mui0
i_mui0=data.frame(cbind(params_matrix_Sint$P_mui0, params_matrix_Sint$Relative_change))
i_mui0_g=melt(i_mui0,id='X1')
i_mui0_gp<-ggplot(i_mui0_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("mui0")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_mui0_gp

#h
i_h=data.frame(cbind(params_matrix_Sint$P_h, params_matrix_Sint$Relative_change))
i_h_g=melt(i_h,id='X1')
i_h_gp<-ggplot(i_h_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("h")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_h_gp

#j
i_j=data.frame(cbind(params_matrix_Sint$P_j, params_matrix_Sint$Relative_change))
i_j_g=melt(i_j,id='X1')
i_j_gp<-ggplot(i_j_g, aes(x=X1, y=value)) + 
  geom_point()+
  geom_smooth(method = 'lm',se=FALSE,col='darkred')+
  xlab("j")+
  ylim(0,100)+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
i_j_gp




i_GG<-ggarrange(i_cdr_b_gp,i_q_gp,i_d_gp,
                 i_beta2_gp,i_nu_gp, i_epsilon_gp,
               i_kappa_gp,i_gamma2_gp,i_gamma1_gp,
             i_gamma0_gp,i_mui2_gp,i_mui0_gp,
             i_mui1_gp,i_j_gp,i_h_gp,
              ncol=3, nrow=5, common.legend = TRUE,
              legend="bottom")

annotate_figure(i_GG,left = text_grob("Relative incidence %",rot = 90))
#


