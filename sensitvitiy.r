
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(lhs)
library(sensitivity)

source("function_tool_kit.R")

# SENSITIVITY ANALYSIS

# Model Func
Baseline_model <- function(current_timepoint, state_values, parameters)
  {#inst #what is "inst" James?Error in func(time, state, parms, ...) :
      #object 'inst' not found
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

# Time in each compartment
T_L1=1/4#time in L1 in years
T_L2=20#time in L2 in years
T_I=3#Time in I in years

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
                           p2 = list(min = 0, max = 0))

# Latin hypercube sampling
z <- 1000 # choose number of points to simulate
set.seed(6242015) # random number generator
# To map these points in the unit cube to our parameters, we need minimum and maximum values for each.

lhs <- maximinLHS(z, length(param_value_limits)) # simulate h= number of simulations

#Now we can generate a ?parameter set? by rescaling our simulated latin hypercube sample
params.set_o <- cbind(
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
  N= adjust_lhs_to_range(lhs[, 24], "N", param_value_limits))

View(params.set_o)

# create matrix to save whole info
params_matrix = data.frame(params.set_o)
# add colums for parameters 
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

# add colums for delta base line
delta0_b = data.frame('delta0_b'=rep(NA,z))
delta1_b = data.frame('delta1_b'=rep(NA,z))
delta2_b = data.frame('delta2_b'=rep(NA,z))

params_matrix_equi=cbind(params_matrix,beta1,epsilon0,epsilon1,epsilon2,
                         mui0,mui1,mui2,kappa, nu0,nu1,nu2,gamma0,gamma1,gamma2,h,j,
                         delta0_b,delta1_b,delta2_b)
View(params_matrix_equi)

# create columns for incidences
Equi_incidence = data.frame('Equi_incidence'=rep(NA,z))
output_matrix_equi = cbind(params_matrix_equi, Equi_incidence) #add incidence column

#View(output_matrix_equi)
#compute beta1
output_matrix_equi$beta1=(output_matrix_equi$beta2)*(output_matrix_equi$alpha)
#compute epsilon
output_matrix_equi$epsilon0=log(1-(output_matrix_equi$P_epsilon*output_matrix_equi$a))/(-T_L1)
output_matrix_equi$epsilon1=log(1-(output_matrix_equi$P_epsilon*output_matrix_equi$b))/(-T_L1)
output_matrix_equi$epsilon2=log(1-(output_matrix_equi$P_epsilon*output_matrix_equi$c))/(-T_L1)
#compute kappa
output_matrix_equi$kappa=log(1-(output_matrix_equi$P_kappa))/(-T_L1)
#compute nu
output_matrix_equi$nu0=log(1-(output_matrix_equi$P_nu*output_matrix_equi$a))/(-T_L2)
output_matrix_equi$nu1=log(1-(output_matrix_equi$P_nu*output_matrix_equi$b))/(-T_L2)
output_matrix_equi$nu2=log(1-(output_matrix_equi$P_nu*output_matrix_equi$c))/(-T_L2)
#compute gamma
output_matrix_equi$gamma0=log(1-(output_matrix_equi$P_gamma0))/(-T_I)
output_matrix_equi$gamma1=log(1-(output_matrix_equi$P_gamma1))/(-T_I)
output_matrix_equi$gamma2=log(1-(output_matrix_equi$P_gamma2))/(-T_I)
#compute mui
output_matrix_equi$mui0=log(1-(output_matrix_equi$P_mui0))/(-T_I)
output_matrix_equi$mui1=log(1-(output_matrix_equi$P_mui1))/(-T_I)
output_matrix_equi$mui2=log(1-(output_matrix_equi$P_mui2))/(-T_I)
#compute h and j
output_matrix_equi$h=log(1-(output_matrix_equi$P_h))/(-T_I)
output_matrix_equi$j=log(1-(output_matrix_equi$P_j))/(-T_I)
#compute baseline delta
output_matrix_equi$delta0_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b, 
                      output_matrix_equi$gamma0 + output_matrix_equi$mui0 +
                        output_matrix_equi$mu + output_matrix_equi$h, 
                      output_matrix_equi$s)
output_matrix_equi$delta1_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b,
                      output_matrix_equi$gamma1 + output_matrix_equi$mui1 + 
                        output_matrix_equi$mu + output_matrix_equi$j, 
                      output_matrix_equi$s)
output_matrix_equi$delta2_b <-
  find_delta_from_cdr(output_matrix_equi$cdr_b,
                      output_matrix_equi$gamma2 + output_matrix_equi$mui2 + 
                        output_matrix_equi$mu,
                      output_matrix_equi$s)

View(output_matrix_equi)
## Initial values for sub population: 
a=0.58 #proportion for I0
b=0.31 #proportion for I1
c=0.11 #proportion for I2
A=1 #Fully susceptible hosts
B=0 #Early latent hosts
C=0 #Late latent hosts
D=1e-6*a #active TB hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # active TB Normal spreaders
F=1e-6*c # active TB Super-spreader 10% of TB patients are super-spreaders
G=0 # diagnosed 

#View(output_matrix_equi)
for(i in 1:z){
  #run baseline 
  
  initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
  params <- as.list(c(output_matrix_equi[i,]))
  times=seq(0, 5000, by = 1)
  B_out <- as.data.frame(lsoda(initial_values, times, Baseline_model, params))
  
  #Record Baseline equilibrium incidence
  Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
  
  incidence_B_out=(diff(B_out$inc)/Nq)*100000
  
  plot(incidence_B_out)
  #max(incidence_B_out)
  B_Incidence_time_n = incidence_B_out[4999] #the model reaches equilibrium at time around 500
  output_matrix_equi$Equi_incidence[i] = B_Incidence_time_n
  
} 

#View(output_matrix_equi) #now we have incidence and each parameters in one matrix

#save output as csv

write.csv(x=output_matrix_equi,file='..//heterogeneity/output/output_matrix_equi.csv')

#test significance of correlation, selcet parameters from csv manually 
Test_corelation<-read.csv(file.choose())#open csv with selected params and incidence
bonferroni.alpha <- 0.05/14 #14 parametrs

prcc <- pcc(Test_corelation[,1:14], Test_corelation[,15], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
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
for(i in 1:14) lines(c(i,i),c(Corl[i,4], Corl[i,5]))
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
  
  xlab(expression(nu))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_nu_gp
#kappa
e_kappa=data.frame(cbind(output_matrix_equi$P_kappa, output_matrix_equi$Equi_incidence))
e_kappa_g=melt(e_kappa,id='X1')
e_kappa_gp<-ggplot(e_kappa_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  xlab(expression(kappa))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_kappa_gp 

#beta2

e_beta2=data.frame(cbind(output_matrix_equi$beta2,output_matrix_equi$Equi_incidence))
e_beta2_g=melt(e_beta2,id='X1')
e_beta2_gp<-ggplot(e_beta2_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(beta))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_beta2_gp
#epsiolon
e_epsilon=data.frame(cbind(output_matrix_equi$P_epsilon,output_matrix_equi$Equi_incidence))
e_epsilon_g=melt(e_epsilon,id='X1')
e_epsilon_gp<-ggplot(e_epsilon_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(epsilon))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_epsilon_gp

#mu
e_mu=data.frame(cbind(output_matrix_equi$mu,output_matrix_equi$Equi_incidence))

e_mu_g=melt(e_mu,id='X1')
e_mu_gp<-ggplot(e_mu_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(mu))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mu_gp

#gamma2
e_gamma2=data.frame(cbind(output_matrix_equi$P_gamma2,output_matrix_equi$Equi_incidence))
e_gamma2_g=melt(e_gamma2,id='X1')
e_gamma2_gp<-ggplot(e_gamma2_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(gamma[2]))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_gamma2_gp
#gamma1
e_gamma1=data.frame(cbind(output_matrix_equi$P_gamma1,output_matrix_equi$Equi_incidence))

e_gamma1_g=melt(e_gamma1,id='X1')
e_gamma1_gp<-ggplot(e_gamma1_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(gamma[1]))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_gamma1_gp

#gamma0
e_gamma0=data.frame(cbind(output_matrix_equi$P_gamma0,output_matrix_equi$Equi_incidence))
e_gamma0_g=melt(e_gamma0,id='X1')
e_gamma0_gp<-ggplot(e_gamma0_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(gamma[0]))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_gamma0_gp
#mui2
e_mui2=data.frame(cbind(output_matrix_equi$P_mui2,output_matrix_equi$Equi_incidence))
e_mui2_g=melt(e_mui2,id='X1')
e_mui2_gp<-ggplot(e_mui2_g, aes(x=X1, y=value)) + 
  geom_point(color="coral4")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(mu[i2]))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mui2_gp

#mui1
e_mui1=data.frame(cbind(output_matrix_equi$P_mui1,output_matrix_equi$Equi_incidence))

e_mui1_g=melt(e_mui1,id='X1')
e_mui1_gp<-ggplot(e_mui1_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(mu[i1]))+
  #theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_mui1_gp
#mui0
e_mui0=data.frame(cbind(output_matrix_equi$P_mui0,output_matrix_equi$Equi_incidence))
e_mui0_g=melt(e_mui0,id='X1')
e_mui0_gp<-ggplot(e_mui0_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  
  xlab(expression(mu[i0]))+
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

#j
e_j=data.frame(cbind(output_matrix_equi$P_j,output_matrix_equi$Equi_incidence))
e_j_g=melt(e_j,id='X1')
e_j_gp<-ggplot(e_j_g, aes(x=X1, y=value)) + 
  geom_point(color="chocolate1")+
  #geom_smooth(method = lm,se=FALSE,col='darkred')+
  xlab("j")+
 # theme(plot.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))+
  ylab("")
e_j_gp

GG<-ggarrange(e_cdr_b_gp,e_nu_gp,e_beta2_gp,e_kappa_gp,
              e_epsilon_gp,e_gamma2_gp,e_mui2_gp,e_mu_gp,
              e_mui0_gp,e_mui1_gp,e_gamma0_gp,e_gamma1_gp,
              e_j_gp,e_h_gp,
              ncol=4, nrow=4, common.legend = TRUE,
              legend="bottom")

annotate_figure(GG,left = text_grob("Equilibrium incidence per 100,000popn",rot = 90))
#











