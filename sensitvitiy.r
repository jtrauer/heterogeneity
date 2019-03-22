
library(deSolve)
library(reshape2)
library(ggplot2)
#library(ggpubr)
require(lhs)


library(sensitivity)

source("function_tool_kit.R")

# SENSITIVITY ANALYSIS

# Model Func
Baseline_model <- function(current_timepoint, state_values, parameters)
  {inst
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
# Time in each compartment
T_L1=1/4#time in L1 in years
T_L2=20#time in L2 in years
T_I=3#Time in I in years


#Latin hypercube sampling
z <- 10 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(z,24) #simulate h= number of simulations, 35=number of parameters
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
                           #here
                           r = list(min = 0.21, max = 0.21),
                           beta0 = list(min = 0, max = 0),
                           beta2 = list(min = 30, max = 60),
                           P_epsilon = list(min = 0.074, max = 0.128),
                           P_kappa = list(min = 0.54, max = 0.721),
                           
                           P_gamma0 = list(min = 0.548, max = 0.595),
                           P_gamma1 = list(min = 0.548, max = 0.548),
                           P_gamma2 = list(min = 0.155, max = 0.466),
                           P_nu = list(min = 0.018, max = 0.077),
                          
                           cdr_b = list(min = 0.5, max = 0.8),#baseline CDR
                           
                           s = list(min = 0.8, max = 0.8),#treatment success
                           
                           P_h = list(min = 0.0, max = 0.058),
                           P_j = list(min = 0.0, max = 0.058),
                           p1 = list(min = 0, max = 0),
                           p2 = list(min = 0, max = 0))


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


#creat matrix to save whole info
params_matrix = data.frame(params.set_o)
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
#add colums for delta base line
delta0_b = data.frame('delta0_b'=rep(NA,z))
delta1_b = data.frame('delta1_b'=rep(NA,z))
delta2_b = data.frame('delta2_b'=rep(NA,z))

params_matrix_equi=cbind(params_matrix,beta1,epsilon0,epsilon1,epsilon2,
                         mui0,mui1,mui2,kappa, nu0,nu1,nu2,gamma0,gamma1,gamma2,h,j,
                         delta0_b,delta1_b,delta2_b)
View(params_matrix_equi)
#creat columns for incidences
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

#######################

write.csv(x=output_matrix_equi,file='..//Enchik.com/output_matrix_equi.csv')

#creat matrix to save whole info

#which parameter affects more,inspecting the partial correlation
#Partial rank correlations can be computed using the pcc function in the R package sensitivity 
#install.packages('sensitivity')

bonferroni.alpha <- 0.05/35

prcc <- pcc(output_matrix_equi[,1:35], output_matrix_equi[,36], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')

#We can view a table of the resulting partial correlation coefficients. if none of the (penalized)
#confidence intervals contains zero, we conclude that all are significant and produce a plot showing their
#relative magnitudes.

load('prcc.Rdata')
summary <- print(prcc)
Corl=data.frame(summary) 
View(Corl)

write.csv(x=Corl,file='..//Enchik.com//corl.csv')
#plote the partial corelation coeeficient

par(mar=c(9,4,4,2)+0.1)
plot(Corl$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',
     axes=FALSE)
axis(2)
axis(1, at=seq(1:35), labels=row.names(Corl), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:35) lines(c(i,i),c(Corl[i,4], Corl[i,5]))
abline(h=0)
#tornado plot


#We can plot each simulated value as a point
par(mfrow=c(4,3))
#plot(output_matrix$beta0, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
# xlab='beta0',
#ylab='Incidence',main = 'PRCC= 0.017')

plot(output_matrix_equi$nu2, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.57, col='blue',
     xlab=expression(nu2),
     ylab='Equilibrium incidence',main = 'PRCC= 0.8456')
plot(output_matrix_equi$beta2, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(beta2),
     ylab='Equilibrium incidence',main = 'PRCC= 0.8332')
plot(output_matrix_equi$nu1, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(nu1),
     ylab='Equilibrium incidence',main = 'PRCC= 0.8262')
plot(output_matrix_equi$kappa, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(kappa),
     ylab='Equilibrium incidence',main = 'PRCC= -0.8178')

plot(output_matrix_equi$nu0, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(nu0),
     ylab='Equilibrium incidence',main = 'PRCC= 0.7070')

plot(output_matrix_equi$beta1, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab='beta1',
     ylab='Equilibrium incidence',main = 'PRCC= 0.6678')

plot(output_matrix_equi$epsilon2, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(epsilon2),
     ylab='Equilibrium incidence',main = 'PRCC= 0.016')

plot(output_matrix_equi$epsilon1, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=.8,pch=19, cex=0.7, col='blue',
     xlab=expression(epsilon1),
     ylab='Equilibrium incidence',main = 'PRCC= 0.5039')

plot(output_matrix_equi$epsilon0, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='blue',
     xlab=expression(epsilon0),
     ylab='Equilibrium incidence',main = 'PRCC= 0.32084')

plot(output_matrix_equi$mu, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='blue',
     xlab=expression(mu),
     ylab='Equilibrium incidence',main = 'PRCC= -0.1620')

plot(output_matrix_equi$gamma2, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='blue',
     xlab=expression(gamma2),
     ylab='Equilibrium incidence',main = 'PRCC= -0.1367')

plot(output_matrix_equi$mui2, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='blue',
     xlab=expression(mui2),
     ylab='Equilibrium incidence',main = 'PRCC= -0.1292')

plot(output_matrix_equi$mui0, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(mui0),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0718')

plot(output_matrix_equi$mui1, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(mui1),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0612')
plot(output_matrix_equi$gamma0, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(gamma0),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0393')

plot(output_matrix_equi$gamma1, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(gamma1),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0269')
plot(output_matrix_equi$h, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(h),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0144')

plot(output_matrix_equi$j, output_matrix_equi$Equi_incidence, type = 'p' ,lwd=0.8,pch=19, cex=0.7, col='coral4',
     xlab=expression(j),
     ylab='Equilibrium incidence',main = 'PRCC= -0.0037')







