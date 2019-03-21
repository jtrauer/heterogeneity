
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
require(lhs)
library(sensitivity)
#SENSITIVITY ANALYSIS
#Model Func
Baseline_model=function(current_timepoint, state_values, parameters)
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


#Latin hypercube sampling

z <- 100 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(z,32) #simulate h= number of simulations, 35=number of parameters
#To map these points in the unit cube to our parameters, we need minimum and maximum values for each.

a.min=0.58
a.max=0.58
b.min=0.31
b.max=0.31
c.min=0.11
c.max=0.11
alpha.min=0.22
alpha.max=0.22

mu.min= 0.0133
mu.max= 0.0182
mui0.min=0.0167
mui0.max=0.032
mui1.min=0.0167
mui1.max=0.032
mui2.min=0.109
mui2.max=0.262
r.min=0.21
r.max=0.21

beta0.min<-0
beta0.max<-0
beta2.min<-30
beta2.max<-60
beta1.min<-6.6
beta1.max<-13.2
epsilon0.min<- 0.1779498
epsilon0.max<-0.3177646
epsilon1.min<-0.0951111
epsilon1.max<-0.1698397
epsilon2.min<-0.0337491
epsilon2.max<- 0.0602657
kappa.min<-3.0104625
kappa.max<-5.1135

gamma0.min<-0.265
gamma0.max<-0.301
gamma1.min<-0.265
gamma1.max<-0.301
gamma2.min<-0.056
gamma2.max<-0.209

nu0.min<-0.0005278
nu0.max<-0.002330295
nu1.min<-0.0002821
nu1.max<-0.001245502
nu2.min<-0.0001001
nu2.max<-0.0004419525
#baseline CDR 0.7
cdr0_b.min=0.7
cdr0_b.max=0.7
cdr1_b.min=0.7
cdr1_b.max=0.7
cdr2_b.min=0.7
cdr2_b.max=0.7

#treatment success rate 
s0.min=0.8 
s0.max=0.8
s1.min=0.8
s1.max=0.8
s2.min=0.8
s2.max=0.8

h.min<-0.0
h.max<-0.02
j.min<-0.0
j.max<-0.02
p1.min<-0
p1.max<-0
p2.min<-0
p2.max<-0

#Now we can generate a “parameter set” by rescaling our simulated latin hypercube sample
params.set_o <- cbind(
  a = lhs[,1]*(a.max-a.min)+a.min,
  b = lhs[,2]*(b.max-b.min)+b.min,
  c = lhs[,3]*(c.max-c.min)+c.min,
  alpha = lhs[,4]*(alpha.max-alpha.min)+alpha.min,
  beta0 = lhs[,5]*(beta0.max-beta0.min)+beta0.min,
  beta1 = lhs[,6]*(beta1.max-beta1.min)+beta1.min,
  beta2 = lhs[,7]*(beta2.max-beta2.min)+beta2.min,
  epsilon0 =lhs[,8]*(epsilon0.max-epsilon0.min)+epsilon0.min,
  epsilon1 =lhs[,9]*(epsilon1.max-epsilon1.min)+epsilon1.min,
  epsilon2 =lhs[,10]*(epsilon2.max-epsilon2.min)+epsilon2.min,
  kappa = lhs[,11]*(kappa.max-kappa.min)+kappa.min,
  gamma0 = lhs[,12]*(gamma0.max-gamma0.min)+gamma0.min,
  gamma1 = lhs[,13]*(gamma1.max-gamma1.min)+gamma1.min,
  gamma2 = lhs[,14]*(gamma2.max-gamma2.min)+gamma2.min,
  nu0 = lhs[,15]*(nu0.max-nu0.min)+nu0.min,
  nu1 = lhs[,16]*(nu1.max-nu1.min)+nu1.min,
  nu2 = lhs[,17]*(nu2.max-nu2.min)+nu2.min,
  mu = lhs[,18]*(mu.max-mu.min)+mu.min,
  mui0 = lhs[,19]*(mui0.max-mui0.min)+mui0.min,
  mui1 = lhs[,20]*(mui1.max-mui1.min)+mui1.min,
  mui2 = lhs[,21]*(mui2.max-mui2.min)+mui2.min,
  r = lhs[,22]*(r.max-r.min)+r.min,
  s0 = lhs[,23]*(s0.max-s0.min)+s0.min,
  s1 = lhs[,24]*(s1.max-s1.min)+s1.min,
  s2 = lhs[,25]*(s2.max-s2.min)+s2.min,
  cdr0_b= lhs[,26]*(cdr0_b.max-cdr0_b.min)+cdr0_b.min,
  cdr1_b= lhs[,27]*(cdr1_b.max-cdr1_b.min)+cdr1_b.min,
  cdr2_b= lhs[,28]*(cdr2_b.max-cdr2_b.min)+cdr2_b.min,

  h = lhs[,29]*(h.max-h.min)+h.min,
  j = lhs[,30]*(j.max-j.min)+j.min,
  p1 = lhs[,31]*(p1.max-p1.min)+p1.min,
  p2 = lhs[,32]*(p2.max-p2.min)+p2.min)

#View(params.set_o)
#class(params.set_o)

#creat matrix to save whole info
params_matrix = data.frame(params.set_o)
#add colums for delta base line
delta0_b = data.frame('delta0_b'=rep(NA,z))
delta1_b = data.frame('delta1_b'=rep(NA,z))
delta2_b = data.frame('delta2_b'=rep(NA,z))

params_matrix_equi=cbind(params_matrix,delta0_b,delta1_b,delta2_b)
#View(params_matrix_equi)
#creat columns for incidences
Equi_incidence = data.frame('Equi_incidence'=rep(NA,z))
output_matrix_equi = cbind(params_matrix_equi, Equi_incidence) #add incidence column

View(output_matrix_equi)

#compute baseline delta
output_matrix_equi$delta0_b=(output_matrix_equi$cdr0_b*(output_matrix_equi$gamma0+output_matrix_equi$mui0+output_matrix_equi$mu+output_matrix_equi$h)/(1-output_matrix_equi$cdr0_b))*output_matrix_equi$s0
output_matrix_equi$delta1_b=(output_matrix_equi$cdr1_b*(output_matrix_equi$gamma1+output_matrix_equi$mui1+output_matrix_equi$mu+output_matrix_equi$j)/(1-output_matrix_equi$cdr1_b))*output_matrix_equi$s1
output_matrix_equi$delta2_b=(output_matrix_equi$cdr2_b*(output_matrix_equi$gamma2+output_matrix_equi$mui2+output_matrix_equi$mu)/(1-output_matrix_equi$cdr2_b))*output_matrix_equi$s2
View(output_matrix_equi)
##Initial values for sub population: 
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

View(output_matrix_equi) #now we have incidence and each parameters in one matrix

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







