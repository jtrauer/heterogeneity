#high and low value of parameters, mannualy save incidence at each run
#one way sensitivity of equi incidence to params value
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
require(lhs)
library(sensitivity)
#SENSITIVITY ANALYSIS
#Model Func
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
T_L1=1/4#time in L1
T_L2=20#time in L2
T_I=3#time in I
N=1


a=0.58
b=0.31
c=0.11
alpha=0.22
r=0.21

beta0<-0
beta2<-45
beta1<-beta2*alpha

P_epsilon=0.095
epsilon0<-(log(1-P_epsilon)/(-T_L1))*a
epsilon1<-(log(1-P_epsilon)/(-T_L1))*b
epsilon2<-(log(1-P_epsilon)/(-T_L1))*c

P_kappa=0.599
kappa<-(log(1-P_kappa)/(-T_L1))

P_nu=0.039
nu0<-(log(1-P_nu)/(-T_L2))*a
nu1<-(log(1-P_nu)/(-T_L2))*b
nu2<-(log(1-P_nu)/(-T_L2))*c

mu= 1/65

P_mui0=0.074
P_mui1=0.074
P_mui2=0.414
mui0=log(1-P_mui0)/(-T_I)
mui1=log(1-P_mui1)/(-T_I)
mui2=log(1-P_mui2)/(-T_I)

P_gamma0=0.583
P_gamma1=0.583
P_gamma2=0.343
gamma0<-(log(1-P_gamma0)/(-T_I))
gamma1<-(log(1-P_gamma1)/(-T_I))
gamma2<-(log(1-P_gamma2)/(-T_I))

#baseline CDR 0.7
cdr_b=0.7



#treatment success rate 
s0=0.8 
s1=0.8
s2=0.8
P_h=0.044
P_j=0.044


h<-(log(1-P_h)/(-T_I))
j<-(log(1-P_j)/(-T_I))
p1<-0
p2<-0



#compute baseline delta
delta0_b=(cdr_b*(gamma0+mui0+mu+h)/(1-cdr_b))*s0
delta1_b=(cdr_b*(gamma1+mui1+mu+j)/(1-cdr_b))*s1
delta2_b=(cdr_b*(gamma2+mui2+mu)/(1-cdr_b))*s2

parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, beta2=beta2,epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                 kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                 mui0=mui0,mui1=mui1,mui2=mui2,delta0_b=delta0_b,delta1_b=delta1_b,delta2_b=delta2_b,r=r,
                 h=h,j=j,p1=p1,p2=p2) 
##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6*a    #active TB hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # active TB Normal spreaders
F=1e-6*c # active TB Super-spreader 10% of TB patients are super-spreaders
G=0 # diagnosed 


#View(output_matrix_equi)

initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)

times=seq(0, 2000, by = 1)
B_out <- as.data.frame(lsoda(initial_values, times, Baseline_model, parameter_list))

#Record Baseline equilibrium incidence
Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2

incidence_B_out=(diff(B_out$inc)/Nq)*100000

plot(incidence_B_out)
#max(incidence_B_out)
B_Incidence_time_n = incidence_B_out[1999] #the model reaches equilibrium at time around 500
print(B_Incidence_time_n)


















