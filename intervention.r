
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
require(lhs)
library(sensitivity)


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

#Latin hypercube sampling

z <- 100 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(z,28) #simulate h= number of simulations, 35=number of parameters
#To map these points in the unit cube to our parameters, we need minimum and maximum values for each.
#proportions for flows
a.min=0.58
a.max=0.58
b.min=0.31
b.max=0.31
c.min=0.11
c.max=0.11

#natural mortality
mu.min= 0.0133
mu.max= 0.0182
#TB related mortality
mui0.min=0.0167
mui0.max=0.032
mui1.min=0.0167
mui1.max=0.032
mui2.min=0.109
mui2.max=0.262
# reduced force re-infection 
r.min=0.21
r.max=0.21
#transmission parametters
beta0.min<-0
beta0.max<-0
beta2.min<-30
beta2.max<-60
beta1.min<-6.6
beta1.max<-13.2
#early progressions
epsilon0.min<- 0.1779498
epsilon0.max<-0.3177646
epsilon1.min<-0.0951111
epsilon1.max<-0.1698397
epsilon2.min<-0.0337491
epsilon2.max<- 0.0602657
#progression to late latency
kappa.min<-3.0104625
kappa.max<-5.1135
#self recovery
gamma0.min<-0.265
gamma0.max<-0.301
gamma1.min<-0.265
gamma1.max<-0.301
gamma2.min<-0.056
gamma2.max<-0.209
#late reactivation
nu0.min<-0.0005278
nu0.max<-0.002330295
nu1.min<-0.0002821
nu1.max<-00.001245502
nu2.min<-0.0001001
nu2.max<-0.0004419525
#conversion from non-spreaders to low-spreaders
h.min<-0.0
h.max<-0.02
#conversion from low-spreaders to super-spreaders
j.min<-0.0
j.max<-0.02
#IPT for contacts of low-spreader
p1.min<-0
p1.max<-0
#IPT for contacts of super-spreader
p2.min<-0
p2.max<-0
#baseline CDR
CDR_b.min<-0.2
CDR_b.max<-0.7
#proportion for targeting super-spreaders
d.min<-0
d.max<-1
#proportion of active case finding missed cases
q.min<-0
q.max<-1

#Now we can generate a “parameter set” by rescaling our simulated latin hypercube sample
params.set_uj <- cbind(
  
  a = lhs[,1]*(a.max-a.min)+a.min,
  b = lhs[,2]*(b.max-b.min)+b.min,
  c = lhs[,3]*(c.max-c.min)+c.min,
  
  beta0 = lhs[,4]*(beta0.max-beta0.min)+beta0.min,
  beta1 = lhs[,5]*(beta1.max-beta1.min)+beta1.min,
  beta2 = lhs[,6]*(beta2.max-beta2.min)+beta2.min,
  
  epsilon0 =lhs[,7]*(epsilon0.max-epsilon0.min)+epsilon0.min,
  epsilon1 =lhs[,8]*(epsilon1.max-epsilon1.min)+epsilon1.min,
  epsilon2 =lhs[,9]*(epsilon2.max-epsilon2.min)+epsilon2.min,
  
  kappa = lhs[,10]*(kappa.max-kappa.min)+kappa.min,
  
  gamma0 = lhs[,11]*(gamma0.max-gamma0.min)+gamma0.min,
  gamma1 = lhs[,12]*(gamma1.max-gamma1.min)+gamma1.min,
  gamma2 = lhs[,13]*(gamma2.max-gamma2.min)+gamma2.min,
  
  nu0 = lhs[,14]*(nu0.max-nu0.min)+nu0.min,
  nu1 = lhs[,15]*(nu1.max-nu1.min)+nu1.min,
  nu2 = lhs[,16]*(nu2.max-nu2.min)+nu2.min,
  
  mu = lhs[,17]*(mu.max-mu.min)+mu.min,
  
  mui0 = lhs[,18]*(mui0.max-mui0.min)+mui0.min,
  mui1 = lhs[,19]*(mui1.max-mui1.min)+mui1.min,
  mui2 = lhs[,20]*(mui2.max-mui2.min)+mui2.min,
  
  r = lhs[,21]*(r.max-r.min)+r.min,
  h = lhs[,22]*(h.max-h.min)+h.min,
  j = lhs[,23]*(j.max-j.min)+j.min,
  p1 = lhs[,24]*(p1.max-p1.min)+p1.min,
  p2 = lhs[,25]*(p2.max-p2.min)+p2.min,
  CDR_b = lhs[,26]*(CDR_b.max-CDR_b.min)+CDR_b.min,
  d=lhs[,27]*(d.max-d.min)+d.min,
  q=lhs[,28]*(q.max-q.min)+q.min)
#View(params.set_uj)
##Initial values for sub population:

N=1
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6*a    #active TB hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # active TB Normal spreaders
F=1e-6*c # active TB Super-spreader 10% of TB patients are super-spreaders
G=0 # diagnosed 

#creat additional columns 
CDR_inter=data.frame('CDR_inter'=rep(NA,z))#intervention level CDR
B_incidence = data.frame('B_incidence'=rep(NA,z))#baseline equilibrium incidence
Int_incidence = data.frame('Int_incidence'=rep(NA,z))#10years intervention incidence
Change_incidence= data.frame('Change_incidence'=rep(NA,z))#absolute incidence difference
Relative_change=data.frame('Relative_change'=rep(NA,z))#relative incidence difference

output_matrix_j = cbind(params.set_uj, CDR_inter,B_incidence,Int_incidence,Change_incidence,Relative_change) #add incidence column
#View(output_matrix_j)



for(i in 1:z){
  #run baseline 
  params <- as.list(c(output_matrix_j[i,]))
  cdr0_b=params$CDR_b#baseline CDR
  cdr1_b=params$CDR_b
  cdr2_b=params$CDR_b
  
  CDR_i_posibl=params$CDR_b+(1-params$CDR_b)*(params$q)#possible intervention size
  CDR_max=0.9#maximum possible CDR that the health system acan do
  CDR_i=min(CDR_i_posibl,CDR_max)#intervention CDR that can be used in model
  output_matrix_j$CDR_inter[i]=CDR_i#fill the parameter matrix
  
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

#plot each parameter with relative change in incidence

#q
q=cbind(output_matrix_j$q,output_matrix_j$Relative_change)
write.csv(x=q,file = '..//Enchik.com/q.csv')
gg_q=read.csv(file.choose())
gg_q=melt(gg_q,id='q')
gp_q<-ggplot(gg_q, aes(x=q, y=value)) + 
  geom_point()+
  geom_smooth()+
  ylim(0,100)+
  xlab("proportion of missed cases detected")+
  ylab("Relative incidence reduction %")
gp_q

#d
d=cbind(output_matrix_j$d,output_matrix_j$Relative_change)
write.csv(x=q,file = '..//Enchik.com/d.csv')
gg_d=read.csv(file.choose())
gg_d=melt(gg_d,id='d')

gp_d<-ggplot(gg_d, aes(x=d, y=value)) + 
  geom_point()+
  geom_smooth()+
  ylim(0,100)+
  xlab("proportion of missed cases detected")+
  ylab("Relative incidence reduction %")
gp_d
#beta1
beta1=cbind(output_matrix_j$beta1,output_matrix_j$Relative_change)
View(beta1)
write.csv(x=beta1,file='..//Enchik.com/beta1.csv')
gg_beta1=read.csv(file.choose())
gg_beta1<-melt(gg_beta1,id='beta1')


gp_beta1<-ggplot(gg_beta1, aes(x=beta1, y=value)) + 
  geom_point()+
  geom_smooth()+
  ylim(0,100)+
  xlab("beta1")+
  ylab("Relative incidence reduction %")
gp_beta1

#beta2
beta2=cbind(output_matrix_j$beta2,output_matrix_j$Relative_change)
View(beta2)
write.csv(x=beta2,file='..//Enchik.com/beta2.csv')
gg_beta2=read.csv(file.choose())
gg_beta2<-melt(gg_beta2,id='beta2')


gp_beta2<-ggplot(gg_beta2, aes(x=beta2, y=value)) + 
  geom_point()+
  geom_smooth()+
  ylim(0,100)+
  xlab("beta2")+
  ylab("Relative incidence reduction %")
gp_beta2
ggarrange(gp_beta1,gp_beta2,gp_q,gp_d,
          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")





