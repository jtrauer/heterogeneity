#Scenarios
library(deSolve)
library(reshape2)
library(ggplot2)
library(ggpubr)
threeIIC_model=function(current_timepoint, state_values, parameters)
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
      dS=(mu*N+mui0*I0+mui1*I1+mui2*I2)+(delta0*I0+delta1*I1+delta2*I2)+((r*(p1*beta1*I1+p2*beta2*I2)/N)*L2)-((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S-mu*S
      dL1=((beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*S+(r*(beta0*I0+(1-p1)*beta1*I1+(1-p2)*beta2*I2)/N)*L2-(epsilon0+epsilon1+epsilon2+mu+kappa)*L1
      dL2=kappa*L1+(gamma0*I0+gamma1*I1+gamma2*I2)-(r*(beta0*I0+beta1*I1+beta2*I2)/N)*L2-(nu0+nu1+nu2+mu)*L2
      dI0=(epsilon0*L1+nu0*L2)-(mui0+mu+gamma0+delta0+h)*I0
      dI1=(epsilon1*L1+nu1*L2+h*I0)-(mui1+mu+gamma1+delta1+j)*I1
      dI2=(epsilon2*L1+nu2*L2+j*I1)-(mui2+mu+gamma2+delta2)*I2
      
      dinc=(epsilon0*L1+nu0*L2)+(epsilon1*L1+nu1*L2)+(epsilon2*L1+nu2*L2)
      
      #combine results
      results=c(dS, dL1, dL2, dI0,dI1,dI2, dinc)
      list(results)
    }
  )
}
# Set parameters 
#proportions from Vic contact tracing data
N=1

a=0.58
b=0.31
c=0.11
alpha=0.22
r=0.21
Time_L1=1/4#time in L1 
Time_L2=20#time in L2
Time_I=3#time in I
mu= 1/65#

beta0<-0
beta2<-45#
P_epsilon=0.1004#
P_nu=0.039#
P_mui0=0.077#
P_mui1= 0.077#
P_mui2=0.534#
P_h=0.045#
P_j=0.045#
beta1<-beta2*alpha
epsilon=P_epsilon*(1/Time_L1)
epsilon0=epsilon*a
epsilon1=epsilon*b
epsilon2=epsilon*c
kappa=(1/Time_L1)-(epsilon+mu) 
nu=P_nu*(1/Time_L2)
nu0=nu*a
nu1=nu*b
nu2=nu*c

mui0=P_mui0*(1/Time_I)
mui1=P_mui1*(1/Time_I)
mui2=P_mui2*(1/Time_I)

h=P_h*(1/Time_I)
j=P_j*(1/Time_I)

gamma0=(1/Time_I)-(mu+h+mui0)
gamma1=(1/Time_I)-(mu+j+mui1)
gamma2=(1/Time_I)-( mu+mui2)
#baeseline scenario
CDR_b=0.5
#Baseline intervention
#d=0
#q=0
#mass intervention
#d=0
#q=0.1
#targeted intervention
d=0.8
q=0.1

CDR_i_posibl=CDR_b+(1-CDR_b)*(q)#possible intervention size
CDR_max=0.9#maximum possible CDR that the health system acan do
CDR_i=min(CDR_i_posibl,CDR_max)#intervention CDR that can be used in model

cdr2_imax_possible=(CDR_i-(a+b)*CDR_b)/c#super-spreader CDR intervention 
cdr2_imax=min(cdr2_imax_possible,CDR_max)#the limit that SS CDR can take

cdr2_i=CDR_i+((cdr2_imax-CDR_i)*d)# Super-Spreaders CDR to be used in the model
cdr1_i=(CDR_i-(c*cdr2_i))/(a+b)#low-spreaders CDR to be used in the model
cdr0_i=cdr1_i#non-spreaders CDR




#scenario 1 a 10% increase in total CDR=.77 untargeted
#cdr0=0.77
#cdr1=0.77
#cdr2=0.77
#scenario 2, extent of targeting
#cdr0=0.77
#cdr1=0.77
#cdr2=0.77

# code for all CDRS channgins as cdr2 changes and this change deltas
delta0=(cdr0_i*(gamma0+mui0+mu+h)/(1-cdr0_i))*s0
delta1=(cdr1_i*(gamma1+mui1+mu+j)/(1-cdr1_i))*s1
delta2=(cdr2_i*(gamma2+mui2+mu)/(1-cdr2_i))*s2
##Initial values for sub population: 
A=1 #Fully susceptible hosts
B=0        #Early latent hosts
C=0        #Late latent hosts
D=1e-6*a    #active TB hosts proportion of extrapulmonary 15% of all infectious 
E=1e-6*b # active TB Normal spreaders
F=1e-6*c # active TB Super-spreader 10% of TB patients are super-spreaders
G=0 # diagnosed 

parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, beta2=beta2,epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                 kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                 mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,
                 h=h,j=j,p1=p1,p2=p2)  



#Initial state values for differential equasions
#initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
#at Equilibrium
#Scenario1
#initial_values=c(S=  0.5914698, L1=0.002143275, L2=  0.4052267, I0= 0.0006635413,I1= 0.0003614614,I2= 0.0001352274,inc= 0.0006635413+0.0003614614+ 0.0001352274)#done
#Scenario2
#initial_values=c(S= 0.1794266, L1=0.007160292, L2=  0.8060518, I0=0.004158092,I1= 0.002322004,I2= 0.0008811804,inc=0.004158092+ 0.002322004+0.0008811804)#done
#Scenario3
#initial_values=c(S= 0.936477, L1=0.0002921968, L2=   0.06306177, I0=9.666564e-05,I1= 5.265821e-05,I2= 1.970011e-05,inc=  9.666564e-05+ 5.265821e-05+ 1.970011e-05)#done
#Scenario4
initial_values=c(S= 0.3176218, L1=0.00430439, L2=  0.6730307, I0=0.002848669,I1= 0.001590783,I2=  0.0006036881,inc= 0.002848669+ 0.001590783+  0.0006036881)#done

## Output Time points
times=seq(0, 4000, by = 1)

s3ICoutput=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
#View(soutput)
### fill vectors with results
s3ICvS=s3ICoutput$S
s3ICvL1=s3ICoutput$L1
s3ICvL2=s3ICoutput$L2
s3ICvI0=s3ICoutput$I0
s3ICvI1=s3ICoutput$I1
s3ICvI2=s3ICoutput$I2
s3ICvtime=s3ICoutput$time
s3ICvinc=s3ICoutput$inc
Nq=s3ICvS+s3ICvL1+s3ICvL2+s3ICvI0+s3ICvI1+s3ICvI2
plot(Nq,ylim=c(0,1))

incidence=(diff(s3ICvinc)/Nq)*100000
max(incidence)


write.csv(x=s3ICoutput,file='..//Enchik.com/output.csv')
plot(incidence)
plot(s3ICvS,xlim = c(0,10000))
s3ICvS[4000]
plot(s3ICvL1)
plot(s3ICvL2)
plot(s3ICvI0)
plot(s3ICvI1)
plot(s3ICvI2)

s3ICvS[4000]
s3ICvL1[4000]
s3ICvL2[4000]
s3ICvI0[4000]
s3ICvI1[4000]
s3ICvI2[4000]

min(s3ICvS)
max(s3ICvL1)
max(s3ICvL2)
max(s3ICvI0)
max(s3ICvI1)
max(s3ICvI2)

#Bseline
Baseline_Incidence=(diff(s3ICvinc)/Nq)*100000
max(Baseline_Incidence)
plot(Baseline_Incidence)

#scenario 1 high beta,high CDR

Baseline1_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Mass1_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Targeted1_Incidence=(diff(s3ICvinc)/s3ICvN)*100000

Scenario1=cbind(s3ICvtime,Baseline1_Incidence,Mass1_Incidence,Targeted1_Incidence)
View(Scenario1)
write.csv(x=Scenario1,file = '..//Enchik.com/Scenario1.csv')
#scenario 2 High beta low CDR
Baseline2_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Mass2_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Targeted2_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Scenario2=cbind(s3ICvtime,Baseline2_Incidence,Mass2_Incidence,Targeted2_Incidence)
View(Scenario2)
write.csv(x=Scenario2,file = '..//Enchik.com/Scenario2.csv')

#scenario 3 low beta high CDR
Baseline3_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
plot(Baseline3_Incidence,xlim = c(0,1000))
Mass3_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Targeted3_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Scenario3=cbind(s3ICvtime,Baseline3_Incidence,Mass3_Incidence,Targeted3_Incidence)
View(Scenario3)
write.csv(x=Scenario3,file = '..//Enchik.com/Scenario3.csv')
#scenario 4 low beta low CDR
Baseline4_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Mass4_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Targeted4_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
Scenario4=cbind(s3ICvtime,Baseline4_Incidence,Mass4_Incidence,Targeted4_Incidence)
View(Scenario4)
write.csv(x=Scenario4,file = '..//Enchik.com/Scenario4.csv')

#plot
#1
g_Scenario1=read.csv(file.choose())
g_Scenario1=melt(g_Scenario1,id='Time')
g_Sc1<-ggplot(g_Scenario1,aes(x=Time,y=value,colour=variable))+
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,450)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("High transmission and high detection")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

#2
g_Scenario2=read.csv(file.choose())
g_Scenario2=melt(g_Scenario2,id='Time')
g_Sc2<-ggplot(g_Scenario2,aes(x=Time,y=value,colour=variable))+
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,450)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("High transmission and low detection")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
#3
g_Scenario3=read.csv(file.choose())
g_Scenario3=melt(g_Scenario3,id='Time')

g_Sc3<-ggplot(g_Scenario3,aes(x=Time,y=value,colour=variable))+
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,450)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("Low transmission and High detection")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
#4
g_Scenario4=read.csv(file.choose())
g_Scenario4=melt(g_Scenario4,id='Time')

g_Sc4<-ggplot(g_Scenario4,aes(x=Time,y=value,colour=variable))+
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,450)+
  ylab("")+
  xlab("Time in years")+
  
  ggtitle("Low transmission and low detection")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))
g_Sc4

F0<-ggarrange( g_Sc1,g_Sc2,g_Sc3,g_Sc4,
               ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
annotate_figure(F0,left = text_grob("Incidence per 100,000 population",rot = 90))
#


#plot %
r_sc1<-read.csv(file.choose() ,header = TRUE,sep = ",",dec=".")
r_sc1=melt(r_sc1,id='Time')
View(r_sc1)
g_r_sc1<-ggplot(r_sc1, aes( x=Time,y=value,colour=variable )) +
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,50)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("Scenario 1")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


r_sc2<-read.csv(file.choose() ,header = TRUE,sep = ",",dec=".")
r_sc2=melt(r_sc2,id='Time')
View(r_sc2)
g_r_sc2<-ggplot(r_sc2, aes( x=Time,y=value,colour=variable )) +
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,50)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("Scenario 2")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

g_r_sc2


r_sc3<-read.csv(file.choose() ,header = TRUE,sep = ",",dec=".")
r_sc3=melt(r_sc3,id='Time')
View(r_sc3)
g_r_sc3<-ggplot(r_sc3, aes( x=Time,y=value,colour=variable )) +
  geom_line(size=1.5)+
  xlim(0,20)+
  ylim(0,50)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("Scenario 3")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

g_r_sc3


r_sc4<-read.csv(file.choose() ,header = TRUE,sep = ",",dec=".")
r_sc4=melt(r_sc4,id='Time')
View(r_sc4)
g_r_sc4<-ggplot(r_sc4, aes( x=Time,y=value,colour=variable )) +
  geom_line(size=1.5)+
  xlim(0,50)+
  ylim(0,50)+
  ylab("")+
  xlab("Time in years")+
  ggtitle("Scenario 4")+
  theme(legend.title = element_blank())+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

g_r_sc4

F1<-ggarrange(g_r_sc1,g_r_sc2,g_r_sc3,g_r_sc4,ncol = 4,nrow = 1,
              common.legend = TRUE, legend="bottom")
annotate_figure(F1,left = text_grob("Relative incidence reduction in %",rot = 90))































interventions_scenarios<-cbind(s3ICvtime,Baseline_Incidence,
                               Scenario1_Incidence,
                               Scenario2_Incidence,
                               Scenario3_Incidence,Scenario4_Incidence)



View(interventions_3I)

write.csv(x=interventions_scenarios,file='..//Enchik.com/interventions_scenarios.csv')

interventions_3IC2<-read.csv(file.choose())

plot (interventions_3IC2$s3ICvtime,interventions_3IC2$treeIC_equi_Incidence, type='l',col = 'red',lwd=3,main="Incidence",
      ylab="incidence per 100,000 population",xlab="Time in years",xlim = c(0,10),ylim = c(0,200))

lines(interventions_3IC2$treeIC_Incidence_sc2, col = "darkgreen",lwd=3)
#lines(interventions_3IC2$treeIC_Incidence_diff_sc2, col = "darkgreen",lwd=3)
#lines(s3ICvtime,treeIIC_combined_Incidemce_1_2, col = "blue",lwd=3)
legend(0,130, box.lty=0,c("Baseline","Mass increase CDR"),
       lty=c(1,1),lwd=c(3,3),
       col=c("red","darkgreen")) 
int_3Ic<-melt(interventions_3IC2,id="s3ICvtime")
ggplot(data=int_3Ic,mapping=aes(x=s3ICvtime,y=value,colour=variable))+
  geom_line(size=1.5)+
  ylim(0,220)+
  xlim(0,100)+
  xlab("CDR towards to Super-spreaders")+
  ylab("Incidence per 100,000 populations")

############## untargeted
#overall_cdr = 0.77

CDRs = c()
incidences = c()
for(CDR in seq(from=0.7, to=0.77, by=0.001)){
  cdr0=CDR
  cdr1=CDR
  cdr2=CDR
  delta0=(cdr0*(gamma0+mui0+mu+h)/(1-cdr0))*s0
  delta1=(cdr1*(gamma1+mui1+mu+j)/(1-cdr1))*s1
  delta2=(cdr2*(gamma2+mui2+mu)/(1-cdr2))*s2
  
  initial_values=c(S= 0.515578, L1=0.002687582, L2=  0.4797668, I0=0.001166967,I1= 0.0006412711,I2= 0.0002410786,inc= 0.001166967+ 0.0006412711+ 0.0002410786)
  ## Output Time points
  times=seq(0, 10, by = 1)
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  output=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  Nq=output$S+output$L1+output$L2+output$I0+output$I1+output$I2
  
  Incidence_all=100000*(diff(output$inc)/Nq)
  time_index = length(Incidence_all)
  Incidence_time_n = Incidence_all[10]
  incidences = c(incidences, Incidence_time_n)
  CDRs = c(CDRs, CDR)
}
out_unt<-cbind(CDRs,incidences)
write.csv(x=out_unt,file='..//Enchik.com/out_unt.7-77.csv')
out_unt<-read.csv(file.choose())
out_unt_g<-melt(out_unt,id='CDRs')
p_unt_CDR<-ggplot(data=out_unt_g,aes(x=CDRs,y=value))+
  geom_line(size=2)+
  
  xlab("Mass intervention")+
  ylab("Incidence per 100,000 populations")+
  ggtitle("A")

p_unt_CDR=p_unt_CDR+theme_bw()

p_unt_CDR
#############Targeted 
overall_cdr = 0.77

cdr2s = c()
incidences = c()
for(cdr2 in seq(from=0.77, to=0.99, by=0.01)){
  
  cdr1 = (overall_cdr-c*cdr2)/(a+b)
  cdr0 = cdr1
  
  delta0=(cdr0*(gamma0+mui0+mu+h)/(1-cdr0))*s0
  delta1=(cdr1*(gamma1+mui1+mu+j)/(1-cdr1))*s1
  delta2=(cdr2*(gamma2+mui2+mu)/(1-cdr2))*s2
  
  initial_values=c(S= 0.515578, L1=0.002687582, L2=  0.4797668, I0=0.001166967,I1= 0.0006412711,I2= 0.0002410786,inc= 0.001166967+ 0.0006412711+ 0.0002410786)
  ## Output Time points
  times=seq(0, 10, by = 1)
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  output=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  Nq=output$S+output$L1+output$L2+output$I0+output$I1+output$I2
  
  Incidence_all=100000*(diff(output$inc)/Nq)
  time_index = length(Incidence_all)
  Incidence_time_n = Incidence_all[10]
  incidences = c(incidences, Incidence_time_n)
  cdr2s = c(cdr2s, cdr2)
}
out_int<-cbind(cdr2s,incidences)
write.csv(x=out_int,file='..//Enchik.com/out_int.SS_77_99.csv')
out_int<-read.csv(file.choose())

g_out<-melt(out_int,id='cdr2s')


p_tar_CDR<-ggplot(data=g_out,aes(x=cdr2s,y=value))+
  geom_line(size=2)+
  
  xlab("CDR, Super-spreaders ")+
  ylab("Incidence per 100,000 populations")+
  ggtitle("B")

p_tar_CDR=p_tar_CDR+theme_bw()
p_tar_CDR

ggarrange(p_unt_CDR,p_tar_CDR,  ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

scale_colour_manual(values=c("red","blue","purple"),name="Intervention",
                    
                    labels=c("Baseline","Mass","Targeted"))
###################
#############Targe NEW

intervention_cdr = 0.77#intervention cdr
cdr2_max=0.95# the maximum possible case detection rate which is min[((intervention_cdr -(a+b)*baseline_cdr))/c,Maximum_cdr]
#we consider we can go beyond 95%

ds = c()
incidences = c()
for(d in seq(from=0.0, to=1, by=0.01)){
  cdr2=intervention_cdr+(cdr2_max-intervention_cdr )*d
  cdr1 = (intervention_cdr-c*cdr2)/(a+b)
  cdr0 = cdr1
  
  delta0=(cdr0*(gamma0+mui0+mu+h)/(1-cdr0))*s0
  delta1=(cdr1*(gamma1+mui1+mu+j)/(1-cdr1))*s1
  delta2=(cdr2*(gamma2+mui2+mu)/(1-cdr2))*s2
  
  initial_values=c(S= 0.515578, L1=0.002687582, L2=  0.4797668, I0=0.001166967,I1= 0.0006412711,I2= 0.0002410786,inc= 0.001166967+ 0.0006412711+ 0.0002410786)
  ## Output Time points
  times=seq(0, 10, by = 1)
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  output=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  Nq=output$S+output$L1+output$L2+output$I0+output$I1+output$I2
  
  Incidence_all=100000*(diff(output$inc)/Nq)
  time_index = length(Incidence_all)
  Incidence_time_n = Incidence_all[10]
  incidences = c(incidences, Incidence_time_n)
  ds = c(ds, d)
}
out_int<-cbind(ds,incidences)
plot(out_int$ds,out_int$incidences)
write.csv(x=out_int,file='..//Enchik.com/out_int.SS_d.csv')
out_int<-read.csv(file.choose())

g_out<-melt(out_int,id='ds')


p_tar_CDR<-ggplot(data=g_out,aes(x=ds,y=value))+
  geom_point(size=2,color="dodgerblue2")+
  
  xlab("proportion of intervention targeted")+
  ylab("Incidence per 100,000 populations")+
  ggtitle("B")

p_tar_CDR=p_tar_CDR+theme_bw()
p_tar_CDR
##########
#for LHS min and max value
intervention_cdr = 0.77#intervention cdr
cdr2_max=0.95# 
delta0is=c()
delta1is=c()
delta2is=c()
for(d in seq(from=0.0, to=1, by=0.01)){
  cdr2=intervention_cdr+(cdr2_max-intervention_cdr )*d
  cdr1 = (intervention_cdr-c*cdr2)/(a+b)
  cdr0 = cdr1
  
  delta0i=(cdr0*(gamma0+mui0+mu+h)/(1-cdr0))*s0
  delta1i=(cdr1*(gamma1+mui1+mu+j)/(1-cdr1))*s1
  delta2i=(cdr2*(gamma2+mui2+mu)/(1-cdr2))*s2
  delta0is=c(delta0is,delta0i)
  delta1is=c(delta1is,delta1i)
  delta2is=c(delta2is,delta2i)
  
}
deltais=cbind(delta0is,delta1is,delta2is)
write.csv(x=deltais,file = '..//Enchik.com/deltais.csv')
View(deltais)# see max and min
####################IPT intervention Untargeted
#50% coverage=b*I1+c*I2=0.5

#overall_coverage = 0.5

IPTs = c()
incidences = c()
for(IPT in seq(from=0.0, to=0.25, by=0.01)){
  
  cov_I1=IPT
  cov_I2=IPT
  
  IPT_effectivness=0.8
  p1=cov_I1*IPT_effectivness#rate_of_infections_from_I1_prevented_by_IPT_over_1_year_
  p2=cov_I2*IPT_effectivness
  
  initial_values=c(S= 0.515578, L1=0.002687582, L2=  0.4797668, I0=0.001166967,I1= 0.0006412711,I2= 0.0002410786,inc= 0.001166967+ 0.0006412711+ 0.0002410786)
  ## Output Time points
  times=seq(0, 10, by = 1)
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  output=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  Nq=output$S+output$L1+output$L2+output$I0+output$I1+output$I2
  
  Incidence_all=100000*(diff(output$inc)/Nq)
  time_index = length(Incidence_all)
  Incidence_time_n = Incidence_all[10]
  incidences = c(incidences, Incidence_time_n)
  IPTs = c(IPTs, IPT)
}
out_ipt<-cbind(IPTs,incidences)
View(out_ipt)
write.csv(x=out_ipt,file='..//Enchik.com/out_ipt_25%.csv')
out_ipt<-read.csv(file.choose())



G_ipt=melt(out_ipt,id='IPTs')
p_unt_IPT<-ggplot(data=G_ipt,aes(x=IPTs,y=value))+
  geom_line(size=1.5)+
  
  xlab("IPT mass coverage")+
  ylab("Incidence per 100,000 populations")+
  ggtitle("A")

p_unt_IPT=p_unt_IPT+theme_bw()
p_unt_IPT
####################IPT Targeted intervention
#50% coverage=b*I1+c*I2=0.5

overall_coverage = 0.25

cov_I2s = c()
incidences = c()
for(cov_I2 in seq(from=0.25, to=0.9, by=0.01)){
  
  cov_I1= (overall_coverage*(b+c)-(cov_I2*c))/b
  
  IPT_effectivness=0.8
  p1=cov_I1*IPT_effectivness#rate_of_infections_from_I1_prevented_by_IPT_over_1_year_
  p2=cov_I2*IPT_effectivness
  
  initial_values=c(S= 0.515578, L1=0.002687582, L2=  0.4797668, I0=0.001166967,I1= 0.0006412711,I2= 0.0002410786,inc= 0.001166967+ 0.0006412711+ 0.0002410786)
  ## Output Time points
  times=seq(0, 10, by = 1)
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  output=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  Nq=output$S+output$L1+output$L2+output$I0+output$I1+output$I2
  
  Incidence_all=100000*(diff(output$inc)/Nq)
  time_index = length(Incidence_all)
  Incidence_time_n = Incidence_all[10]
  incidences = c(incidences, Incidence_time_n)
  cov_I2s = c(cov_I2s, cov_I2)
}
out_ipt_Target<-cbind(cov_I2s,incidences)
View(out_ipt_Target)
write.csv(x=out_ipt_Target,file='..//Enchik.com/out_ipt_SS.csv')
out_ipt_Ta_SS<-read.csv(file.choose())

G_ipt_SS=melt(out_ipt_Ta_SS,id='cov_I2s')
p_tar_IPT<-ggplot(data=G_ipt_SS,aes(x=cov_I2s,y=value))+
  geom_line(size=1.5)+
  
  xlab("IPT Super-spreaders coverage")+
  ylab("Incidence per 100,000 populations")+
  ggtitle("B")

p_tar_IPT=p_tar_IPT+theme_bw()
p_tar_IPT

###################
library(ggpubr) 
ggarrange(p_unt_IPT, p_tar_IPT, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

scale_colour_manual(values=c("red","blue","purple"),name="Intervention",
                    
                    labels=c("Baseline","Mass","Targeted"))



#Simulate the S1S2L1L2IT transmission 
s3ICoutput=as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
#View(soutput)
### fill vectors with results
s3ICvS=s3ICoutput$S
s3ICvL1=s3ICoutput$L1
s3ICvL2=s3ICoutput$L2
s3ICvI0=s3ICoutput$I0
s3ICvI1=s3ICoutput$I1
s3ICvI2=s3ICoutput$I2
s3ICvtime=s3ICoutput$time
s3ICvinc=s3ICoutput$inc


s3ICvN=s3vS+s3vL1+s3vL2+s3vI0+s3vI1+s3vI2
plot(s3ICvtime,s3ICvinc
     ,ylim = c(0,5))
plot(s3ICvS)
plot(s3ICvL1)
plot(s3ICvL2)
plot(s3ICvI0)
plot(s3ICvI1)
plot(s3ICvI2)


min(s3ICvS)
max(s3ICvL1)
max(s3ICvL2)
max(s3ICvI0)
max(s3ICvI1)
max(s3ICvI2)

s3ICIncidence_tri=(s3ICvinc)*100000
max(s3ICIncidence_tri)#218.0602

#baseline 
#Intervention comparion by models
#Bseline
treeIC_equi_Incidence=(diff(s3ICvinc)/s3ICvN)*100000
max(treeIC_equi_Incidence)
#scenario 1
treeIC_Incidence_sc1=(diff(s3ICvinc)/s3ICvN)*100000
max(treeIC_Incidence_md)
# ###sxcenario 2
treeIC_Incidence_diff_sc2=(diff(s3ICvinc)/s3ICvN)*100000
plot(treeIC_equi_Incidence)


interventions_3IC<-cbind(s3ICvtime,treeIC_equi_Incidence,treeIC_Incidence_sc1,treeIC_Incidence_diff_sc2)
View(interventions_3I)
write.csv(x=interventions_3IC,file='..//Enchik.com/interventions_3I.csv')
interventions_3IC2<-read.csv(file.choose())

plot (interventions_3IC2$s3ICvtime,interventions_3IC2$treeIC_equi_Incidence, type='l',col = 'red',lwd=3,main="Incidence",
      ylab="incidence per 100,000 population",xlab="Time in years",xlim = c(0,10),ylim = c(190,220))

lines(interventions_3IC2$treeIC_Incidence_sc1, col = "brown",lwd=3)
lines(interventions_3IC2$treeIC_Incidence_diff_sc2, col = "darkgreen",lwd=3)
#lines(s3ICvtime,treeIIC_combined_Incidemce_1_2, col = "blue",lwd=3)
legend(0,200, box.lty=0,c("Baseline","Scenario 1","Scenario 2"),
       lty=c(1,1,1),lwd=c(3,3,3),
       col=c("red","brown","darkgreen")) 

#Ro
#=====================================================================
#SENSITIVITY ANALYSIS
#Latin hypercube sampling
#install.packages('lhs')

a=0.58 #proportion of Io
b=0.31#proportion of I1
c=0.11##proportion of I2
require(lhs) #add the lhs library
z <- 1 #choose number of points to simulate
set.seed(6242015)#random number generator
lhs<-maximinLHS(z,25) #simulate h= number of simulations, 19=number of parameters
#To map these points in the unit cube to our parameters, we need minimum and maximum values for each.
mu.min= 0.01538462
mu.max= 0.01538462
mui0.min=0.0256
mui0.max=0.0256
mui1.min=0.0256
mui1.max=0.0256
mui2.min=0.178
mui2.max=0.178
r.min=0.21
r.max=0.21
beta0.min<-0
beta0.max<-0
beta1.min<-12.54
beta1.max<-12.54
beta2.min<-57
beta2.max<-57
epsilon0.min<-0.402*a
epsilon0.max<-0.402*a
epsilon1.min<-0.402*b
epsilon1.max<-0.402*b
epsilon2.min<-0.402*c
epsilon2.max<-0.402*c
kappa.min<-3.65
kappa.max<-3.65

gamma0.min<-0.292
gamma0.max<-0.292
gamma1.min<-0.292
gamma1.max<-0.292
gamma2.min<-0.14
gamma2.max<-0.14

nu0.min<-0.002*a
nu0.max<-0.002*a
nu1.min<-0.002*b
nu1.max<-0.002*b
nu2.min<-0.002*c
nu2.max<-0.002*c

delta0.min<-0.6495713#0.9319936#cdr=77%
delta0.max<-0.6495713#0.9319936
delta1.min<-0.6495713#0.9319936
delta1.max<-0.6495713#0.9319936
delta2.min<-0.9471179#1.358908#
delta2.max<-0.9471179#1.358908#

h.min<-0.015
h.max<-0.015
j.min<-0.015
j.max<-0.015
p1.min<-0
p1.max<-0
p2.min<-0
p2.max<-0



#Now we can generate a “parameter set” by rescaling our simulated latin hypercube sample
params.set <- cbind(
  
  beta0 = lhs[,1]*(beta0.max-beta0.min)+beta0.min,
  beta1 = lhs[,2]*(beta1.max-beta1.min)+beta1.min,
  beta2 = lhs[,3]*(beta2.max-beta2.min)+beta2.min,
  epsilon0 =lhs[,4]*(epsilon0.max-epsilon0.min)+epsilon0.min,
  epsilon1 =lhs[,5]*(epsilon1.max-epsilon1.min)+epsilon1.min,
  epsilon2 =lhs[,6]*(epsilon2.max-epsilon2.min)+epsilon2.min,
  kappa = lhs[,7]*(kappa.max-kappa.min)+kappa.min,
  gamma0 = lhs[,8]*(gamma0.max-gamma0.min)+gamma0.min,
  gamma1 = lhs[,9]*(gamma1.max-gamma1.min)+gamma1.min,
  gamma2 = lhs[,10]*(gamma2.max-gamma2.min)+gamma2.min,
  nu0 = lhs[,11]*(nu0.max-nu0.min)+nu0.min,
  nu1 = lhs[,12]*(nu1.max-nu1.min)+nu1.min,
  nu2 = lhs[,13]*(nu2.max-nu2.min)+nu2.min,
  mu = lhs[,14]*(mu.max-mu.min)+mu.min,
  mui0 = lhs[,15]*(mui0.max-mui0.min)+mui0.min,
  mui1 = lhs[,16]*(mui1.max-mui1.min)+mui1.min,
  mui2 = lhs[,17]*(mui2.max-mui2.min)+mui2.min,
  r = lhs[,18]*(r.max-r.min)+r.min,
  delta0 = lhs[,19]*(delta0.max-delta0.min)+delta0.min,
  delta1 = lhs[,20]*(delta1.max-delta1.min)+delta1.min,
  delta2 = lhs[,21]*(delta2.max-delta2.min)+delta2.min,
  h = lhs[,22]*(h.max-h.min)+h.min,
  j = lhs[,23]*(j.max-j.min)+j.min,
  p1 = lhs[,24]*(p1.max-p1.min)+p1.min,
  p2 = lhs[,25]*(p2.max-p2.min)+p2.min)

View(params.set)
#Test
params <- as.list(c(params.set[1,]))
initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
times=seq(0, 2000, by = 1)
B_out <- as.data.frame(lsoda(initial_values, times, threeIIC_model, params))
plot(B_out$time,B_out$S)
Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
incidence_B_out=(diff(B_out$inc)/Nq)*100000
max(incidence_B_out)
#creat matrix to save whole info
output_matrix = data.frame(params.set)
B_incidence = data.frame('B_incidence'=rep(NA,z))
Int_incidence = data.frame('Int_incidence'=rep(NA,z))
Change_incidence= data.frame('Change_incidence'=rep(NA,z))
output_matrix = cbind(output_matrix, B_incidence,Int_incidence,Change_incidence) #add incidence column
View(output_matrix)
#(#THese are all the simulated parameters we need to exlore how these parameters affect the incidence
#To facilitate later plotting, we decide how many different levels of I to consider.
#levels <- 15
#Even though we have simulated 100 points, to speed up computing we will use only a fraction of these
#for this demonstration.
#h2 <-30

#we set up a nested loop, first to cycle through different values of I, then to cycle through
#different simulated parameter sets. Note the pre-allocated data frame and use of the counter j, also to
#speed up evaluation

for(i in 1:z){
  #run baseline 
  
  
  initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
  cdr0=0.7
  cdr1=0.7
  cdr2=0.7
  #scenario 2, extent of targeting
  # code for all CDRS channgins as cdr2 changes and this change deltas
  delta0=(cdr0*(gamma0+mui0+mu+h)/(1-cdr0))*s0
  delta1=(cdr1*(gamma1+mui1+mu+j)/(1-cdr1))*s1
  delta2=(cdr2*(gamma2+mui2+mu)/(1-cdr2))*s2
  
  parameter_list=c(N=N,mu=mu,beta0=beta0, beta1=beta1, epsilon0=epsilon0,epsilon1=epsilon1,epsilon2=epsilon2, 
                   kappa=kappa,gamma0=gamma0,gamma1=gamma1,gamma2=gamma2,nu0=nu0,nu1=nu1,nu2=nu2, 
                   mui0=mui0,mui1=mui1,mui2=mui2,delta0=delta0,delta1=delta1,delta2=delta2,r=r,h=h,j=j,p1=p1,p2=p2)  
  
  times=seq(0, 2000, by = 1)
  B_out <- as.data.frame(lsoda(initial_values, times, threeIIC_model, parameter_list))
  #Record Baseline equilibrium incidence
  Nq=B_out$S+B_out$L1+B_out$L2+B_out$I0+B_out$I1+B_out$I2
  incidence_B_out=(diff(B_out$inc)/Nq)*100000
  B_Incidence_time_n = incidence_B_out[2000] #the model reaches equilibrium at time around 500
  output_matrix$B_incidence[i] = B_Incidence_time_n
  #record compartment values at equilibrium
  Equi_initial_values=c(S=min(B_out$S), L1=max(B_out$L1), L2=max(B_out$L2), I0=max(B_out$I0),I1=max(B_out$I1),I2=max(B_out$I2),inc=max(B_out$I0)+max(B_out$I1)+max(B_out$I2))
  #run intervention model, the intervention is mass 10% increase CDR i.e delta0= 0.9319936,delta1=0.9319936,delta2=1.358908
  
  params <- as.list(c(params.set[i,]))
  
  times=seq(0, 2000, by = 1)
  Int_out <- as.data.frame(lsoda(initial_values, times, threeIIC_model, params))
  Nq=Int_out$S+Int_out$L1+Int_out$L2+Int_out$I0+Int_out$I1+Int_out$I2
  Int_incidence_out=(diff( Int_out$inc)/Nq)*100000
  #considering a 50 years impact of intervention
  Int_Incidence_time_n = Int_incidence_out[2000]
  output_matrix$Int_incidence[i] = Int_Incidence_time_n
  #calculate the difference between baseline incidence and intervention incidence
  output_matrix$Change_incidence[i]= (output_matrix$B_incidence[i])-(output_matrix$Int_incidence[i])
  
} 


View(output_matrix) #now we have incidence and each parameters in one matrix










for(i in 1:z){
  
  #initial_values=c(S=0.5418456958, L1=0.0197116035, L2=0.4348963153, I=0.0026741092, T=0.0008722763,inc=0.0026741092,notif=0.0008722763)#equilibrium
  initial_values=c(S=A-(D+E+F), L1=B, L2=C, I0=D,I1=E,I2=F,inc=D+E+F)
  
  params <- as.list(c(params.set[i,]))
  times=seq(0, 2000, by = 1)
  out <- as.data.frame(lsoda(initial_values, times, threeIIC_model, params))
  Nq=out$S+out$L1+out$L2+out$I0+out$I1+out$I2
  incidence_out=(diff(out$inc)/Nq)*100000
  
  Incidence_time_n = incidence_out[1000]
  output_matrix$B_incidence[i] = Incidence_time_n
}
View(Incidence_time_n)
View(output_matrix) #now we have incidence and each parameters in one matrix

plot(output_matrix$h,output_matrix$incidence,type = 'p')

#Impact of parameters on interevention 




















#Plot impact of h and j

#The general approach is to convert the data to long format (using melt() from package reshape or reshape2) or gather() from the tidyr package:
write.csv(x=output_matrix,file='..//Enchik.com/impact_h_j.csv')
impact<-read.csv(file.choose())
View(impact)
plot(impact$h,impact$incidence,type = 'l',col='red',
     ylab = "incidence per 100,000 population",xlab = "conversion parameter",ylim = c(100,250))
lines(impact$j,impact$incidence,type = 'l',col='blue')
legend(0.001,180, box.lty=0,c("h","j"),lty=c(1,1),lwd=c(3,3),
       col=c("red","blue"))

library("reshape2")
library("ggplot2")
test_hj<-melt(impact,id="sample")
ggplot(data = test_hj,aes(x=h,y=incidence,colours=variable))+
  geom_line()




write.csv(x=output_matrix,'..//Enchik.com/sensitivity.csv')
plot(output_matrix$delta2,output_matrix$incidence,type="l")

#which parameter affects more,inspecting the partial correlation
#Partial rank correlations can be computed using the pcc function in the R package sensitivity 
#install.packages('sensitivity')
library(sensitivity)
bonferroni.alpha <- 0.05

prcc <- pcc(output_matrix[,1:19], output_matrix[,20], nboot = 1000, rank=TRUE, conf=1-bonferroni.alpha)
save(prcc, file='prcc.Rdata')

#We can view a table of the resulting partial correlation coefficients. if none of the (penalized)
#confidence intervals contains zero, we conclude that all are significant and produce a plot showing their
#relative magnitudes.
library(sensitivity)
load('prcc.Rdata')
summary <- print(prcc)
Corl=data.frame(summary) 
View(Corl)

write.csv(x=Corl,file='..//Enchik.com//corl6.1.csv')
#plote the partial corelation coeeficient

par(mar=c(9,4,4,2)+0.1)
plot(Corl$original, main='Partial rank correlation coefficients', ylim=c(-1,1),
     xlab='', ylab='Coefficient',
     axes=FALSE)
axis(2)
axis(1, at=seq(1:19), labels=row.names(Corl), las=2)
mtext(text='Parameter', side=1, line=4.5)
box()
for(i in 1:19) lines(c(i,i),c(Corl[i,4], Corl[i,5]))
abline(h=0)
#tornado plot


#We can plot each simulated value as a point
par(mfrow=c(4,4))
#plot(output_matrix$beta0, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
# xlab='beta0',
#ylab='Incidence',main = 'PRCC= 0.017')

plot(output_matrix$beta1, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(beta1),
     ylab='Incidence',main = 'PRCC= 0.486')
plot(output_matrix$beta2, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(beta2),
     ylab='Incidence',main = 'PRCC= 0.940')
plot(output_matrix$epsilon0, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(epsilon0),
     ylab='Incidence',main = 'PRCC= 0.420')
plot(output_matrix$epsilon1, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(epsilon1),
     ylab='Incidence',main = 'PRCC= 0.509')

plot(output_matrix$epsilon2, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(epsilon2),
     ylab='Incidence',main = 'PRCC= 0.769')


plot(output_matrix$kappa, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab='kappa',
     ylab='Incidence',main = 'PRCC= -0.013')

plot(output_matrix$gamma, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='blue',
     xlab=expression(gamma),
     ylab='Incidence',main = 'PRCC= 0.016')

plot(output_matrix$nu0, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='coral4',
     xlab=expression(nu0),
     ylab='Incidence',main = 'PRCC= 0.030')

plot(output_matrix$nu1, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='coral4',
     xlab=expression(nu1),
     ylab='Incidence',main = 'PRCC= 0.104')

plot(output_matrix$nu2, output_matrix$incidence, type = 'p' ,lwd=2,pch=19, cex=0.9, col='coral4',
     xlab=expression(nu2),
     ylab='Incidence',main = 'PRCC= 0.190')

