

df_parameter=data.frame("Battery_ID" = c(rep("A1",4),rep("A2",4),rep("A3",4)), "Parameter" = rep(c('a','b','c','d'),3), 
                        "Low_bound" = c(-1.66*10^-3,2.068*10^-2,9.079*10^-1,-1.21*10^-3,-2.007*10^-6,5.283*10^-2,8.931*10^-1,-9.007*10^-4,-3.788*10^-5,5.398*10^-2,8.631*10^-1,-1.188*10^-3),
                        "Mean"=c(-1.042*10^-3,2.268*10^-2,9.190*10^-1,-1.035*10^-3,-9.860*10^-7,5.752*10^-2,8.983*10^-1,-8.34*10^-4,-1.53*10^-5,6.296*10^-2,8.757*10^-1,-9.4*10^-4),
                        "Upper_Bound"=c(-4.24*10^-4,2.467*10^-2,9.301*10^-1,-8.6*10^-4,3.442*10^-8,6.221*10^-2,9.035*10^-1,-7.67*10^-4,7.272*10^-6,7.193*10^-2,8.883*10^-1,-6.92*10^-4))

#add sigma column, for the variability of each parameter

possible_sigma=data.frame("with_inf"=(df_parameter$Mean-df_parameter$Low_bound)/2,"with_sup"=(df_parameter$Upper_Bound-df_parameter$Mean)/2)
df_parameter=cbind(df_parameter, "variability"=apply(possible_sigma,1,mean))

df_splitted =split(df_parameter, df_parameter$Battery_ID)
generate_data <- function(sub_df,n_cycles=800) {
  data_by_battery=data.frame("cycle"=seq(0,n_cycles-1), "a"=rnorm(n_cycles,sub_df[1,4],sub_df[1,6]),"b"=rnorm(n_cycles,sub_df[2,4],sub_df[2,6]),
                             "c"=rnorm(n_cycles,sub_df[3,4],sub_df[3,6]),"d"=rnorm(n_cycles,sub_df[4,4],sub_df[4,6]))
  
  data_by_battery=cbind(data_by_battery,"Q"=data_by_battery$a*exp(data_by_battery$cycle*data_by_battery$b)+data_by_battery$c*exp(data_by_battery$cycle*data_by_battery$d))
  return (data_by_battery)
}


n_cycles=200
#we reduce variability of table 3, to avoid too unstable values
df_splitted[[3]][2,6]=df_splitted[[3]][2,6]/1000
df_splitted[[3]][4,6]=df_splitted[[3]][4,6]/1000
liste_curve_data=lapply(X = df_splitted,generate_data,n_cycles=n_cycles)



library(nlstools)
par(mfrow=c(3,1))
capacity_formula =as.formula(Q ~ a*exp(cycle*b) +c*exp(cycle*d))
for (num_plot in seq(1,3)) {
  data=liste_curve_data[[num_plot]][,c("cycle","Q")]
  #preview(capacity_formula, data = data,
  #        start = list(a = df_splitted[[num_plot]][1,4], b =df_splitted[[num_plot]][2,4], c =df_splitted[[num_plot]][3,4],d=df_splitted[[num_plot]][4,4]))
  fitted_curb= nls(capacity_formula,start = list(a = df_splitted[[num_plot]][1,4], b =df_splitted[[num_plot]][2,4], c =df_splitted[[num_plot]][3,4],d=df_splitted[[num_plot]][4,4]), data = data)
  overview(fitted_curb)
  plotfit(fitted_curb, smooth = TRUE,las=1)
  #plot(liste_curve_data[[num_plot]]$Q~liste_curve_data[[num_plot]]$cycle,type="p",col=num_plot+1,xlab="Cycle",ylab="Capacity(Ah)",main=names(liste_curve_data)[num_plot],ylim=c(0.5,1),xlim=c(0,n_cycles),pch=3,cex=0.5,las=1)
}

#function for interval inclusion
is.included <- function(ref_interval,temp_interval) {
  if (temp_interval[1]>=ref_interval[1] & temp_interval[2]<=ref_interval[2]) {
    return (TRUE)
  }
  else {return (FALSE)}
}


#computation of dempster-shafer theory

#compuation of belief measure

compute_belief_measure=function(df_splitted) {
  number_row=length(df_splitted)
  number_column=nrow(df_splitted[[1]])
  belief_measure=matrix(1, nrow = number_row, ncol = number_column,dimnames=list(names(df_splitted),as.character(df_splitted[[1]][,"Parameter"])))
  for (row in 1:number_row) {
    for (column in 1:number_column) {
      bel_sum=0
      ref_interval=c(df_splitted[[row]][column,3], df_splitted[[row]][column,5])
      paste("intervalle de refrence est ",ref_interval)
      for (interval in 1:3) {
        temp_interval=c(df_splitted[[interval]][column,3], df_splitted[[interval]][column,5])
        if (is.included(ref_interval,temp_interval)) {
          bel_sum=bel_sum+1/number_row
        }
      }
      belief_measure[row,column]=bel_sum
    }
  }
  return (belief_measure)
  }

belief_measure=compute_belief_measure(df_splitted)
#computation of basic belief

compute_basic_belief=function(belief_df,df_splitted) {
  number_row=nrow(belief_df)
  number_column=ncol(belief_df)
  basic_belief=matrix(1, nrow = number_row, ncol = number_column,dimnames = list(rownames(belief_df),colnames(belief_df)))
  for (row in 1:number_row) {
    for (column in 1:number_column) {
      basic_bel_sum=0
      ref_interval=c(df_splitted[[row]][column,3], df_splitted[[row]][column,5])
      for (interval in 1:number_row) {
        temp_interval=c(df_splitted[[interval]][column,3], df_splitted[[interval]][column,5])
        if (is.included(ref_interval,temp_interval)) {
          basic_bel_sum=basic_bel_sum+belief_df[interval,column]
        }
      }
      basic_belief[row,column]=basic_bel_sum
    }
  }
  
  basic_belief=apply(basic_belief,2,function(x) x/sum(x)) 
  return (basic_belief)
}
basic_belief=compute_basic_belief(belief_measure,df_splitted)

#computation of initial parameters of mc carlo
parameter_table= matrix(df_parameter$Mean,nrow=3,ncol=4,byrow = TRUE,dimnames=list(rownames(basic_belief),colnames(basic_belief)))
initial_parameters=matrix(apply(basic_belief*parameter_table,2,sum),nrow=1,ncol=ncol(parameter_table),dimnames = list("1",colnames(parameter_table)))



generate_sequential_data=function(sigma_vector,init_parameter,num_cycle=200,var_capacity=0.0002) {
  df_generated=data.frame(init_parameter)
  
  for (cycle in 2:num_cycle) {
    df_generated=rbind(df_generated,df_generated[cycle-1,]+rnorm(4,0,sd=sigma_vector))
  }
  
  df_generated=cbind(df_generated,"cycle"=seq(0,num_cycle-1))
  df_generated=cbind(df_generated,"Q"=df_generated$a*exp(df_generated$cycle*df_generated$b)+df_generated$c*exp(df_generated$cycle*df_generated$d)+rnorm(num_cycle,0,var_capacity))
  rownames(df_generated)=1:num_cycle
  return (df_generated)
  }


df_batteryA4=generate_sequential_data(rep(0.000002,4),initial_parameters)


#build bootstrapfilter with nimble
library(nimble)


## define model with a priori logit
modele_decay_battery <- nimbleCode({
  
  #use of logit p, to be more uniform on the interval
  logit_sigA~dnorm(-log(10),sd=log(10)/2)
  logit_sigB~dnorm(-log(10),sd=log(10)/2)
  logit_sigC~dnorm(-log(10),sd=log(10)/2)
  logit_sigD~dnorm(-log(10),sd=log(10)/2)
  logit_sigQ~dnorm(-log(10),sd=log(10)/2)
  
  sigA<-exp(logit_sigA)/exp(1+logit_sigA)
  sigB<-exp(logit_sigB)/exp(1+logit_sigB)
  sigC<-exp(logit_sigC)/exp(1+logit_sigC)
  sigD<-exp(logit_sigD)/exp(1+logit_sigD)
  sigQ<-exp(logit_sigQ)/exp(1+logit_sigQ)
  
  #with x0 the parameters determined with DST method
  x[1,1]~ dnorm(x0[1] , sd = 1e-02)
  x[1,2]~ dnorm(x0[2] , sd = 1e-02)
  x[1,3]~ dnorm(x0[3] , sd = 1e-02)
  x[1,4]~ dnorm(x0[4] , sd = 1e-02)
  
  for (i in 2:num_cycle) {
    x[i,1] ~ dnorm(x[i - 1,1] , sd = sigA)
    x[i,2] ~ dnorm(x[i - 1,2] , sd = sigB)
    x[i,3] ~ dnorm(x[i - 1,3] , sd = sigC)
    x[i,4] ~ dnorm(x[i - 1,4] , sd = sigD)
    
    Q[i] ~ dnorm(x[i,1]*exp(x[i,2]*(i-1))+x[i,3]*exp(x[i,4]*(i-1)), sd = sigQ)
  }
})


## define data, constants, and initial values  
data <- list(  Q = df_batteryA4$Q)


constants <- list(num_cycle = 200,x0=as.vector(initial_parameters))
inits <- list(logit_sigA = log(1e-02/(1-1e-02)), logit_sigB = log(1e-02/(1-1e-02)), logit_sigC = log(1e-02/(1-1e-02)),logit_sigD = log(1e-02/(1-1e-02)),logit_sigQ = log(1e-02/(1-1e-02)))

## build the model
stateSpaceModel <- nimbleModel(modele_decay_battery, data = data, constants = constants,inits = inits,check = TRUE)

### build bootstrap filter and compile model and filter
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = 'x')
compiledList <- compileNimble(stateSpaceModel, bootstrapFilter,showCompilerOutput = TRUE)

#run with 100 particles
compiledList$bootstrapFilter$run(10000)

posteriorSamples <- as.matrix(compiledList$bootstrapFilter$mvEWSamples)

#define mean Q for several charges
#with monte carlo resampling, all weights are already balanced
Q_a_posteriori=posteriorSamples[,1]*exp(200*posteriorSamples[,2])+posteriorSamples[,3]*exp(200*posteriorSamples[,4])


#define model with uniform a priori

modele_decay_battery_uniform <- nimbleCode({
  
  sigA~dunif(1e-04,1)
  sigB~dunif(1e-04,1)
  sigC~dunif(1e-04,1)
  sigD~dunif(1e-04,1)
  sigQ~dunif(1e-04,1)
  
  #with x0 the parameters determined with DST method
  x[1,1]~ dnorm(x0[1] , sd = 1e-02)
  x[1,2]~ dnorm(x0[2] , sd = 1e-02)
  x[1,3]~ dnorm(x0[3] , sd = 1e-02)
  x[1,4]~ dnorm(x0[4] , sd = 1e-02)
  
  for (i in 2:num_cycle) {
    x[i,1] ~ dnorm(x[i - 1,1] , sd = sigA)
    x[i,2] ~ dnorm(x[i - 1,2] , sd = sigB)
    x[i,3] ~ dnorm(x[i - 1,3] , sd = sigC)
    x[i,4] ~ dnorm(x[i - 1,4] , sd = sigD)
    
    Q[i] ~ dnorm(x[i,1]*exp(x[i,2]*(i-1))+x[i,3]*exp(x[i,4]*(i-1)), sd = sigQ)
  }
})


other_init=c(-0.00040,0.075,0.75,-0.00050)
constants <- list(num_cycle = 100,x0=other_init)
inits <- list(sigA = 1e-02, sigB = 1e-02,sigC = 1e-02,sigD = 1e-02,sigQ = 1e-02)

## build the model
stateSpaceModel <- nimbleModel(modele_decay_battery_uniform, data = data, constants = constants,inits = inits,check = TRUE)

### build bootstrap filter and compile model and filter
bootstrapFilter <- buildBootstrapFilter(stateSpaceModel, nodes = 'x')
compiledList <- compileNimble(stateSpaceModel, bootstrapFilter,showCompilerOutput = TRUE)

#run with 100 particles
compiledList$bootstrapFilter$run(10000)

posteriorSamples <- as.matrix(compiledList$bootstrapFilter$mvEWSamples)

#define mean Q for several charges
#with monte carlo resampling, all weights are already balanced
Q_a_posteriori_uniform=posteriorSamples[,1]*exp(100*posteriorSamples[,2])+posteriorSamples[,3]*exp(100*posteriorSamples[,4])




dag1 = model2network("[A][B|A][C|A]")
dag2 = model2network("[A|B:C][B][C]")

shd(dag1,dag2)
















#use MCMC to determine parameters
## create MCMC specification for the state space model
stateSpaceMCMCconf <- configureMCMC(stateSpaceModel, nodes = NULL)

## add a block pMCMC sampler for sigma Q
stateSpaceMCMCconf$addSampler(target = c('sigQ','sigA'),
                              type = 'RW_PF_block', control = list(latents = 'x'))

## build and compile pMCMC sampler
stateSpaceMCMC <- buildMCMC(stateSpaceMCMCconf)
compiledListMCMC <- compileNimble(stateSpaceModel, stateSpaceMCMC,'showCompilerOutput = TRUE', resetFunctions = TRUE)

## run compiled sampler for 200 iterations
compiledListMCMC$stateSpaceMCMC$run(100)
library('coda')
posteriorSamps <- as.mcmc(as.matrix(compiledListMCMC$stateSpaceMCMC$mvSamples))


traceplot(posteriorSamps[,'sigQ'], ylab = 'sigQ')


