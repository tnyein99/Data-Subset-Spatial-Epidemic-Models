library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(EpiILM)

##input: population size (sam), cut-off infection time (sub_time)
##output: absolute biases
temp_sub<-function(sam,sub_time){
  x<- runif(sam,0,sqrt(sam))
  y<- runif(sam,0,sqrt(sam))
  lambda <- rep(3, sam)
  
  alphapar2 <- c(1, 1)
  betapar2 <- c(1, 1)
  
  #generate the epidemic data until the convergence to the posterior distribution has reached
  repeat{
    sir_mod0 <- epidata(type = "SIR", n = sam, tmax = 15, sus.par = 0.2, beta = 2, infperiod = lambda,
                        x = x, y = y)
    
    mcmc_out100p<-epimcmc(object = sir_mod0, tmin = min(sir_mod0$inftime), tmax = max(sir_mod0$inftime),
                          niter = 50000, sus.par.ini = 1, beta.ini = 1,
                          Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                          prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                          prior.sus.par = alphapar2, prior.beta.par = betapar2,
                          adapt = FALSE, acc.rate = NULL)
    test<-as.mcmc(mcmc_out100p[[3]])
    gt<-geweke.diag(test, frac1=0.1, frac2=0.5)
    pval<-pnorm(abs(gt$z),lower.tail=FALSE)*2 
    
    if(any(is.na(pval))==FALSE)
      break
    if (length(pval[pval > 0.05])<1)
      break 
  }
  
  
  
  
  sum_100p_spatial<-summary(mcmc_out100p, start = 10001)
  
  
  ##extract the posterior mean and sd of entire epidemic data
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_sd_100p<-sum_100p_spatial$statistics[1,2]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_sd_100p<-sum_100p_spatial$statistics[2,2]
  ##run the mcmc with maximum infection time as desired cut-off points 
  mcmc_sub<-epimcmc(object = sir_mod0, tmin = 1, tmax = sub_time,
                    niter = 50000, sus.par.ini = 1, beta.ini = 1,
                    Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                    prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                    prior.sus.par = alphapar2, prior.beta.par = betapar2,
                    adapt = FALSE, acc.rate = NULL)
  
  mcmc_sub_sum<-summary(mcmc_sub)
  ##extract the parameter estimates for the subset data
  alpha_sub<-mcmc_sub_sum$statistics[1,1]
  alpha_sd_sub<-mcmc_sub_sum$statistics[1,2]
  beta_sub<-mcmc_sub_sum$statistics[2,1]
  beta_sd_sub<-mcmc_sub_sum$statistics[2,2]
  ##calculate absolute biases 
  alpha_abs_bias<-abs(alpha_100p-alpha_sub)
  alpha_sd_abs_bias<-abs(alpha_sd_100p-alpha_sd_sub)
  beta_abs_bias<-abs(beta_100p-beta_sub)
  beta_sd_abs_bias<-abs(beta_sd_100p-beta_sd_sub)
  result<-c(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias,
            alpha_sub,alpha_sd_sub,beta_sub,beta_sd_sub)
  return(result)
}


##setting up for parallel computing 
library(foreach)
library(doParallel)
##detecting the number of cores
##number of cores must match up the number of CPUs requested in your slurm script 
cores=detectCores()
cl <- makeCluster(8) #not to overload your computer
registerDoParallel(cl)


sim<-100
sub_time<-2
n<-20
##running foreach loop via parallel computing
result<-foreach(i=1:n) %dopar%{
  library(EpiILM)
  res<-temp_sub(sim,sub_time)
  alpha_abs_bias<-res[1]
  alpha_sd_abs_bias<-res[2]
  beta_abs_bias<-res[3]
  beta_sd_abs_bias <-res[4]
  alpha_sub<-res[5]
  alpha_sd_sub<-res[6]
  beta_sub<-res[7]
  beta_sd_sub<-res[8]
  cbind(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias,
        alpha_sub,alpha_sd_sub,beta_sub,beta_sd_sub)
}

##extracting the absolute biases in dataframe format
alpha_abs_bias<-c(rep(NA,n))
beta_abs_bias<-c(rep(NA,n))
alpha_sd_abs_bias<-c(rep(NA,n))
beta_sd_abs_bias<-c(rep(NA,n))
alpha_sub<-c(rep(NA,n))
alpha_sd_sub<-c(rep(NA,n))
beta_sub<-c(rep(NA,n))
beta_sd_sub<-c(rep(NA,n))

for(i in 1:n){
  out<-result[[i]]
  alpha_abs_bias[i]<-out[1]
  alpha_sd_abs_bias[i]<-out[2]
  beta_abs_bias[i]<-out[3]
  beta_sd_abs_bias[i]<-out[4]
  alpha_sub[i]<-out[5]
  alpha_sd_sub[i]<-out[6]
  beta_sub[i]<-out[7]
  beta_sd_sub[i]<-out[8]
}


##data frame for parameter estiamtes 
df_estimates=data.frame(alpha_sub,alpha_sd_sub,beta_sub,beta_sd_sub)
write.csv(df_estimates,"/home/thethtetchan.nyein/subt2_100_estimates.csv", row.names = FALSE)

##dataframe for absolute biases
df_abs_biases=data.frame(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
write.csv(df_abs_biases,"/home/thethtetchan.nyein/subt2_100_abs_biases.csv", row.names = FALSE)




