library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(EpiILM)

##input: population size or n (sam)
##output: absolute biases 
temporal_median<-function(sam){
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
  
  ##subset 1: minimum infection time to 3
  ##subset 2: 4 to 7
  ##subset 3: 8 to 11
  ##subset 4: 12 to 15
  mcmc_out_t1<-epimcmc(object = sir_mod0, tmin = min(sir_mod0$inftime), tmax = 3,
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_t1_spatial<-summary(mcmc_out_t1,start=10001)
  
  
  
  mcmc_out_t2<-epimcmc(object = sir_mod0, tmin = 4, tmax =7,
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_t2_spatial<-summary(mcmc_out_t2,start=10001)
  
  
  mcmc_out_t3<-epimcmc(object = sir_mod0, tmin = 8, tmax = 11,
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_t3_spatial<-summary(mcmc_out_t3,start=10001)
  
  
  mcmc_out_t4<-epimcmc(object = sir_mod0, tmin = 12, tmax = 15,
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_t4_spatial<-summary(mcmc_out_t4,start=10001)
  
  ##find median by parameter estimates
  alpha_t_median<-median(c(sum_t1_spatial$statistics[1,1],sum_t2_spatial$statistics[1,1],
                       sum_t3_spatial$statistics[1,1],sum_t4_spatial$statistics[1,1]))
  
  alpha_sd_t_median<-median(c(sum_t1_spatial$statistics[1,2],sum_t2_spatial$statistics[1,2],
                          sum_t3_spatial$statistics[1,2],sum_t4_spatial$statistics[1,2]))
  
  beta_t_median<-median(c(sum_t1_spatial$statistics[2,1],sum_t2_spatial$statistics[2,1],
                      sum_t3_spatial$statistics[2,1],sum_t4_spatial$statistics[2,1]))
  
  beta_sd_t_median<-median(c(sum_t1_spatial$statistics[2,2],sum_t2_spatial$statistics[2,2],
                         sum_t3_spatial$statistics[2,2],sum_t4_spatial$statistics[2,2]))
  
  
  #extract the parameter estimates for the subset data
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_sd_100p<-sum_100p_spatial$statistics[1,2]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_sd_100p<-sum_100p_spatial$statistics[2,2]
  #find the absolute bias 
  alpha_abs_bias<-abs(alpha_t_median-alpha_100p)
  alpha_sd_abs_bias<-abs(alpha_sd_t_median-alpha_sd_100p)
  beta_abs_bias<-abs(beta_t_median-beta_100p)
  beta_sd_abs_bias<-abs(beta_sd_t_median-beta_sd_100p)
  ##create result vector
  result<-c(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
  return(result)
}

#setting up parallel computing 
library(foreach)
library(doParallel)
##detecting the number of cores
##number of cores must match up the number of CPUs requested in your slurm script 
cores=detectCores()
cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)

sim<-100
n<-20

result<-foreach(i=1:n) %dopar%{
  library(EpiILM)
  res<-temporal_median(sim)
  alpha_abs_bias<-res[1]
  alpha_sd_abs_bias<-res[2]
  beta_abs_bias<-res[3]
  beta_sd_abs_bias <-res[4]
  cbind(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
  
}


##extracting the absolute biases in dataframe format
alpha_abs_bias<-c(rep(NA,n))
alpha_sd_abs_bias<-c(rep(NA,n))
beta_abs_bias<-c(rep(NA,n))
beta_sd_abs_bias<-c(rep(NA,n))


for(i in 1:n){
  out<-data.frame(result[[i]])
  alpha_abs_bias[i]<-out$alpha_abs_bias
  alpha_sd_abs_bias[i]<-out$alpha_sd_abs_bias
  beta_abs_bias[i]<-out$beta_abs_bias
  beta_sd_abs_bias[i]<-out$beta_sd_abs_bias
}


##dataframe for absolute biases
df_abs_biases=data.frame(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
write.csv(df_abs_biases,"/home/thethtetchan.nyein/temporal_median4_100_abs_biases.csv", row.names = FALSE)







