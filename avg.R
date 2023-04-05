library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(EpiILM)

##function for subset and finding absolute biases 
##input=sample size (n)
##output= absolute biases of all parameter estimates
avg<-function(sam){
  x<- runif(sam,0,sqrt(sam))
  y<- runif(sam,0,sqrt(sam))
  lambda <- rep(3, sam)
  alphapar2 <- c(1, 1)
  betapar2 <- c(1, 1)
  
  
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
  #for 25 percent of the data
  xl_25p<-2
  xu_25p<-2+sqrt(sam*0.25)
  yl_25p<-2
  yu_25p<-2+sqrt(sam*0.25)
  
  #extracting x_coordinates
  ind_vec_x_25p<-c()
  #extracting indices of the x coordinates which are less than xl_25p and xu_25p
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p[1]>xl_25p){
      if(test_vec_25p[1]<xu_25p){
        ind_vec_x_25p<-c(ind_vec_x_25p,i)
      }
    }
  }
  #extracting indices of the y coordinates which are less than yl_25p and yu_25p
  ind_vec_y_25p<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p[2]>yl_25p){
      if(test_vec_25p[2]<yu_25p){
        ind_vec_y_25p<-c(ind_vec_y_25p,i)
      }
    }
  }
  
  x_final_ind_25p<-ind_vec_x_25p[which(ind_vec_x_25p %in% ind_vec_y_25p)]
  XY_25p_dat<-sir_mod0$XYcoordinates[x_final_ind_25p,]
  x_vec_25p<-XY_25p_dat[,1]
  y_vec_25p<-XY_25p_dat[,2]
  inftime_25p<-sir_mod0[["inftime"]][x_final_ind_25p]
  remtime_25p<-sir_mod0[["remtime"]][x_final_ind_25p]
  lambda_25p<-rep(3,length(x_final_ind_25p))
  sir_25p<-as.epidata(type = "SIR", n = length(x_vec_25p),x = x_vec_25p, y = y_vec_25p, 
                      inftime =inftime_25p,infperiod = lambda_25p)
  mcmc_out25p<-epimcmc(object = sir_25p, tmin = min(inftime_25p), tmax = max(inftime_25p),
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_25p_spatial<-summary(mcmc_out25p, start = 10001)
  #for 36 percent of data
  
  xl_36p<-1
  xu_36p<-1+sqrt(sam*0.36)
  
  yl_36p<-1
  yu_36p<-1+sqrt(sam*0.36)
  
  #extracting x_coordinates
  ind_vec_x_36p<-c()
  #extracting indices of the x coordinates which are less than xl_36p and xu_36p
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p[1]>xl_36p){
      if(test_vec_36p[1]<xu_36p){
        ind_vec_x_36p<-c(ind_vec_x_36p,i)
      }
    }
  }
  
  
  #extracting indices of the y coordinates which are less than yl_36p and yu_36p
  
  ind_vec_y_36p<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p[2]>yl_36p){
      if(test_vec_36p[2]<yu_36p){
        ind_vec_y_36p<-c(ind_vec_y_36p,i)
      }
    }
  }
  
  x_final_ind_36p<-ind_vec_x_36p[which(ind_vec_x_36p %in% ind_vec_y_36p)]
  XY_36p_dat<-sir_mod0$XYcoordinates[x_final_ind_36p,]
  x_vec_36p<-XY_36p_dat[,1]
  y_vec_36p<-XY_36p_dat[,2]
  
  inftime_36p<-sir_mod0[["inftime"]][x_final_ind_36p]
  remtime_36p<-sir_mod0[["remtime"]][x_final_ind_36p]
  lambda_36p<-rep(3,length(x_final_ind_36p))
  sir_36p<-as.epidata(type = "SIR", n = length(x_vec_36p),x = x_vec_36p, y = y_vec_36p, 
                      inftime =inftime_36p,infperiod = lambda_36p)
  
  mcmc_out36p<-epimcmc(object = sir_36p, tmin = min(inftime_36p), tmax = max(inftime_36p),
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  
  sum_36p_spatial<-summary(mcmc_out36p, start = 10001)
  
  
  
  #for 49 percent of the data
  xl_49p<-1.5
  xu_49p<-1.5+sqrt(sam*0.49)
  
  yl_49p<-1.5
  yu_49p<-1.5+sqrt(sam*0.49)
  
  
  
  
  #extracting x_coordinates
  ind_vec_x_49p<-c()
  #extracting indices of the x coordinates which are less than xl_49p and xu_49p
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p[1]>xl_49p){
      if(test_vec_49p[1]<xu_49p){
        ind_vec_x_49p<-c(ind_vec_x_49p,i)
      }
    }
  }
  
  
  #extracting indices of the y coordinates which are less than yl_49p and yu_49p
  
  ind_vec_y_49p<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p[2]>yl_49p){
      if(test_vec_49p[2]<yu_49p){
        ind_vec_y_49p<-c(ind_vec_y_49p,i)
      }
    }
  }
  
  x_final_ind_49p<-ind_vec_x_49p[which(ind_vec_x_49p %in% ind_vec_y_49p)]
  XY_49p_dat<-sir_mod0$XYcoordinates[x_final_ind_49p,]
  x_vec_49p<-XY_49p_dat[,1]
  y_vec_49p<-XY_49p_dat[,2]
  
  inftime_49p<-sir_mod0[["inftime"]][x_final_ind_49p]
  remtime_49p<-sir_mod0[["remtime"]][x_final_ind_49p]
  lambda_49p<-rep(3,length(x_final_ind_49p))
  sir_49p<-as.epidata(type = "SIR", n = length(x_vec_49p),x = x_vec_49p, y = y_vec_49p, 
                      inftime =inftime_49p,infperiod = lambda_49p)
  
  mcmc_out49p<-epimcmc(object = sir_49p, tmin = min(inftime_49p), tmax = max(inftime_49p),
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  
  sum_49p_spatial<-summary(mcmc_out49p, start = 10001)
  
  
  #for 64 percent of the data
  xl_64p<-1
  xu_64p<-1+sqrt(sam*0.64)
  
  yl_64p<-1
  yu_64p<-1+sqrt(sam*0.64)
  
  
  #extracting x_coordinates
  ind_vec_x_64p<-c()
  #extracting indices of the x coordinates which are less than xl_64p and xu_64p
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p[1]>xl_64p){
      if(test_vec_64p[1]<xu_64p){
        ind_vec_x_64p<-c(ind_vec_x_64p,i)
      }
    }
  }
  
  #extracting indices of the y coordinates which are less than yl_64p and yu_64p
  
  ind_vec_y_64p<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p[2]>yl_64p){
      if(test_vec_64p[2]<yu_64p){
        ind_vec_y_64p<-c(ind_vec_y_64p,i)
      }
    }
  }
  
  x_final_ind_64p<-ind_vec_x_64p[which(ind_vec_x_64p %in% ind_vec_y_64p)]
  XY_64p_dat<-sir_mod0$XYcoordinates[x_final_ind_64p,]
  x_vec_64p<-XY_64p_dat[,1]
  y_vec_64p<-XY_64p_dat[,2]
  
  inftime_64p<-sir_mod0[["inftime"]][x_final_ind_64p]
  remtime_64p<-sir_mod0[["remtime"]][x_final_ind_64p]
  lambda_64p<-rep(3,length(x_final_ind_64p))
  sir_64p<-as.epidata(type = "SIR", n = length(x_vec_64p),x = x_vec_64p, y = y_vec_64p, 
                      inftime =inftime_64p,infperiod = lambda_64p)
  
  mcmc_out64p<-epimcmc(object = sir_64p, tmin = min(inftime_64p), tmax = max(inftime_64p),
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
  sum_64p_spatial<-summary(mcmc_out64p, start = 10001)
  
  
  ##extracting the parameter estimate by each subset percentage
  ##taking the mean of those parameter estimates and find the bias 
  alpha_25p<-sum_25p_spatial$statistics[1,1]
  alpha_36p<-sum_36p_spatial$statistics[1,1]
  alpha_49p<-sum_49p_spatial$statistics[1,1]
  alpha_64p<-sum_64p_spatial$statistics[1,1]
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_est<-mean(c(alpha_25p,alpha_36p,alpha_49p,alpha_64p))
  alpha_abs_bias<-abs(alpha_est-alpha_100p)
  
  
  alpha_25p_sd<-sum_25p_spatial$statistics[1,2]
  alpha_36p_sd<-sum_36p_spatial$statistics[1,2]
  alpha_49p_sd<-sum_49p_spatial$statistics[1,2]
  alpha_64p_sd<-sum_64p_spatial$statistics[1,2]
  alpha_100p_sd<-sum_100p_spatial$statistics[1,2]
  alpha_sd_est<-mean(c(alpha_25p_sd,alpha_36p_sd,alpha_49p_sd,alpha_64p_sd))
  alpha_sd_abs_bias<-abs(alpha_sd_est-alpha_100p_sd)
  
  
  beta_25p<-sum_25p_spatial$statistics[2,1]
  beta_36p<-sum_36p_spatial$statistics[2,1]
  beta_49p<-sum_49p_spatial$statistics[2,1]
  beta_64p<-sum_64p_spatial$statistics[2,1]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_est<-mean(c(beta_25p,beta_36p,beta_49p,beta_64p))
  beta_abs_bias<-abs(beta_est-beta_100p)
  
  
  beta_25p_sd<-sum_25p_spatial$statistics[2,2]
  beta_36p_sd<-sum_36p_spatial$statistics[2,2]
  beta_49p_sd<-sum_49p_spatial$statistics[2,2]
  beta_64p_sd<-sum_64p_spatial$statistics[2,2]
  beta_100p_sd<-sum_100p_spatial$statistics[2,2]
  beta_sd_est<-mean(c(beta_25p_sd,beta_36p_sd,beta_49p_sd,beta_64p_sd))
  beta_sd_abs_bias<-abs(beta_sd_est-beta_100p_sd)
  
  result<-c(rep(NA,4))
  result[1]<-alpha_abs_bias
  result[2]<-alpha_sd_abs_bias
  result[3]<-beta_abs_bias
  result[4]<-beta_sd_abs_bias
  return(result)
}

##setting up for parallel computing 
library(foreach)
library(doParallel)
##detecting the number of cores
##number of cores must match up the number of CPUs requested in your slurm script 
cores=detectCores()
cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)

##running foreach loop via parallel computing
sim<-100
n<-20

result<-foreach(i=1:n) %dopar%{
  library(EpiILM)
  res<-avg(sim)
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

##extracting the absolute biases in dataframe format
for(i in 1:n){
  out<-data.frame(result[[i]])
  alpha_abs_bias[i]<-out$alpha_abs_bias
  alpha_sd_abs_bias[i]<-out$alpha_sd_abs_bias
  beta_abs_bias[i]<-out$beta_abs_bias
  beta_sd_abs_bias[i]<-out$beta_sd_abs_bias
}


##dataframe for absolute biases
df_abs_biases=data.frame(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
write.csv(df_abs_biases,"/home/thethtetchan.nyein/avg_100_abs_biases.csv", row.names = FALSE)








