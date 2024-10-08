##procedures for the code
##generate epidemic data positioned on the \sqrt n x \sqrt n square grid
##apply the metropolis-hastings random walk mcmc
##if the resulting the parameter estimates fail the gweke's diagnostic test or fail to reach the stationary distribution,
##generate new epidemics until the stationary distribution is reached
##generate coordinates for the subsets based on the methods
##create subset areas
##apply MCMC to subsets
## find the overall means of the sample and find the absolute difference between sample mean - mean from the overall data
####repeat the same procedure for twenty times
##find the mean and sd of the absolute biases

library(dplyr)
library(reshape2)
library(purrr)
library(EpiILM)
library(ggplot2)

avg_center<-function(sam,pv){
  
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
  
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_sd_100p<-sum_100p_spatial$statistics[1,2]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_sd_100p<-sum_100p_spatial$statistics[2,2]
  
  
  alpha_sub<-c(rep(NA,length(pv)))
  alpha_sd_sub<-c(rep(NA,length(pv)))
  beta_sub<-c(rep(NA,length(pv)))
  beta_sd_sub<-c(rep(NA,length(pv)))
  
  
  for(k in 1:length(pv)){
    percent=pv[k]
    ## take multiple subset and take the average of the estimates
    ## for a certain percentage of data
    xl_p<-(sqrt(sam)-sqrt(sam*percent))/2
    xu_p<-(sqrt(sam)+sqrt(sam*percent))/2
    yl_p<-(sqrt(sam)-sqrt(sam*percent))/2
    yu_p<-(sqrt(sam)+sqrt(sam*percent))/2
    
    #extracting x_coordinates
    ind_vec_x_p<-c()
    
    for(i in 1:nrow(sir_mod0$XYcoordinates)){
      test_vec_p<-c(sir_mod0$XYcoordinates[i,])
      if(test_vec_p[1]>xl_p){
        if(test_vec_p[1]<xu_p){
          ind_vec_x_p<-c(ind_vec_x_p,i)
        }
      }
    }
    
    ind_vec_y_p<-c()
    for(i in 1:nrow(sir_mod0$XYcoordinates)){
      test_vec_p<-c(sir_mod0$XYcoordinates[i,])
      if(test_vec_p[2]>yl_p){
        if(test_vec_p[2]<yu_p){
          ind_vec_y_p<-c(ind_vec_y_p,i)
        }
      }
    }
    
    x_final_ind_p<-ind_vec_x_p[which(ind_vec_x_p %in% ind_vec_y_p)]
    XY_p_dat<-sir_mod0$XYcoordinates[x_final_ind_p,]
    x_vec_p<-XY_p_dat[,1]
    y_vec_p<-XY_p_dat[,2]
    inftime_p<-sir_mod0[["inftime"]][x_final_ind_p]
    remtime_p<-sir_mod0[["remtime"]][x_final_ind_p]
    lambda_p<-rep(3,length(x_final_ind_p))
    sir_p<-as.epidata(type = "SIR", n = length(x_vec_p),x = x_vec_p, y = y_vec_p, 
                      inftime =inftime_p,infperiod = lambda_p)
    mcmc_outp<-epimcmc(object = sir_p, tmin = min(inftime_p), tmax = max(inftime_p),
                       niter = 50000, sus.par.ini = 1, beta.ini = 1,
                       Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                       prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                       prior.sus.par = alphapar2, prior.beta.par = betapar2,
                       adapt = FALSE, acc.rate = NULL)
    sum_p_spatial<-summary(mcmc_outp, start = 10001)
    
    alpha_sub[k]<-sum_p_spatial$statistics[1,1]
    
    alpha_sd_sub[k]<-sum_p_spatial$statistics[1,2]
    
    
    beta_sub[k]<-sum_p_spatial$statistics[2,1]
    
    beta_sd_sub[k]<-sum_p_spatial$statistics[2,2]
  }
  
  alpha_abs_bias<-abs(mean(alpha_sub)-alpha_100p)
  alpha_sd_abs_bias<-abs(mean(alpha_sd_sub)-alpha_sd_100p)
  beta_abs_bias<-abs(mean(beta_sub)-beta_100p)
  beta_sd_abs_bias<-abs(mean(beta_sd_sub)-beta_sd_100p)
  result<-c(alpha_abs_bias,alpha_sd_abs_bias,
            beta_abs_bias,beta_sd_abs_bias)
  return(result)
}


library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(10) #not to overload your computer
registerDoParallel(cl)
n=20
sam=100
pv=c(0.25,0.36,0.49,0.64)

result<-foreach(i=1:n) %dopar%{
  library(EpiILM)
  res<-avg_center(sam,pv)
  alpha_abs_bias<-res[1]
  alpha_sd_abs_bias<-res[2]
  beta_abs_bias<-res[3]
  beta_sd_abs_bias <-res[4]
  cbind(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
}

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
write.csv(df_abs_biases,"/home/thethtetchan.nyein/avg_center_100.csv", row.names = FALSE)




