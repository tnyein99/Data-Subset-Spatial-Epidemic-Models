library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(EpiILM)

##function for subset and finding absolute biases 
##input=sample size (n), subset percentage (percent)
##output= absolute biases of all parameter estimates
random_subset<-function(sam,percent){
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
  
  
  ## take multiple subset and take the average of the estimates (randomization factor included)
  ## for a certain percentage of data
  repeat{
    xl_p1<-runif(1,0,sqrt(sam))
    xu_p1<-xl_p1+sqrt(sam*percent)
    ##making sure upper limit stays within the square grid 
    if(xu_p1<=sqrt(sam))
      break
  }
  
  repeat{
    yl_p1<-runif(1,0,sqrt(sam))
    yu_p1<-yl_p1+sqrt(sam*percent)
    ##making sure upper limit stays within the square grid 
    if(yu_p1<=sqrt(sam))
      break
  }
  
  
  
  #extracting x_coordinates
  ind_vec_x_p1<-c()
  #extracting indices of the x coordinates which are less than xl_p1 and xu_p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p1[1]>xl_p1){
      if(test_vec_p1[1]<xu_p1){
        ind_vec_x_p1<-c(ind_vec_x_p1,i)
      }
    }
  }
  
  ind_vec_y_p1<-c()
  #extracting indices of the x coordinates which are less than yl_p1 and yu_p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p1[2]>yl_p1){
      if(test_vec_p1[2]<yu_p1){
        ind_vec_y_p1<-c(ind_vec_y_p1,i)
      }
    }
  }
  ##extracting the subsets which include in both x and y subset
  x_final_ind_p1<-ind_vec_x_p1[which(ind_vec_x_p1 %in% ind_vec_y_p1)]
  XY_p1_dat<-sir_mod0$XYcoordinates[x_final_ind_p1,]
  x_vec_p1<-XY_p1_dat[,1]
  y_vec_p1<-XY_p1_dat[,2]
  inftime_p1<-sir_mod0[["inftime"]][x_final_ind_p1]
  remtime_p1<-sir_mod0[["remtime"]][x_final_ind_p1]
  lambda_p1<-rep(3,length(x_final_ind_p1))
  sir_p1<-as.epidata(type = "SIR", n = length(x_vec_p1),x = x_vec_p1, y = y_vec_p1, 
                     inftime =inftime_p1,infperiod = lambda_p1)
  mcmc_outp1<-epimcmc(object = sir_p1, tmin = min(inftime_p1), tmax = max(inftime_p1),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_p1_spatial<-summary(mcmc_outp1, start = 10001)
  
  
  ## subset 2
  
  repeat{
    xl_p2<-runif(1,0,sqrt(sam))
    xu_p2<-xl_p2+sqrt(sam*percent)
    
    if(xu_p2<=sqrt(sam))
      break
  }
  
  
  repeat{
    yl_p2<-runif(1,0,sqrt(sam))
    yu_p2<-yl_p2+sqrt(sam*percent)
    if(yu_p2<=sqrt(sam))
      break
  }
  #extracting x_coordinates
  ind_vec_x_p2<-c()
  #extracting indices of the x coordinates which are less than xl_p2 and xu_p2
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p2[1]>xl_p2){
      if(test_vec_p2[1]<xu_p2){
        ind_vec_x_p2<-c(ind_vec_x_p2,i)
      }
    }
  }
  #extracting indices of the y coordinates which are less than yl_p2 and yu_p2
  ind_vec_y_p2<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p2[2]>yl_p2){
      if(test_vec_p2[2]<yu_p2){
        ind_vec_y_p2<-c(ind_vec_y_p2,i)
      }
    }
  }
  
  x_final_ind_p2<-ind_vec_x_p2[which(ind_vec_x_p2 %in% ind_vec_y_p2)]
  XY_p2_dat<-sir_mod0$XYcoordinates[x_final_ind_p2,]
  x_vec_p2<-XY_p2_dat[,1]
  y_vec_p2<-XY_p2_dat[,2]
  inftime_p2<-sir_mod0[["inftime"]][x_final_ind_p2]
  remtime_p2<-sir_mod0[["remtime"]][x_final_ind_p2]
  lambda_p2<-rep(3,length(x_final_ind_p2))
  sir_p2<-as.epidata(type = "SIR", n = length(x_vec_p2),x = x_vec_p2, y = y_vec_p2, 
                     inftime =inftime_p2,infperiod = lambda_p2)
  mcmc_outp2<-epimcmc(object = sir_p2, tmin = min(inftime_p2), tmax = max(inftime_p2),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_p2_spatial<-summary(mcmc_outp2, start = 10001)
  
  repeat{
    xl_p3<-runif(1,0,sqrt(sam))
    xu_p3<-xl_p3+sqrt(sam*percent)
    
    if(xu_p3<=sqrt(sam))
      break
  }
  
  
  repeat{
    yl_p3<-runif(1,0,sqrt(sam))
    yu_p3<-yl_p3+sqrt(sam*percent)
    if(yu_p3<=sqrt(sam))
      break
  }
  
  
  #extracting x_coordinates
  ind_vec_x_p3<-c()
  #extracting indices of the x coordinates which are less than xl_p3 and xu_p3
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p3[1]>xl_p3){
      if(test_vec_p3[1]<xu_p3){
        ind_vec_x_p3<-c(ind_vec_x_p3,i)
      }
    }
  }
  #extracting indices of the y coordinates which are less than yl_p3 and yu_p3
  ind_vec_y_p3<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p3[2]>yl_p3){
      if(test_vec_p3[2]<yu_p3){
        ind_vec_y_p3<-c(ind_vec_y_p3,i)
      }
    }
  }
  
  x_final_ind_p3<-ind_vec_x_p3[which(ind_vec_x_p3 %in% ind_vec_y_p3)]
  XY_p3_dat<-sir_mod0$XYcoordinates[x_final_ind_p3,]
  x_vec_p3<-XY_p3_dat[,1]
  y_vec_p3<-XY_p3_dat[,2]
  inftime_p3<-sir_mod0[["inftime"]][x_final_ind_p3]
  remtime_p3<-sir_mod0[["remtime"]][x_final_ind_p3]
  lambda_p3<-rep(3,length(x_final_ind_p3))
  sir_p3<-as.epidata(type = "SIR", n = length(x_vec_p3),x = x_vec_p3, y = y_vec_p3, 
                     inftime =inftime_p3,infperiod = lambda_p3)
  mcmc_outp3<-epimcmc(object = sir_p3, tmin = min(inftime_p3), tmax = max(inftime_p3),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_p3_spatial<-summary(mcmc_outp3, start = 10001)
  
  repeat{
    xl_p4<-runif(1,0,sqrt(sam))
    xu_p4<-xl_p4+sqrt(sam*percent)
    if(xu_p4<=sqrt(sam))
      break
  }
  
  repeat{
    yl_p4<-runif(1,0,sqrt(sam))
    yu_p4<-yl_p4+sqrt(sam*percent)
    if(yu_p4<=sqrt(sam))
      break
  }
  
  
  #extracting x_coordinates
  ind_vec_x_p4<-c()
  #extracting indices of the x coordinates which are less than xl_p4 and xu_p4
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p4[1]>xl_p4){
      if(test_vec_p4[1]<xu_p4){
        ind_vec_x_p4<-c(ind_vec_x_p4,i)
      }
    }
  }
  #extracting indices of the y coordinates which are less than yl_p4 and yu_p4
  ind_vec_y_p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_p4[2]>yl_p4){
      if(test_vec_p4[2]<yu_p4){
        ind_vec_y_p4<-c(ind_vec_y_p4,i)
      }
    }
  }
  
  x_final_ind_p4<-ind_vec_x_p4[which(ind_vec_x_p4 %in% ind_vec_y_p4)]
  XY_p4_dat<-sir_mod0$XYcoordinates[x_final_ind_p4,]
  x_vec_p4<-XY_p4_dat[,1]
  y_vec_p4<-XY_p4_dat[,2]
  inftime_p4<-sir_mod0[["inftime"]][x_final_ind_p4]
  remtime_p4<-sir_mod0[["remtime"]][x_final_ind_p4]
  lambda_p4<-rep(3,length(x_final_ind_p4))
  sir_p4<-as.epidata(type = "SIR", n = length(x_vec_p4),x = x_vec_p4, y = y_vec_p4, 
                     inftime =inftime_p4,infperiod = lambda_p4)
  mcmc_outp4<-epimcmc(object = sir_p4, tmin = min(inftime_p4), tmax = max(inftime_p4),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_p4_spatial<-summary(mcmc_outp4, start = 10001)
  
  ##extracting the postetior parameter estimates of entire epidemic data
  alpha_p<-mean(c(sum_p1_spatial$statistics[1,1],sum_p2_spatial$statistics[1,1],
                  sum_p3_spatial$statistics[1,1],sum_p4_spatial$statistics[1,1]))
  
  alpha_sd_p<-mean(c(sum_p1_spatial$statistics[1,2],sum_p2_spatial$statistics[1,2],
                     sum_p3_spatial$statistics[1,2],sum_p4_spatial$statistics[1,2]))
  
  beta_p<-mean(c(sum_p1_spatial$statistics[2,1],sum_p2_spatial$statistics[2,1],
                 sum_p3_spatial$statistics[2,1],sum_p4_spatial$statistics[2,1]))
  
  beta_sd_p<-mean(c(sum_p1_spatial$statistics[2,2],sum_p2_spatial$statistics[2,2],
                    sum_p3_spatial$statistics[2,2],sum_p4_spatial$statistics[2,2]))
  
  
  ## calculating absolute biases 
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_sd_100p<-sum_100p_spatial$statistics[1,2]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_sd_100p<-sum_100p_spatial$statistics[2,2]
  ## calculating absolute biases  
  alpha_abs_bias<-abs(alpha_p-alpha_100p)
  alpha_sd_abs_bias<-abs(alpha_sd_p-alpha_sd_100p)
  beta_abs_bias<-abs(beta_p-beta_100p)
  beta_sd_abs_bias<-abs(beta_sd_p-beta_sd_100p)
  ##create result vector
  
  result<-c(alpha_abs_bias,alpha_sd_abs_bias,beta_abs_bias,beta_sd_abs_bias)
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


##running foreach loop via parallel computing
sim<-100
n<-20
per<-0.1

result<-foreach(i=1:n) %dopar%{
  library(EpiILM)
  res<-random_subset(sim,per)
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
write.csv(df_abs_biases,"/home/thethtetchan.nyein/percent10_random_100_abs_biases.csv", row.names = FALSE)


