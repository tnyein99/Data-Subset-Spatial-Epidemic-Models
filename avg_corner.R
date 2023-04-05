library(dplyr)
library(reshape2)
library(ggplot2)
library(purrr)
library(EpiILM)
##function for subset and finding absolute biases 
##input=sample size (n)
##output= absolute biases of all parameter estimates
avg_corner<-function(sam){
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
  
  
  ## take multiple subset and take the average of the estimates
  ## for a certain percentage of data
  xl_25p1<-0
  xu_25p1<-sqrt(sam*0.25)
  yl_25p1<-0
  yu_25p1<-sqrt(sam*0.25)
  
  #extracting x_coordinates
  ind_vec_x_25p1<-c()
  #extracting indices of the x coordinates which are between xl_25p1 and xu_25p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p1[1]>xl_25p1){
      if(test_vec_25p1[1]<xu_25p1){
        ind_vec_x_25p1<-c(ind_vec_x_25p1,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_25p1 and yu_25p1
  ind_vec_y_25p1<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p1[2]>yl_25p1){
      if(test_vec_25p1[2]<yu_25p1){
        ind_vec_y_25p1<-c(ind_vec_y_25p1,i)
      }
    }
  }
  
  x_final_ind_25p1<-ind_vec_x_25p1[which(ind_vec_x_25p1 %in% ind_vec_y_25p1)]
  XY_25p1_dat<-sir_mod0$XYcoordinates[x_final_ind_25p1,]
  x_vec_25p1<-XY_25p1_dat[,1]
  y_vec_25p1<-XY_25p1_dat[,2]
  inftime_25p1<-sir_mod0[["inftime"]][x_final_ind_25p1]
  remtime_25p1<-sir_mod0[["remtime"]][x_final_ind_25p1]
  lambda_25p1<-rep(3,length(x_final_ind_25p1))
  sir_25p1<-as.epidata(type = "SIR", n = length(x_vec_25p1),x = x_vec_25p1, y = y_vec_25p1, 
                       inftime =inftime_25p1,infperiod = lambda_25p1)
  mcmc_outp1<-epimcmc(object = sir_25p1, tmin = min(inftime_25p1), tmax = max(inftime_25p1),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_25p1_spatial<-summary(mcmc_outp1, start = 10001)
  
  
  ## subset 2
  
  xl_25p2<-sqrt(sam)-sqrt(sam*0.25)
  xu_25p2<-sqrt(sam)
  yl_25p2<-0
  yu_25p2<-sqrt(sam*0.25)
  
  #extracting x_coordinates
  ind_vec_x_25p2<-c()
  #extracting indices of the x coordinates which are between xl_25p2 and xu_25p2
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p2[1]>xl_25p2){
      if(test_vec_25p2[1]<xu_25p2){
        ind_vec_x_25p2<-c(ind_vec_x_25p2,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_25p2 and xu_25p2
  ind_vec_y_25p2<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p2[2]>yl_25p2){
      if(test_vec_25p2[2]<yu_25p2){
        ind_vec_y_25p2<-c(ind_vec_y_25p2,i)
      }
    }
  }
  
  x_final_ind_25p2<-ind_vec_x_25p2[which(ind_vec_x_25p2 %in% ind_vec_y_25p2)]
  XY_25p2_dat<-sir_mod0$XYcoordinates[x_final_ind_25p2,]
  x_vec_25p2<-XY_25p2_dat[,1]
  y_vec_25p2<-XY_25p2_dat[,2]
  inftime_25p2<-sir_mod0[["inftime"]][x_final_ind_25p2]
  remtime_25p2<-sir_mod0[["remtime"]][x_final_ind_25p2]
  lambda_25p2<-rep(3,length(x_final_ind_25p2))
  sir_25p2<-as.epidata(type = "SIR", n = length(x_vec_25p2),x = x_vec_25p2, y = y_vec_25p2, 
                       inftime =inftime_25p2,infperiod = lambda_25p2)
  mcmc_outp2<-epimcmc(object = sir_25p2, tmin = min(inftime_25p2), tmax = max(inftime_25p2),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_25p2_spatial<-summary(mcmc_outp2, start = 10001)
  
  xl_25p3<-0
  xu_25p3<-sqrt(sam*0.25)
  
  yl_25p3<-sqrt(sam)-sqrt(sam*0.25)
  yu_25p3<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_25p3<-c()
  #extracting indices of the x coordinates which are between xl_25p3 and xu_25p3
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p3[1]>xl_25p3){
      if(test_vec_25p3[1]<xu_25p3){
        ind_vec_x_25p3<-c(ind_vec_x_25p3,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_25p3 and yu_25p3
  ind_vec_y_25p3<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p3[2]>yl_25p3){
      if(test_vec_25p3[2]<yu_25p3){
        ind_vec_y_25p3<-c(ind_vec_y_25p3,i)
      }
    }
  }
  
  x_final_ind_25p3<-ind_vec_x_25p3[which(ind_vec_x_25p3 %in% ind_vec_y_25p3)]
  XY_25p3_dat<-sir_mod0$XYcoordinates[x_final_ind_25p3,]
  x_vec_25p3<-XY_25p3_dat[,1]
  y_vec_25p3<-XY_25p3_dat[,2]
  inftime_25p3<-sir_mod0[["inftime"]][x_final_ind_25p3]
  remtime_25p3<-sir_mod0[["remtime"]][x_final_ind_25p3]
  lambda_25p3<-rep(3,length(x_final_ind_25p3))
  sir_25p3<-as.epidata(type = "SIR", n = length(x_vec_25p3),x = x_vec_25p3, y = y_vec_25p3, 
                       inftime =inftime_25p3,infperiod = lambda_25p3)
  mcmc_outp3<-epimcmc(object = sir_25p3, tmin = min(inftime_25p3), tmax = max(inftime_25p3),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_25p3_spatial<-summary(mcmc_outp3, start = 10001)
  
  
  xl_25p4<-sqrt(sam)-sqrt(sam*0.25)
  xu_25p4<-sqrt(sam)
  
  yl_25p4<-sqrt(sam)-sqrt(sam*0.25)
  yu_25p4<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_25p4<-c()
  #extracting indices of the x coordinates which are between xl_25p4 and xu_25p4
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p4[1]>xl_25p4){
      if(test_vec_25p4[1]<xu_25p4){
        ind_vec_x_25p4<-c(ind_vec_x_25p4,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_25p4 and yu_25p4
  ind_vec_y_25p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_25p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_25p4[2]>yl_25p4){
      if(test_vec_25p4[2]<yu_25p4){
        ind_vec_y_25p4<-c(ind_vec_y_25p4,i)
      }
    }
  }
  
  x_final_ind_25p4<-ind_vec_x_25p4[which(ind_vec_x_25p4 %in% ind_vec_y_25p4)]
  XY_25p4_dat<-sir_mod0$XYcoordinates[x_final_ind_25p4,]
  x_vec_25p4<-XY_25p4_dat[,1]
  y_vec_25p4<-XY_25p4_dat[,2]
  inftime_25p4<-sir_mod0[["inftime"]][x_final_ind_25p4]
  remtime_25p4<-sir_mod0[["remtime"]][x_final_ind_25p4]
  lambda_25p4<-rep(3,length(x_final_ind_25p4))
  sir_25p4<-as.epidata(type = "SIR", n = length(x_vec_25p4),x = x_vec_25p4, y = y_vec_25p4, 
                       inftime =inftime_25p4,infperiod = lambda_25p4)
  mcmc_outp4<-epimcmc(object = sir_25p4, tmin = min(inftime_25p4), tmax = max(inftime_25p4),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_25p4_spatial<-summary(mcmc_outp4, start = 10001)
  
  ##extracting the parameter estimate by each subset percentage
  alpha_25p<-mean(c(sum_25p1_spatial$statistics[1,1],sum_25p2_spatial$statistics[1,1],
                    sum_25p3_spatial$statistics[1,1],sum_25p4_spatial$statistics[1,1]))
  
  alpha_sd_25p<-mean(c(sum_25p1_spatial$statistics[1,2],sum_25p2_spatial$statistics[1,2],
                       sum_25p3_spatial$statistics[1,2],sum_25p4_spatial$statistics[1,2]))
  
  beta_25p<-mean(c(sum_25p1_spatial$statistics[2,1],sum_25p2_spatial$statistics[2,1],
                   sum_25p3_spatial$statistics[2,1],sum_25p4_spatial$statistics[2,1]))
  
  beta_sd_25p<-mean(c(sum_25p1_spatial$statistics[2,2],sum_25p2_spatial$statistics[2,2],
                      sum_25p3_spatial$statistics[2,2],sum_25p4_spatial$statistics[2,2]))
  
  
  xl_36p1<-0
  xu_36p1<-sqrt(sam*0.36)
  yl_36p1<-0
  yu_36p1<-sqrt(sam*0.36)
  
  #extracting x_coordinates
  ind_vec_x_36p1<-c()
  #extracting indices of the x coordinates which are between xl_36p1 and xu_36p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p1[1]>xl_36p1){
      if(test_vec_36p1[1]<xu_36p1){
        ind_vec_x_36p1<-c(ind_vec_x_36p1,i)
      }
    }
  }
  
  ind_vec_y_36p1<-c()
  #extracting indices of the x coordinates which are between yl_36p1 and yu_36p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p1[2]>yl_36p1){
      if(test_vec_36p1[2]<yu_36p1){
        ind_vec_y_36p1<-c(ind_vec_y_36p1,i)
      }
    }
  }
  
  x_final_ind_36p1<-ind_vec_x_36p1[which(ind_vec_x_36p1 %in% ind_vec_y_36p1)]
  XY_36p1_dat<-sir_mod0$XYcoordinates[x_final_ind_36p1,]
  x_vec_36p1<-XY_36p1_dat[,1]
  y_vec_36p1<-XY_36p1_dat[,2]
  inftime_36p1<-sir_mod0[["inftime"]][x_final_ind_36p1]
  remtime_36p1<-sir_mod0[["remtime"]][x_final_ind_36p1]
  lambda_36p1<-rep(3,length(x_final_ind_36p1))
  sir_36p1<-as.epidata(type = "SIR", n = length(x_vec_36p1),x = x_vec_36p1, y = y_vec_36p1, 
                       inftime =inftime_36p1,infperiod = lambda_36p1)
  mcmc_outp1<-epimcmc(object = sir_36p1, tmin = min(inftime_36p1), tmax = max(inftime_36p1),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_36p1_spatial<-summary(mcmc_outp1, start = 10001)
  
  
  ## subset 2
  
  xl_36p2<-sqrt(sam)-sqrt(sam*0.36)
  xu_36p2<-sqrt(sam)
  yl_36p2<-0
  yu_36p2<-sqrt(sam*0.36)
  
  #extracting x_coordinates
  ind_vec_x_36p2<-c()
  #extracting indices of the x coordinates which are between xl_36p2 and xu_36p2
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p2[1]>xl_36p2){
      if(test_vec_36p2[1]<xu_36p2){
        ind_vec_x_36p2<-c(ind_vec_x_36p2,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_36p2 and yu_36p2
  ind_vec_y_36p2<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p2[2]>yl_36p2){
      if(test_vec_36p2[2]<yu_36p2){
        ind_vec_y_36p2<-c(ind_vec_y_36p2,i)
      }
    }
  }
  
  x_final_ind_36p2<-ind_vec_x_36p2[which(ind_vec_x_36p2 %in% ind_vec_y_36p2)]
  XY_36p2_dat<-sir_mod0$XYcoordinates[x_final_ind_36p2,]
  x_vec_36p2<-XY_36p2_dat[,1]
  y_vec_36p2<-XY_36p2_dat[,2]
  inftime_36p2<-sir_mod0[["inftime"]][x_final_ind_36p2]
  remtime_36p2<-sir_mod0[["remtime"]][x_final_ind_36p2]
  lambda_36p2<-rep(3,length(x_final_ind_36p2))
  sir_36p2<-as.epidata(type = "SIR", n = length(x_vec_36p2),x = x_vec_36p2, y = y_vec_36p2, 
                       inftime =inftime_36p2,infperiod = lambda_36p2)
  mcmc_outp2<-epimcmc(object = sir_36p2, tmin = min(inftime_36p2), tmax = max(inftime_36p2),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_36p2_spatial<-summary(mcmc_outp2, start = 10001)
  
  xl_36p3<-0
  xu_36p3<-sqrt(sam*0.36)
  
  yl_36p3<-sqrt(sam)-sqrt(sam*0.36)
  yu_36p3<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_36p3<-c()
  #extracting indices of the x coordinates which are between xl_36p3 and xu_36p3
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p3[1]>xl_36p3){
      if(test_vec_36p3[1]<xu_36p3){
        ind_vec_x_36p3<-c(ind_vec_x_36p3,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_36p3 and yu_36p3
  ind_vec_y_36p3<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p3[2]>yl_36p3){
      if(test_vec_36p3[2]<yu_36p3){
        ind_vec_y_36p3<-c(ind_vec_y_36p3,i)
      }
    }
  }
  
  x_final_ind_36p3<-ind_vec_x_36p3[which(ind_vec_x_36p3 %in% ind_vec_y_36p3)]
  XY_36p3_dat<-sir_mod0$XYcoordinates[x_final_ind_36p3,]
  x_vec_36p3<-XY_36p3_dat[,1]
  y_vec_36p3<-XY_36p3_dat[,2]
  inftime_36p3<-sir_mod0[["inftime"]][x_final_ind_36p3]
  remtime_36p3<-sir_mod0[["remtime"]][x_final_ind_36p3]
  lambda_36p3<-rep(3,length(x_final_ind_36p3))
  sir_36p3<-as.epidata(type = "SIR", n = length(x_vec_36p3),x = x_vec_36p3, y = y_vec_36p3, 
                       inftime =inftime_36p3,infperiod = lambda_36p3)
  mcmc_outp3<-epimcmc(object = sir_36p3, tmin = min(inftime_36p3), tmax = max(inftime_36p3),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_36p3_spatial<-summary(mcmc_outp3, start = 10001)
  
  
  xl_36p4<-sqrt(sam)-sqrt(sam*0.36)
  xu_36p4<-sqrt(sam)
  
  yl_36p4<-sqrt(sam)-sqrt(sam*0.36)
  yu_36p4<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_36p4<-c()
  #extracting indices of the x coordinates which are between xl_36p4 and xu_36p4
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p4[1]>xl_36p4){
      if(test_vec_36p4[1]<xu_36p4){
        ind_vec_x_36p4<-c(ind_vec_x_36p4,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_36p4 and yu_36p4
  ind_vec_y_36p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_36p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_36p4[2]>yl_36p4){
      if(test_vec_36p4[2]<yu_36p4){
        ind_vec_y_36p4<-c(ind_vec_y_36p4,i)
      }
    }
  }
  
  x_final_ind_36p4<-ind_vec_x_36p4[which(ind_vec_x_36p4 %in% ind_vec_y_36p4)]
  XY_36p4_dat<-sir_mod0$XYcoordinates[x_final_ind_36p4,]
  x_vec_36p4<-XY_36p4_dat[,1]
  y_vec_36p4<-XY_36p4_dat[,2]
  inftime_36p4<-sir_mod0[["inftime"]][x_final_ind_36p4]
  remtime_36p4<-sir_mod0[["remtime"]][x_final_ind_36p4]
  lambda_36p4<-rep(3,length(x_final_ind_36p4))
  sir_36p4<-as.epidata(type = "SIR", n = length(x_vec_36p4),x = x_vec_36p4, y = y_vec_36p4, 
                       inftime =inftime_36p4,infperiod = lambda_36p4)
  mcmc_outp4<-epimcmc(object = sir_36p4, tmin = min(inftime_36p4), tmax = max(inftime_36p4),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_36p4_spatial<-summary(mcmc_outp4, start = 10001)
  
  ##extracting the parameter estimate by each subset percentage
  alpha_36p<-mean(c(sum_36p1_spatial$statistics[1,1],sum_36p2_spatial$statistics[1,1],
                    sum_36p3_spatial$statistics[1,1],sum_36p4_spatial$statistics[1,1]))
  
  alpha_sd_36p<-mean(c(sum_36p1_spatial$statistics[1,2],sum_36p2_spatial$statistics[1,2],
                       sum_36p3_spatial$statistics[1,2],sum_36p4_spatial$statistics[1,2]))
  
  beta_36p<-mean(c(sum_36p1_spatial$statistics[2,1],sum_36p2_spatial$statistics[2,1],
                   sum_36p3_spatial$statistics[2,1],sum_36p4_spatial$statistics[2,1]))
  
  beta_sd_36p<-mean(c(sum_36p1_spatial$statistics[2,2],sum_36p2_spatial$statistics[2,2],
                      sum_36p3_spatial$statistics[2,2],sum_36p4_spatial$statistics[2,2]))
  
  
  xl_49p1<-0
  xu_49p1<-sqrt(sam*0.49)
  yl_49p1<-0
  yu_49p1<-sqrt(sam*0.49)
  
  #extracting x_coordinates
  ind_vec_x_49p1<-c()
  #extracting indices of the x coordinates which are between xl_49p1 and xu_49p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p1[1]>xl_49p1){
      if(test_vec_49p1[1]<xu_49p1){
        ind_vec_x_49p1<-c(ind_vec_x_49p1,i)
      }
    }
  }
  
  ind_vec_y_49p1<-c()
  #extracting indices of the x coordinates which are between yl_49p1 and yu_49p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p1[2]>yl_49p1){
      if(test_vec_49p1[2]<yu_49p1){
        ind_vec_y_49p1<-c(ind_vec_y_49p1,i)
      }
    }
  }
  
  x_final_ind_49p1<-ind_vec_x_49p1[which(ind_vec_x_49p1 %in% ind_vec_y_49p1)]
  XY_49p1_dat<-sir_mod0$XYcoordinates[x_final_ind_49p1,]
  x_vec_49p1<-XY_49p1_dat[,1]
  y_vec_49p1<-XY_49p1_dat[,2]
  inftime_49p1<-sir_mod0[["inftime"]][x_final_ind_49p1]
  remtime_49p1<-sir_mod0[["remtime"]][x_final_ind_49p1]
  lambda_49p1<-rep(3,length(x_final_ind_49p1))
  sir_49p1<-as.epidata(type = "SIR", n = length(x_vec_49p1),x = x_vec_49p1, y = y_vec_49p1, 
                       inftime =inftime_49p1,infperiod = lambda_49p1)
  mcmc_outp1<-epimcmc(object = sir_49p1, tmin = min(inftime_49p1), tmax = max(inftime_49p1),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_49p1_spatial<-summary(mcmc_outp1, start = 10001)
  
  
  ## subset 2
  
  xl_49p2<-sqrt(sam)-sqrt(sam*0.49)
  xu_49p2<-sqrt(sam)
  yl_49p2<-0
  yu_49p2<-sqrt(sam*0.49)
  
  #extracting x_coordinates
  ind_vec_x_49p2<-c()
  #extracting indices of the x coordinates which are between xl_49p2 and xu_49p2
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p2[1]>xl_49p2){
      if(test_vec_49p2[1]<xu_49p2){
        ind_vec_x_49p2<-c(ind_vec_x_49p2,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_49p2 and yu_49p2
  ind_vec_y_49p2<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p2[2]>yl_49p2){
      if(test_vec_49p2[2]<yu_49p2){
        ind_vec_y_49p2<-c(ind_vec_y_49p2,i)
      }
    }
  }
  
  x_final_ind_49p2<-ind_vec_x_49p2[which(ind_vec_x_49p2 %in% ind_vec_y_49p2)]
  XY_49p2_dat<-sir_mod0$XYcoordinates[x_final_ind_49p2,]
  x_vec_49p2<-XY_49p2_dat[,1]
  y_vec_49p2<-XY_49p2_dat[,2]
  inftime_49p2<-sir_mod0[["inftime"]][x_final_ind_49p2]
  remtime_49p2<-sir_mod0[["remtime"]][x_final_ind_49p2]
  lambda_49p2<-rep(3,length(x_final_ind_49p2))
  sir_49p2<-as.epidata(type = "SIR", n = length(x_vec_49p2),x = x_vec_49p2, y = y_vec_49p2, 
                       inftime =inftime_49p2,infperiod = lambda_49p2)
  mcmc_outp2<-epimcmc(object = sir_49p2, tmin = min(inftime_49p2), tmax = max(inftime_49p2),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_49p2_spatial<-summary(mcmc_outp2, start = 10001)
  
  xl_49p3<-0
  xu_49p3<-sqrt(sam*0.49)
  
  yl_49p3<-sqrt(sam)-sqrt(sam*0.49)
  yu_49p3<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_49p3<-c()
  #extracting indices of the x coordinates which are between xl_49p3 and xu_49p3
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p3[1]>xl_49p3){
      if(test_vec_49p3[1]<xu_49p3){
        ind_vec_x_49p3<-c(ind_vec_x_49p3,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_49p3 and yu_49p3
  ind_vec_y_49p3<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p3[2]>yl_49p3){
      if(test_vec_49p3[2]<yu_49p3){
        ind_vec_y_49p3<-c(ind_vec_y_49p3,i)
      }
    }
  }
  
  x_final_ind_49p3<-ind_vec_x_49p3[which(ind_vec_x_49p3 %in% ind_vec_y_49p3)]
  XY_49p3_dat<-sir_mod0$XYcoordinates[x_final_ind_49p3,]
  x_vec_49p3<-XY_49p3_dat[,1]
  y_vec_49p3<-XY_49p3_dat[,2]
  inftime_49p3<-sir_mod0[["inftime"]][x_final_ind_49p3]
  remtime_49p3<-sir_mod0[["remtime"]][x_final_ind_49p3]
  lambda_49p3<-rep(3,length(x_final_ind_49p3))
  sir_49p3<-as.epidata(type = "SIR", n = length(x_vec_49p3),x = x_vec_49p3, y = y_vec_49p3, 
                       inftime =inftime_49p3,infperiod = lambda_49p3)
  mcmc_outp3<-epimcmc(object = sir_49p3, tmin = min(inftime_49p3), tmax = max(inftime_49p3),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_49p3_spatial<-summary(mcmc_outp3, start = 10001)
  
  
  xl_49p4<-sqrt(sam)-sqrt(sam*0.49)
  xu_49p4<-sqrt(sam)
  
  yl_49p4<-sqrt(sam)-sqrt(sam*0.49)
  yu_49p4<-sqrt(sam)
  
  #extracting indices of the x coordinates which are between xl_49p4 and xu_49p4
  ind_vec_x_49p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p4[1]>xl_49p4){
      if(test_vec_49p4[1]<xu_49p4){
        ind_vec_x_49p4<-c(ind_vec_x_49p4,i)
      }
    }
  }
  
  #extracting indices of the x coordinates which are between yl_49p4 and yu_49p4
  ind_vec_y_49p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_49p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_49p4[2]>yl_49p4){
      if(test_vec_49p4[2]<yu_49p4){
        ind_vec_y_49p4<-c(ind_vec_y_49p4,i)
      }
    }
  }
  
  x_final_ind_49p4<-ind_vec_x_49p4[which(ind_vec_x_49p4 %in% ind_vec_y_49p4)]
  XY_49p4_dat<-sir_mod0$XYcoordinates[x_final_ind_49p4,]
  x_vec_49p4<-XY_49p4_dat[,1]
  y_vec_49p4<-XY_49p4_dat[,2]
  inftime_49p4<-sir_mod0[["inftime"]][x_final_ind_49p4]
  remtime_49p4<-sir_mod0[["remtime"]][x_final_ind_49p4]
  lambda_49p4<-rep(3,length(x_final_ind_49p4))
  sir_49p4<-as.epidata(type = "SIR", n = length(x_vec_49p4),x = x_vec_49p4, y = y_vec_49p4, 
                       inftime =inftime_49p4,infperiod = lambda_49p4)
  mcmc_outp4<-epimcmc(object = sir_49p4, tmin = min(inftime_49p4), tmax = max(inftime_49p4),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_49p4_spatial<-summary(mcmc_outp4, start = 10001)
  
  ##extracting the parameter estimate by each subset percentage
  alpha_49p<-mean(c(sum_49p1_spatial$statistics[1,1],sum_49p2_spatial$statistics[1,1],
                    sum_49p3_spatial$statistics[1,1],sum_49p4_spatial$statistics[1,1]))
  
  alpha_sd_49p<-mean(c(sum_49p1_spatial$statistics[1,2],sum_49p2_spatial$statistics[1,2],
                       sum_49p3_spatial$statistics[1,2],sum_49p4_spatial$statistics[1,2]))
  
  beta_49p<-mean(c(sum_49p1_spatial$statistics[2,1],sum_49p2_spatial$statistics[2,1],
                   sum_49p3_spatial$statistics[2,1],sum_49p4_spatial$statistics[2,1]))
  
  beta_sd_49p<-mean(c(sum_49p1_spatial$statistics[2,2],sum_49p2_spatial$statistics[2,2],
                      sum_49p3_spatial$statistics[2,2],sum_49p4_spatial$statistics[2,2]))
  
  xl_64p1<-0
  xu_64p1<-sqrt(sam*0.64)
  yl_64p1<-0
  yu_64p1<-sqrt(sam*0.64)
  
  #extracting x_coordinates
  ind_vec_x_64p1<-c()
  #extracting indices of the x coordinates which are between xl_64p1 and xu_64p1
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p1[1]>xl_64p1){
      if(test_vec_64p1[1]<xu_64p1){
        ind_vec_x_64p1<-c(ind_vec_x_64p1,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_64p1 and yu_64p1
  ind_vec_y_64p1<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p1<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p1[2]>yl_64p1){
      if(test_vec_64p1[2]<yu_64p1){
        ind_vec_y_64p1<-c(ind_vec_y_64p1,i)
      }
    }
  }
  
  x_final_ind_64p1<-ind_vec_x_64p1[which(ind_vec_x_64p1 %in% ind_vec_y_64p1)]
  XY_64p1_dat<-sir_mod0$XYcoordinates[x_final_ind_64p1,]
  x_vec_64p1<-XY_64p1_dat[,1]
  y_vec_64p1<-XY_64p1_dat[,2]
  inftime_64p1<-sir_mod0[["inftime"]][x_final_ind_64p1]
  remtime_64p1<-sir_mod0[["remtime"]][x_final_ind_64p1]
  lambda_64p1<-rep(3,length(x_final_ind_64p1))
  sir_64p1<-as.epidata(type = "SIR", n = length(x_vec_64p1),x = x_vec_64p1, y = y_vec_64p1, 
                       inftime =inftime_64p1,infperiod = lambda_64p1)
  mcmc_outp1<-epimcmc(object = sir_64p1, tmin = min(inftime_64p1), tmax = max(inftime_64p1),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_64p1_spatial<-summary(mcmc_outp1, start = 10001)
  
  
  ## subset 2
  
  xl_64p2<-sqrt(sam)-sqrt(sam*0.64)
  xu_64p2<-sqrt(sam)
  yl_64p2<-0
  yu_64p2<-sqrt(sam*0.64)
  
  #extracting x_coordinates
  ind_vec_x_64p2<-c()
  #extracting indices of the x coordinates which are between xl_64p2 and xu_64p2
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p2[1]>xl_64p2){
      if(test_vec_64p2[1]<xu_64p2){
        ind_vec_x_64p2<-c(ind_vec_x_64p2,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between yl_64p1 and yu_64p2
  ind_vec_y_64p2<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p2<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p2[2]>yl_64p2){
      if(test_vec_64p2[2]<yu_64p2){
        ind_vec_y_64p2<-c(ind_vec_y_64p2,i)
      }
    }
  }
  
  x_final_ind_64p2<-ind_vec_x_64p2[which(ind_vec_x_64p2 %in% ind_vec_y_64p2)]
  XY_64p2_dat<-sir_mod0$XYcoordinates[x_final_ind_64p2,]
  x_vec_64p2<-XY_64p2_dat[,1]
  y_vec_64p2<-XY_64p2_dat[,2]
  inftime_64p2<-sir_mod0[["inftime"]][x_final_ind_64p2]
  remtime_64p2<-sir_mod0[["remtime"]][x_final_ind_64p2]
  lambda_64p2<-rep(3,length(x_final_ind_64p2))
  sir_64p2<-as.epidata(type = "SIR", n = length(x_vec_64p2),x = x_vec_64p2, y = y_vec_64p2, 
                       inftime =inftime_64p2,infperiod = lambda_64p2)
  mcmc_outp2<-epimcmc(object = sir_64p2, tmin = min(inftime_64p2), tmax = max(inftime_64p2),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_64p2_spatial<-summary(mcmc_outp2, start = 10001)
  
  xl_64p3<-0
  xu_64p3<-sqrt(sam*0.64)
  
  yl_64p3<-sqrt(sam)-sqrt(sam*0.64)
  yu_64p3<-sqrt(sam)
  
  #extracting x_coordinates
  ind_vec_x_64p3<-c()
  #extracting indices of the x coordinates which are between xl_64p3 and xu_64p3
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p3[1]>xl_64p3){
      if(test_vec_64p3[1]<xu_64p3){
        ind_vec_x_64p3<-c(ind_vec_x_64p3,i)
      }
    }
  }
  #extracting indices of the x coordinates which are between xl_64p3 and xu_64p3
  ind_vec_y_64p3<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p3<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p3[2]>yl_64p3){
      if(test_vec_64p3[2]<yu_64p3){
        ind_vec_y_64p3<-c(ind_vec_y_64p3,i)
      }
    }
  }
  
  x_final_ind_64p3<-ind_vec_x_64p3[which(ind_vec_x_64p3 %in% ind_vec_y_64p3)]
  XY_64p3_dat<-sir_mod0$XYcoordinates[x_final_ind_64p3,]
  x_vec_64p3<-XY_64p3_dat[,1]
  y_vec_64p3<-XY_64p3_dat[,2]
  inftime_64p3<-sir_mod0[["inftime"]][x_final_ind_64p3]
  remtime_64p3<-sir_mod0[["remtime"]][x_final_ind_64p3]
  lambda_64p3<-rep(3,length(x_final_ind_64p3))
  sir_64p3<-as.epidata(type = "SIR", n = length(x_vec_64p3),x = x_vec_64p3, y = y_vec_64p3, 
                       inftime =inftime_64p3,infperiod = lambda_64p3)
  mcmc_outp3<-epimcmc(object = sir_64p3, tmin = min(inftime_64p3), tmax = max(inftime_64p3),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_64p3_spatial<-summary(mcmc_outp3, start = 10001)
  
  
  xl_64p4<-sqrt(sam)-sqrt(sam*0.64)
  xu_64p4<-sqrt(sam)
  
  yl_64p4<-sqrt(sam)-sqrt(sam*0.64)
  yu_64p4<-sqrt(sam)
  
  #extracting indices of the x coordinates which are between xl_64p4 and xu_64p4
  ind_vec_x_64p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p4[1]>xl_64p4){
      if(test_vec_64p4[1]<xu_64p4){
        ind_vec_x_64p4<-c(ind_vec_x_64p4,i)
      }
    }
  }
  
  #extracting indices of the x coordinates which are between yl_64p4 and yu_64p4
  ind_vec_y_64p4<-c()
  for(i in 1:nrow(sir_mod0$XYcoordinates)){
    test_vec_64p4<-c(sir_mod0$XYcoordinates[i,])
    if(test_vec_64p4[2]>yl_64p4){
      if(test_vec_64p4[2]<yu_64p4){
        ind_vec_y_64p4<-c(ind_vec_y_64p4,i)
      }
    }
  }
  
  x_final_ind_64p4<-ind_vec_x_64p4[which(ind_vec_x_64p4 %in% ind_vec_y_64p4)]
  XY_64p4_dat<-sir_mod0$XYcoordinates[x_final_ind_64p4,]
  x_vec_64p4<-XY_64p4_dat[,1]
  y_vec_64p4<-XY_64p4_dat[,2]
  inftime_64p4<-sir_mod0[["inftime"]][x_final_ind_64p4]
  remtime_64p4<-sir_mod0[["remtime"]][x_final_ind_64p4]
  lambda_64p4<-rep(3,length(x_final_ind_64p4))
  sir_64p4<-as.epidata(type = "SIR", n = length(x_vec_64p4),x = x_vec_64p4, y = y_vec_64p4, 
                       inftime =inftime_64p4,infperiod = lambda_64p4)
  mcmc_outp4<-epimcmc(object = sir_64p4, tmin = min(inftime_64p4), tmax = max(inftime_64p4),
                      niter = 50000, sus.par.ini = 1, beta.ini = 1,
                      Sformula = NULL, pro.sus.var = 0.3, pro.beta.var = 0.1,
                      prior.sus.dist = "gamma", prior.beta.dist = "gamma",
                      prior.sus.par = alphapar2, prior.beta.par = betapar2,
                      adapt = FALSE, acc.rate = NULL)
  sum_64p4_spatial<-summary(mcmc_outp4, start = 10001)
  
  ##extracting the parameter estimate by each subset percentage
  alpha_64p<-mean(c(sum_64p1_spatial$statistics[1,1],sum_64p2_spatial$statistics[1,1],
                    sum_64p3_spatial$statistics[1,1],sum_64p4_spatial$statistics[1,1]))
  
  alpha_sd_64p<-mean(c(sum_64p1_spatial$statistics[1,2],sum_64p2_spatial$statistics[1,2],
                       sum_64p3_spatial$statistics[1,2],sum_64p4_spatial$statistics[1,2]))
  
  beta_64p<-mean(c(sum_64p1_spatial$statistics[2,1],sum_64p2_spatial$statistics[2,1],
                   sum_64p3_spatial$statistics[2,1],sum_64p4_spatial$statistics[2,1]))
  
  beta_sd_64p<-mean(c(sum_64p1_spatial$statistics[2,2],sum_64p2_spatial$statistics[2,2],
                      sum_64p3_spatial$statistics[2,2],sum_64p4_spatial$statistics[2,2]))
  
  ##extract the posterior means and sds of entire epidemic data
  alpha_100p<-sum_100p_spatial$statistics[1,1]
  alpha_sd_100p<-sum_100p_spatial$statistics[1,2]
  beta_100p<-sum_100p_spatial$statistics[2,1]
  beta_sd_100p<-sum_100p_spatial$statistics[2,2]
  
  #finding the average of those mean percentages
  
  alpha_est<-mean(c(alpha_25p,alpha_36p,alpha_49p,alpha_64p))
  
  
  
  alpha_sd_est<-mean(c(alpha_sd_25p,alpha_sd_36p,alpha_sd_49p,alpha_sd_64p))
  
  
  beta_est<-mean(c(beta_25p,beta_36p,beta_49p,beta_64p))
  
  
  
  beta_sd_est<-mean(c(beta_sd_25p,beta_sd_36p,beta_sd_49p,beta_sd_64p))
  
  
  
  
  result<-c(rep(NA,4))
  result[1]<-abs(alpha_est-alpha_100p)
  result[2]<-abs(alpha_sd_est-alpha_sd_100p)
  result[3]<-abs(beta_est-beta_100p)
  result[4]<-abs(beta_sd_est-beta_sd_100p)
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
  res<-avg_corner(sim)
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
write.csv(df_abs_biases,"/home/thethtetchan.nyein/avg_corner_100_abs_biases.csv", row.names = FALSE)








