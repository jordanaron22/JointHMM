################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(dplyr)
library(numDeriv)
library(Matrix)
library(Hmisc)

library(survex)
library(tidyverse)

beta_bool <- 0
real_data <- 1
sim_size <- 0
randomize_init <- 1
set_seed <- 1
tobit <- 1
leave_out <- 0
# incl_surv <- 2
incl_light <- 1
incl_act <- 1
load_data <- 0
check_tran <- 0
wake_sleep <- 0
weekend_only <- 0
misspecification <- 0
save_space <- 0
single_day <- 0


epsilon <- 1e-3
lepsilon <- log(epsilon)

vcovar_num <- 2

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(paste("Seed",sim_num))

if (is.na(sim_num)){sim_num <- 11}
if (set_seed){set.seed(sim_num)}

print("Command line arguments:")
print(commandArgs(T))

# command_args <- as.character(commandArgs(TRUE)[1])
# command_args <- as.numeric(unlist(strsplit(command_args, "")))

# mix_num <- command_args[1]
mix_num <- as.numeric(commandArgs(TRUE)[1])
print(paste("mixnum",mix_num))

# incl_surv <- command_args[2]
incl_surv <- as.numeric(commandArgs(TRUE)[2])
print(paste("Include Surv",incl_surv))

# real_data <- command_args[3]
# print(paste("Real Data",real_data))

# sim_size <- command_args[3]
# print(paste("Single day",sim_size))
# 
# misspecification <- command_args[2]
# print(paste("misspecification",misspecification))
# # 
# incl_light <- command_args[3]
# print(paste("Include Light",incl_light))
#
# wake_sleep <- command_args[1]
# print(paste("Cycle Only?",wake_sleep))
#

if(is.na(mix_num)){mix_num <- 2}
if(is.na(incl_surv)){incl_surv <- 2}
# if(is.na(real_data)){real_data <- 1}
# if(is.na(sim_size)){sim_size <- 0}
# if(is.na(misspecification)){misspecification <- 0}
# if(is.na(load_data)){load_data <- 0}

print(paste("Sim Seed:",sim_num,"HMM Num:",mix_num))

model_name <- "JMHMM"
if (!real_data){model_name <- paste0(model_name,"SimSize",sim_size)}
if (!tobit){model_name <- paste0(model_name,"NonTob")}
if (leave_out){model_name <- paste0(model_name,"LeaveOut")}
if (incl_surv==0){model_name <- paste0(model_name,"NoSurv")}
if (incl_surv==1){model_name <- paste0(model_name,"HalfSurv")}
if (!incl_light){model_name <- paste0(model_name,"NoLight")}
if (!incl_act){model_name <- paste0(model_name,"NoAct")}
if (wake_sleep){model_name <- paste0(model_name,"OnlyCycle")}
if (weekend_only){model_name <- paste0(model_name,"WeekendOnly")}
if (misspecification){model_name <- paste0(model_name,"Misspecified",misspecification)}
if (single_day){model_name <- paste0(model_name,"SingleDay",single_day)}

model_name <- paste0(model_name,"Mix",mix_num,"Seed",sim_num,".rda")
print(model_name)
################## Functions ##################
readCpp <- function(path) {
  tryCatch(
    {
      sourceCpp(file = path)
    },
    error = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      message("Done")
    }
  )
}

expit <- function(x){
  to_ret <- exp(x) / (1+exp(x))
  if (is.na(to_ret)){return(1)}
  return(to_ret)
}

logit <- function(x){
  return(log(x/(1-x)))
}

Params2TranVectorT <- function(re_ind,len,params_tran){
  return(t(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,re_ind=re_ind)))
}


Params2TranVectorTresid <- function(re_ind,len,params_tran){
  return(t(sapply(c(1:(len)),FUN = Params2Tran,params_tran = params_tran,re_ind=re_ind)))
}

TranByTimeVec <- function(re_ind, params_tran,time_vec){
  return(lapply(time_vec, Params2Tran, params_tran = params_tran,re_ind=re_ind))
}

Param2TranHelper <- function(p12,p21){
  tran <- matrix(0,2,2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1- tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  return(tran)
}

Params2Tran <- function(params_tran,time,re_ind){
  
  param_matrix <- matrix(params_tran[re_ind,],ncol=3,nrow=2, byrow = T)
  tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,2]*cos(2*pi*time/96)+param_matrix[1,3]*sin(2*pi*time/96),
                           param_matrix[2,1]+param_matrix[2,2]*cos(2*pi*time/96)+param_matrix[2,3]*sin(2*pi*time/96))


  return(tran)
}

#Calculates mu1|2
CalcCondMean <- function(mu1,sig1,mu2,sig2,bivar_corr,obs2){
  return(mu1 + bivar_corr*(sig1/sig2)*(obs2-mu2))
}

#Calculates sig1|2
CalcCondSig <- function(sig1,bivar_corr){
  return(sig1*sqrt(1-bivar_corr^2))
}

Case4 <- function(act_obs,mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD){
  
  mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs)
  sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
  
  lognorm_dens <- dnorm(act_obs,mu_act,sig_act) * 
    pnorm(light_LOD,mu_light_cond,sig_light_cond)
  return(lognorm_dens)
}

CalcLintegralMat <- function(emit_act,emit_light,corr_mat,lod_act,lod_light){
  mix_num <- dim(emit_act)[3]
  if (is.na(mix_num)){mix_num <- 1}
  
  lintegral_mat <- array(NA,dim = c(mix_num,2,2))
  #j is week/weekend
  for (j in 1:2){
    for (i in 1:mix_num){
      
      # lb <- mu_act - 5*sig_act
      # lb <- min(lb,-10)
      
      lintegral_mat[i,1,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[1,1,i,j],emit_act[1,2,i,j],
                                          emit_light[1,1,i,j],emit_light[1,2,i,j],
                                          corr_mat[i,1,j],lod_light)[[1]])
      
      lintegral_mat[i,2,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[2,1,i,j],emit_act[2,2,i,j],
                                          emit_light[2,1,i,j],emit_light[2,2,i,j],
                                          corr_mat[i,2,j],lod_light)[[1]])
    }
  }
  
  lintegral_mat[lintegral_mat == -Inf] <- -9999
  
  return(lintegral_mat)
}

logClassification <- function(time,act,light,mu_act,sig_act,mu_light,sig_light,lod_act,lod_light,bivar_corr,lintegral){

  #CASE 1
  if (act[time] > lod_act & light[time] > lod_light){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      dnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 2 
  if (act[time] > lod_act & light[time] == lod_light){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      pnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 3
  if (act[time] == lod_act & light[time] > lod_light){
    
    mu_act_cond <- CalcCondMean(mu_act,sig_act,mu_light,sig_light,bivar_corr,light[time])
    sig_act_cond <- CalcCondSig(sig_act,bivar_corr)
    
    lognorm_dens <- pnorm(act[time],mu_act_cond,sig_act_cond,log = T) + 
      dnorm(light[time],mu_light,sig_light,log = T)
  }
  
  
  #CASE 4 
  if (act[time] == lod_act & light[time] == lod_light){
    
    lognorm_dens <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                  mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_light)[[1]])
  }

  return(lognorm_dens)
}

ChooseMixture<- function(mixture_vec){
  return(which(mixture_vec == 1))
}

SimulateMC <- function(day_length,init,tran_list_ind,mixture_ind,vcovar_vec){
  hidden_states <- numeric(day_length)
  
  for (i in 1:day_length){
    
    tran <- tran_list_ind[[vcovar_vec[i]]][[(i-1)%%96+1]]
    
    if (i == 1) {
      hidden_states[1] <- rbinom(1,1,init[mixture_ind,2])
    } else {
      hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
    }
  }
  
  
  return(hidden_states)
}

finv <- function(lam,time, randu, xb_ind){
  return((1-pexp(time,rate = lam))^(exp(xb_ind)) - randu)
}

SimSurvival <- function(mixture_vec,beta_vec,beta_age_true,beta_covar_sim,age_vec,surv_covar_sim,lam = 1/20){
  
  num_of_people <- length(mixture_vec)
  
  failure_times <- numeric(num_of_people)
  censor_times <- numeric(num_of_people)
  
  for(i in 1:num_of_people){
    xb <- beta_vec[mixture_vec[i]] + beta_age_true*age_vec[i] + beta_covar_sim[surv_covar_sim[i]]
    
    evnt<-runif(1)
    cens<-runif(1)
    failure_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=evnt, xb_ind=xb, extendInt = "yes")$root
    censor_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=cens, xb_ind=xb, extendInt = "yes")$root
    
  }
  
  time <- pmin(failure_times, censor_times)
  event <- as.integer(failure_times<censor_times)
  
  return(list(time,event))
}

# rev_order <- c(dim(tran_expand_df)[2]:1
# for (i in 1:dim(tran_expand_df)[1]){
#   list_index <- paste0("[[",tran_expand_df[2,rev_order],"]]",collapse="")
# }


GenTranColVecList <- function(params_tran_array,mix_num,vcovar_num){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  
  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (mixture_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      params_tran <- params_tran_array[,,vcovar_ind]
      if (mix_num == 1) {params_tran <- matrix(params_tran,nrow = 1)}
      vcovar_tran_list[[vcovar_ind]] <- Params2TranVectorT(mixture_ind,len,params_tran)
      
    }
    mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
  }
  
  return(mixture_vcovar_tran_list)
}

GenTranList <- function(params_tran_array,time_vec,mix_num,vcovar_num){
  
  
  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (mixture_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      
      params_tran <- params_tran_array[,,vcovar_ind]
      if (mix_num == 1) {params_tran <- matrix(params_tran,nrow = 1)}
      vcovar_tran_list[[vcovar_ind]] <- TranByTimeVec(re_ind = mixture_ind,
                                                      params_tran = params_tran,
                                                      time_vec = time_vec)
      
    }
    mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
  }
  
  return(mixture_vcovar_tran_list)
}


SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,
                        emit_act,emit_light,corr_mat,
                        lod_act,lod_light, nu_mat,beta_vec_true,beta_age_true,beta_covar_sim,
                        missing_perc,lambda_act_mat,lambda_light_mat,
                        misspecification){
  
  mix_num <- dim(emit_act)[3]
  
  age_vec <- floor(runif(num_of_people,10,81))
  surv_covar_sim_wide <- t(rmultinom(num_of_people,1,c(1/3,1/3,1/3)))
  surv_covar_sim <- apply(surv_covar_sim_wide, 1, which.max)
  
  statact_vec <- floor(runif(num_of_people,0,16))
  modact_vec <- rbinom(num_of_people,1,.5)
  # nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2,modact_vec)
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)
  
  pi_l_true <- CalcPi(nu_mat,nu_covar_mat)
  
  if (misspecification == 1){
    pi_l_true <- cbind(pi_l_true[,1],pi_l_true[,1])
    pi_l_true[,1] <- 1
    pi_l_true[,2] <- 0
  } else if (misspecification > 1){
    pi_l_true <- pi_l_true[,1:misspecification]
    pi_l_true <- pi_l_true/ rowSums(pi_l_true)
  }
  if (dim(pi_l_true)[2] > 1){
    mixture_vec <- rMultinom(pi_l_true,1)
  } else {
    mixture_vec <- matrix(1,nrow = dim(pi_l_true)[1],ncol = 1)
  }
  
  
  if (day_length == 192){
    vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48))
  } else if (day_length == 384){
    vcovar_vec <- c(rep(0,96*2),rep(1,96*2))
  } else if (day_length == 672){
    vcovar_vec <- c(rep(0,96*5),rep(1,96*2))
  } else {
    vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48))
    vcovar_vec <- rep(vcovar_vec,day_length/192)
  }
  
  vcovar_mat <- replicate(num_of_people, vcovar_vec)

  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
  
  
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)
    
    mixture_ind <- mixture_vec[ind]
    
    tran_vcovar_list <- tran_list[[mixture_ind]]
    vcovar_vec_ind <- vcovar_mat[,ind]+1
    
    hidden_states <- SimulateMC(day_length,init,tran_vcovar_list,mixture_ind,vcovar_vec_ind)
    
    for (i in 1:day_length){
      
      mu_act <- emit_act[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]] 
      
      sig_act <- emit_act[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]] 
      
      mu_light <- emit_light[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]] 
      
      sig_light <- emit_light[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]] 
      
      bivar_corr <- corr_mat[mixture_ind,hidden_states[i] + 1,vcovar_vec_ind[i]] 
      
      if (tobit){
        sigma_mat <- matrix(c(sig_act^2,bivar_corr * sig_act* sig_light,
                              bivar_corr * sig_act* sig_light,sig_light^2),2,2,byrow = T)
        
        act_light <- mvrnorm(n = 1,
                             mu = c(mu_act,mu_light),
                             Sigma = sigma_mat)
        
        activity[i] <-act_light[1]
        light[i] <-act_light[2]
      } else {
        
        lambda_act <- lambda_act_mat[mixture_ind,hidden_states[i]+1,vcovar_vec_ind[i]]
        lambda_light <- lambda_light_mat[mixture_ind,hidden_states[i]+1,vcovar_vec_ind[i]]
        
        activity[i] <- rnorm(1,mu_act,sig_act)
        light[i] <- rnorm(1,mu_light,sig_light)
        
        if(rbinom(1,1,lambda_act)){activity[i] <- lod_act}
        if(rbinom(1,1,lambda_light)){light[i] <- lod_light}
        
        
      }
      
        
    }
    
    if (ind == 1){
      hidden_states_matrix <- hidden_states
      activity_matrix <- activity
      light_matrix <- light
    } else {
      hidden_states_matrix <- cbind(hidden_states_matrix,hidden_states)
      activity_matrix <- cbind(activity_matrix,activity)
      light_matrix <- cbind(light_matrix,light)
    }
  }
  
  # init_emp <- init
  # params_tran_array_emp <- params_tran_array
  # params_tran_array_fcovar_emp <- params_tran_array_fcovar
  # emit_act_emp <- emit_act
  # emit_light_emp <- emit_light
  # emit_act_fcovar_emp <- emit_act_fcovar
  # emit_light_fcovar_emp <- emit_light_fcovar
  # corr_mat_emp <- corr_mat
  # pi_l_emp <- pi_l
  # 
  # for (re_ind in 1:mix_num){
  #   mixture_indicator <- (mixture_vec == re_ind)
  #   act_mixture <- activity_matrix[,mixture_indicator]
  #   light_mixture <- light_matrix[,mixture_indicator]
  #   mc_wake_mixture <- hidden_states_matrix[,mixture_indicator]==0
  #   mc_sleep_mixture <- hidden_states_matrix[,mixture_indicator]==1
  #   
  #   init_emp[re_ind,1] <- sum(mc_wake_mixture[1,])/length(mc_wake_mixture[1,])
  #   init_emp[re_ind,2] <- 1 - init_emp[re_ind,1] 
  #   
  #   emit_act_emp[1,1,re_ind] <- mean(act_mixture[mc_wake_mixture])
  #   emit_act_emp[1,2,re_ind]  <- sqrt(var(act_mixture[mc_wake_mixture]))
  #   emit_act_emp[2,1,re_ind]  <- mean(act_mixture[mc_sleep_mixture])
  #   emit_act_emp[2,2,re_ind]  <- sqrt(var(act_mixture[mc_sleep_mixture]))
  #   
  #   emit_light_emp[1,1,re_ind] <- mean(light_mixture[mc_wake_mixture])
  #   emit_light_emp[1,2,re_ind] <- sqrt(var(light_mixture[mc_wake_mixture]))
  #   emit_light_emp[2,1,re_ind] <- mean(light_mixture[mc_sleep_mixture])
  #   emit_light_emp[2,2,re_ind] <- sqrt(var(light_mixture[mc_sleep_mixture]))
  #   
  #   corr_mat_emp[re_ind,1] <- cor(act_mixture[mc_wake_mixture],light_mixture[mc_wake_mixture])
  #   corr_mat_emp[re_ind,2] <- cor(act_mixture[mc_sleep_mixture],light_mixture[mc_sleep_mixture])
  #   
  #   pi_l_emp[re_ind] <- sum(mixture_indicator)/num_of_people
  # }
  
  surv_list <- SimSurvival(mixture_vec,beta_vec_true,beta_age_true,beta_covar_sim,age_vec,surv_covar_sim)
  # if (misspecification == 1 | all(mixture_vec== 1) ){
  #   cox_mod <- coxph(Surv(surv_list[[1]], surv_list[[2]]) ~age_vec)
  # } else {
  #   cox_mod <- coxph(Surv(surv_list[[1]], surv_list[[2]]) ~age_vec + as.factor(mixture_vec))
  # }
  # 
  # beta_vec_long_emp <- as.vector(as.vector(cox_mod$coefficients))
  
  
  # 
  # beta_vec_emp <- beta_vec_long_emp[1:mix_num]
  # beta_vec_emp[1] <- 0
  # 
  # beta_vec_fcovar_emp <- beta_vec_long_emp[(mix_num+1):length(beta_vec_long_emp)]
  # beta_vec_fcovar_emp[1] <- 0
  # 
  # if (fcovar_num == 1){
  #   beta_mat_emp <- replicate(fcovar_num, beta_vec_emp)
  # } else {
  #   beta_mat_emp <- replicate(fcovar_num, beta_vec_emp) + t(replicate(mix_num,beta_vec_fcovar_emp))
  # }
  
  
  
  # beta_vec_emp <- c(0,summary(cph_fit)$coefficients[,1])
  
  # emp_params <- list(init_emp,params_tran_array_emp,params_tran_array_fcovar_emp,emit_act_emp,
  #                    emit_light_emp,corr_mat_emp,pi_l_emp,beta_mat_emp)
  

  if (tobit){
    light_matrix[light_matrix<lod_light] <- lod_light
    activity_matrix[activity_matrix<lod_act] <- lod_act
  }
  
  
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)
  
  activity_matrix[act_missing==1] <- NA
  light_missing[light_missing==1] <- NA
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,mixture_vec,
              age_vec,nu_covar_mat,vcovar_mat,surv_list,surv_covar_sim))
}

CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[ind,re_ind])
    }
  }
  
  ind_like_mat <- apply(alpha_beta,c(1,4),logSumExp)
  
  weight_array_wake <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  weight_array_sleep <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  for (ind in 1:dim(alpha_beta)[4]){
    for (t in 1:dim(alpha_beta)[1]){
      weight_array_wake[t,ind,] <- alpha_beta[t,1,,ind] - ind_like_mat[t,ind]
      weight_array_sleep[t,ind,] <- alpha_beta[t,2,,ind] - ind_like_mat[t,ind]
    }
  }
  
  return(list(weight_array_wake,weight_array_sleep))
}

CalcInit <- function(alpha, beta,pi_l,log_sweights_vec){
  
  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_1_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_mat <- matrix(0,dim(pi_l)[2],2)
  
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)
  
  
  for(ind in 1:length(alpha)){ 
    ind_like <- ind_like_vec[ind]
    
    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l[ind,]) - ind_like + log_sweights_vec[ind]
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l[ind,]) - ind_like + log_sweights_vec[ind]

  }
  
  for (re_ind in 1:(dim(pi_l)[2])){
    init_0 <- logSumExp(init_0_vec[,re_ind])
    init_1 <- logSumExp(init_1_vec[,re_ind])
    init_vec <- exp(c(init_0,init_1) - logSumExp(c(init_0,init_1)))
    init_mat[re_ind,] <- init_vec
  }
  
  return(init_mat)
  
}

CalcProbRE <- function(alpha,pi_l){
  
  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)
  
  
  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[ind,re_ind])
    }
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
    
  }
  
  return(re_weights)
  
}

CalcTranHelper <- function(act, light, tran_list_mat, emit_act, emit_light, 
                           ind_like_vec, alpha, beta, lod_act, lod_light, 
                           corr_mat, lintegral_mat, pi_l,vcovar_mat,
                           lambda_act_mat, lambda_light_mat, tobit){
  
  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  
  tran_vals_re_array <- array(NA,c(2,2,len - 1,num_people,num_re))
  
  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))
  
  for(init_state in 1:2){
    for(new_state in 1:2){
      for (clust_i in 1:num_re){
        tran_vals_re_array[init_state,new_state,,,clust_i] <- CalcTranHelperC(init_state = init_state-1, new_state = new_state-1,act = act,
                                                          light = light,tran_list_mat = tran_list_mat,
                                                          emit_act_week = emit_act_week,emit_light_week = emit_light_week,
                                                          emit_act_weekend = emit_act_weekend,emit_light_weekend = emit_light_weekend,
                                                          ind_like_vec = ind_like_vec,
                                                          alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                                          corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                                          clust_i = clust_i-1, vcovar_mat = vcovar_mat,
                                                          lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
      }
    }
  }
  
  return(tran_vals_re_array)
}

Symmetricize <- function(mat){
  mat <- mat + t(mat)
  diag(mat) <- diag(mat)/2
  return(mat)
}

CalcGradHess <- function(gradient,hessian_vec,cos_part_vec,sin_part_vec,cos_sin_part){
  grad_array <- array(0,dim = c(6,mix_num,vcovar_num))
  hess_array <- array(0,dim = c(6,6,mix_num,vcovar_num))
  
  for(re_ind in 1:mix_num){
    for(vcovar_ind in 1:vcovar_num){
      
      grad_array[,re_ind,vcovar_ind] <- as.vector(t(gradient[,,re_ind,vcovar_ind]))
      
      hessian_vec_re <- hessian_vec[,,re_ind,vcovar_ind]
      hess_upper <- matrix(0,3,3)
      hess_lower <- matrix(0,3,3)
      
      #### HESS 1
      
      diag(hess_upper) <- c(hessian_vec_re[1,])
      hess_upper[1,2] <- cos_part_vec[1,re_ind,vcovar_ind]
      hess_upper[1,3] <- sin_part_vec[1,re_ind,vcovar_ind]
      hess_upper[2,3] <- cos_sin_part[1,re_ind,vcovar_ind]
      hess_upper <- Symmetricize(hess_upper)
      
      #### HESS 2
      
      diag(hess_lower) <- c(hessian_vec_re[2,])
      hess_lower[1,2] <- cos_part_vec[2,re_ind,vcovar_ind]
      hess_lower[1,3] <- sin_part_vec[2,re_ind,vcovar_ind]
      hess_lower[2,3] <- cos_sin_part[2,re_ind,vcovar_ind]
      hess_lower <- Symmetricize(hess_lower)
      
      #######
      
      hess_array[1:3,1:3,re_ind,vcovar_ind] <- hess_upper
      hess_array[4:6,4:6,re_ind,vcovar_ind] <- hess_lower
      # params_tran_working[re_ind,] <- params_tran_working[re_ind,] - solve(-hessian_re,-gradient_re)
    }
    
  }
  
  return(list(grad_array,hess_array))
}

  
CalcTranCHelper <- function(alpha,beta,act,light,params_tran_array,emit_act,emit_light,
                      corr_mat,pi_l,lod_act,lod_light,lintegral_mat, vcovar_mat,
                      lambda_act_mat, lambda_light_mat, tobit, check_tran,likelihood){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  params_tran_array_working <- params_tran_array
  
  gradient <- array(0,c(2,3,mix_num,vcovar_num))
  hessian_vec <- array(0,c(2,3,mix_num,vcovar_num))
  cos_part_vec <- array(0,c(2,mix_num,vcovar_num))
  sin_part_vec <- array(0,c(2,mix_num,vcovar_num))
  cos_sin_part <- array(0,c(2,mix_num,vcovar_num))
  
  # tran_list_mat <- lapply(c(1:mix_num),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_list_mat <- GenTranColVecList(params_tran_array,mix_num,vcovar_num)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))

  tran_vals_re_array <- CalcTranHelper(act = act,
                                       light = light,tran_list_mat = tran_list_mat,
                                       emit_act = emit_act,emit_light = emit_light,
                                       ind_like_vec = ind_like_vec,
                                       alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                       corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                       vcovar_mat = vcovar_mat[-1,],lambda_act_mat, lambda_light_mat, tobit)
  
  
  cos_vec <- cos(2*pi*c(2:(len))/96)
  sin_vec <- sin(2*pi*c(2:(len))/96)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals <- tran_vals_re_array[init_state,new_state,,,]
      if (mix_num  == 1){dim(tran_vals)[3] <- 1}
      
      for (ind in 1:length(alpha)){
        
        vcovar_vec <- vcovar_mat[-1,ind]
        vcovar_vecR <- vcovar_vec + 1
        
        for(re_ind in 1:mix_num){
          
          if (vcovar_num == 1){
            tran_mat <- tran_list_mat[[re_ind]][[1]]
          } else if (vcovar_num == 2){
            tran_mat_week <- tran_list_mat[[re_ind]][[1]]
            tran_mat_weekend <- tran_list_mat[[re_ind]][[2]]
            tran_mat <- tran_mat_week * (1-vcovar_vec) + tran_mat_weekend * vcovar_vec
          }
            
          
          if(init_state == 1 & new_state == 1){
            # tran_prime <- -tran[1,2]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- -tran_mat[,3]
            tran_prime_prime <- -tran_mat[,3]*tran_mat[,1]
            
          } else if(init_state == 1 & new_state == 2){ 
            # tran_prime <- tran[1,1]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- tran_mat[,1]
            tran_prime_prime <- -tran_mat[,3]*tran_mat[,1]
            
          } else if(init_state == 2 & new_state == 2){ 
            # tran_prime <- -tran[2,1]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- -tran_mat[,2]
            tran_prime_prime <- -tran_mat[,2] * tran_mat[,4]
            
          } else if(init_state == 2 & new_state == 1){ 
            # tran_prime <- tran[2,2]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- tran_mat[,4]
            tran_prime_prime <- -tran_mat[,2] * tran_mat[,4]
          }
          
          
          for (vcovar_ind in 1:vcovar_num){

            # Double derivate of b's
            gradient[init_state,1,re_ind,vcovar_ind] <- gradient[init_state,1,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*(vcovar_vecR==vcovar_ind))
            gradient[init_state,2,re_ind,vcovar_ind] <- gradient[init_state,2,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            gradient[init_state,3,re_ind,vcovar_ind] <- gradient[init_state,3,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            hessian_vec[init_state,1,re_ind,vcovar_ind] <- hessian_vec[init_state,1,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,2,re_ind,vcovar_ind] <- hessian_vec[init_state,2,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,3,re_ind,vcovar_ind] <- hessian_vec[init_state,3,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2*(vcovar_vecR==vcovar_ind))
            
            cos_part_vec[init_state,re_ind,vcovar_ind] <- cos_part_vec[init_state,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            sin_part_vec[init_state,re_ind,vcovar_ind] <- sin_part_vec[init_state,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            cos_sin_part[init_state,re_ind,vcovar_ind] <- cos_sin_part[init_state,re_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec*(vcovar_vecR==vcovar_ind))
            
          }
        }
      }
    }
  }
  
  
  grad_hess_list <- CalcGradHess(gradient,hessian_vec,cos_part_vec,sin_part_vec,cos_sin_part)
  grad_array <- grad_hess_list[[1]]
  hess_array <- grad_hess_list[[2]]
  
  return(list(grad_array,hess_array))
  
  # return(LM(grad_array,hess_array,params_tran_array,check_tran,likelihood,pi_l))

}

LM <- function(grad_array,hess_array,params_tran_array,check_tran,likelihood,pi_l, step_size = 1){
  params_tran_array_new <- params_tran_array
  new_likelihood <- -Inf
  # while (new_likelihood - likelihood < -1e-5 ){
  
  for (re_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind],grad_array[,re_ind,vcovar_ind])
      
      step_size <- 1
      while(max(abs(inf_fact)) > 3 | all(inf_fact == -1)){
        step_fact <- matrix(0,6,6)
        diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) * step_size
        if(all(inf_fact == -1)){diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) + step_size}
        
        inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind]+step_fact,grad_array[,re_ind,vcovar_ind])
        step_size <- step_size * 10
      }
      params_tran_array_new[re_ind,,vcovar_ind] <- params_tran_array_new[re_ind,,vcovar_ind] - inf_fact
    }
  }
  
  if (check_tran){
    tran_list <- GenTranList(params_tran_array_new,c(1:day_length),mix_num,vcovar_num)
    alpha <- Forward(act = act,light = light,
                     init = init,tran_list = tran_list,
                     emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
                     emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
                     lod_act = lod_act, lod_light = lod_light, 
                     corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
                     event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                     lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                     surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                     lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit,incl_surv = incl_surv*beta_bool)
    new_like <- CalcLikelihood(alpha,pi_l)
    
    if (new_like < likelihood){
      return(list(params_tran_array,alpha))
    }
    
  }
  
  return(params_tran_array_new)
}
# mc_state <- 1
# re_ind <- 1
# vcovar_ind <- 1
# weights_array <- weights_array_wake
# EmitLogLike(act = act, light = light,
#             mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
#             sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
#             mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
#             sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
#             # bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
#             bivar_corr = 1,
#             lod_act = lod_act, lod_light = lod_light,
#             vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
#             weights_mat = as.vector(weights_array[,,re_ind]))

EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,vcovar_mat,vcovar_ind,weights_mat){
  
  #lower should be -Inf, but had some divergence issuer
  lb <- mu_act - 5*sig_act
  lb <- min(lb,-10)
  lintegral <- log(integrate(Case4,lower = lb,upper = lod_act,
                             mu_act,
                             sig_act,
                             mu_light,
                             sig_light,
                             bivar_corr,lod_light)[[1]])
  
  if (lintegral == -Inf){
    lintegral <- -9999
  }
  
  vcovar_vec <- as.vector(vcovar_mat)
  vcovar_vec_indicator <- vcovar_vec == vcovar_ind
  
  act_vec <- as.vector(act)[vcovar_vec_indicator]
  light_vec <- as.vector(light)[vcovar_vec_indicator]
  
  log_like <- logClassificationCTobit(act_vec,light_vec,
                                 mu_act,
                                 sig_act,
                                 mu_light,
                                 sig_light,
                                 lod_act,lod_light,bivar_corr,lintegral)
  
  log_like[log_like == -Inf] <- -9999

  return(-sum(log_like * weights_mat[vcovar_vec_indicator]))
}

CalcActMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(-10,10), act = act, light = light, 
                     # mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind], 
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcActSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(0.1,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     # sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind], 
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}


CalcLightMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(-30,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     # mu_light = emit_light[mc_state,1,re_ind,vcovar_ind], 
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  
  return(mu_act)
}

CalcLightSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(0.01,20), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind], 
                     #sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}



CalcBivarCorr <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(-.999,.999), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind], 
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     # bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}


UpdateNorm <- function(FUN,mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat){
  opt_param_mat <- matrix(0,mix_num,2)
  
  if(mc_state == 1){weights_array <- weights_array_wake}
  if(mc_state == 2){weights_array <- weights_array_sleep}
  
  for (re_ind in 1:dim(emit_act)[3]){
    for(vcovar_ind in 1:2){
      opt_param_mat[re_ind,vcovar_ind] <- FUN(mc_state,vcovar_ind,
                                              act,light,
                                              emit_act,emit_light,
                                              corr_mat,lod_act,lod_light,
                                              weights_array,re_ind,vcovar_mat)
    }
  }
  
  return(opt_param_mat)
  
}

SolveCatch <- function(hess,grad) {
  tryCatch(
    {
      solve(hess,grad,tol=1e-50)
    },
    error = function(cond) {
      # message("Non-Invertible Matrix")
      numeric(length(grad))-1
    },
    warning = function(cond) {
      NULL
    },
    finally = {}
  )
}


SolveLM <- function(hess,grad,step_size) {
  tryCatch(
    {
      diag(hess) <- diag(hess) * step_size

      solve(hess,grad,tol=1e-50)
    },
    error = function(cond) {
      # message("Non-Invertible Matrix Surv")
      numeric(length(grad))
    },
    warning = function(cond) {
      NULL
    },
    finally = {}
  )
}

CalcLikelihood <- function(alpha,pi_l){
  return(sum(CalcLikelihoodIndVec(alpha,pi_l)))
}

CalcLikelihoodIndVec <- function(alpha,pi_l){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  #i is number of people
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,pi_l,i,num_obs)))
    like_vec <- c(like_vec,ind_like)
  }
  return(like_vec)
}

SumOverREIndTime <- function(fb,pi_l,ind,time, add_re = T){
  
  fb_ind <- fb[[ind]]
  
  fb_sum <- numeric(2)
  if (add_re){
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(pi_l[ind,])))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(pi_l[ind,])))
  } else {
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,]))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,]))
  }
  
  return(fb_sum)
}

IndLike <- function(alpha,pi_l,ind,len){
  likelihood <- logSumExp(SumOverREIndTime(alpha,pi_l,ind,len))
  return(likelihood)
}

CalcBIC <- function(new_likelihood,mix_num,act,light){
  
  #init tran emit pi surv
  num_of_param <- mix_num + 
    mix_num*6*2 +
    mix_num*4*2*2+ 
    (mix_num-1)*2 + 
    mix_num-1 + 1
  
  bic <- num_of_param * log(sum(!is.na(act)) + sum(!is.na(light))) - (2 * new_likelihood)
  return(bic)
}

list1 <- list(matrix(1:4, nrow = 2, ncol = 2), matrix(5:8, nrow = 2, ncol = 2))
list2 <- list(matrix(9:12, nrow = 2, ncol = 2), matrix(13:16, nrow = 2, ncol = 2))

# Use mapply to multiply corresponding matrices


SurvCovarRiskScore <- function(surv_coef,surv_covar,set){
  risk_score <- 0
  if (length(set) < 2){
    if (set == 0){set <- c(1:num_of_people)}
  }
  
  for (i in 1:length(surv_coef)){
    
    if (length(surv_coef[[i]]) == 1){
      risk_score <- risk_score + surv_coef[[i]] * surv_covar[[i]][set]
    } else {
      risk_score <- risk_score + surv_covar[[i]][set,]%*%surv_coef[[i]] 
    }
      
  }
  return(risk_score)
}


SurvCovarRiskVec <- function(surv_covar,surv_coef){
  surv_covar[[1]] <- matrix(surv_covar[[1]],ncol = 1)
  surv_coef[[1]] <- matrix(surv_coef[[1]],ncol = 1)
  surv_covar_risk_vec <- rowSums(mapply(function(x, y) x %*% y,  surv_covar,surv_coef, SIMPLIFY = T))
  return(surv_covar_risk_vec)
}

CalcBLHaz <- function(surv_coef,beta_vec, re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar){
  n <- length(surv_event)
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)
  
  for (time_ind in 1:n){
    risk_set <- surv_time >= surv_time[time_ind]
    # linear_surv_covar_risk <- SurvCovarRiskScore(surv_coef,surv_covar,risk_set)
    linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]
    
    if (dim(re_prob)[2] != 1){
      denom <- sum((re_prob[risk_set,] %*% exp(beta_vec)) * exp(linear_surv_covar_risk))
    } else {
      denom <- sum(exp(linear_surv_covar_risk))
    }
      
    if (denom > .Machine$double.xmax){denom <- .Machine$double.xmax}
    bline_vec[time_ind] <- surv_event[time_ind]/denom
  }
  
  for(time_ind in 1:n){
    anti_risk_set <- surv_time <= surv_time[time_ind]
    cbline_vec[time_ind] <- sum(bline_vec[anti_risk_set])
  }
  
  return(list(bline_vec,cbline_vec))
}

SurvLike <- function(beta_vec,surv_covar_risk_vec,surv_coef){
  bhaz_vec <<- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
  bline_vec <<- bhaz_vec[[1]]
  cbline_vec <<- bhaz_vec[[2]]
  
  
  loglike <- sum(log(bline_vec^surv_event) + 
                   ((re_prob %*% beta_vec)+surv_covar_risk_vec) * surv_event - 
                   cbline_vec*exp((re_prob %*% beta_vec)+surv_covar_risk_vec))
  
  return(-loglike)
}


# CalcBetaHelper <- function(beta_vec_long){
#   
#   #NEED TO REDO THIS IF USING
#   # gradient <- grad(SurvLike, beta_vec_long)
#   # hessian <- hessian(SurvLike, beta_vec_long)
#   
#   nr_fact <- SolveCatch(hessian,gradient)
#   step_size <- 1
#   step_mat <- matrix(0,dim(hessian)[1],dim(hessian)[2])
#   # while(max(abs(nr_fact)) > 1 | all(nr_fact == -1)){
#   while(all(nr_fact == -1)){
#     diag(step_mat) <- (diag(hessian) * (step_size-1)) + .01
#     nr_fact <- SolveCatch(hessian + step_mat,gradient)
#     step_size <- step_size * 10
#   }
#   
#   beta_vec_long <- beta_vec_long - nr_fact
#   return(beta_vec_long)
# }

OutofBetaSurvCoef <- function(beta_surv_coef,surv_coef_len){
  beta_vec <- c(0,beta_surv_coef[2:mix_num])
  surv_coef_new <- list()
  surv_coef_new <- append(surv_coef_new,beta_surv_coef[[1]])
  
  
  surv_coef_len_alt <- cumsum(surv_coef_len-1)
  
  
  for (i in 1:(length(surv_coef_len)-1)) {
    coef_vec <- c(0,beta_surv_coef[(mix_num+1+surv_coef_len_alt[i]):(mix_num+surv_coef_len_alt[i+1])])
    surv_coef_new <- append(surv_coef_new,list(coef_vec))
  }
  
  return(list(beta_vec,surv_coef_new))
}


IntoBetaSurvCoef <- function(beta_vec,surv_coef){
  beta_surv_coef_len <- mix_num + length(unlist(surv_coef)) - length(surv_coef)
  beta_surv_coef <- numeric(beta_surv_coef_len)
  
  beta_surv_coef[1] <- surv_coef[[1]]
  beta_surv_coef[2:mix_num] <- beta_vec[-1]
  
  altered_surv_coef <- surv_coef[-1]
  for (i in 1:length(altered_surv_coef)){
    altered_surv_coef[[i]] <- altered_surv_coef[[i]][-1]
  }
  
  beta_surv_coef[(mix_num+1):beta_surv_coef_len] <- unlist(altered_surv_coef)
  
  return(beta_surv_coef)
}

# beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef)

CalcBetaManual <- function(beta_surv_coef,surv_covar_risk_vec,stop_crit){
  l2norm <- 101
  while (l2norm > stop_crit){

    beta_vec <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len)[[1]]
    surv_coef <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len)[[2]]
    
    old_slike <- SurvLike(beta_vec,surv_covar_risk_vec,surv_coef)
    slike_diff <- -1


    grad <- numeric(mix_num + length(unlist(surv_coef)) - length(surv_coef))
    hess <- matrix(0,mix_num + length(unlist(surv_coef)) - length(surv_coef),mix_num + length(unlist(surv_coef)) - length(surv_coef))
    
    for (ind in which(surv_event == 1)){
      risk_set <- surv_time >= surv_time[ind]
      
      
      linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]
      ageadj_risk <- exp(linear_surv_covar_risk) %x% t(exp(beta_vec))
      
      num_list <- list()
      
      #First survival covariate must be age
      #Age is only continuous surv covar so its treated differently
      num0 <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set] %x% t(numeric(mix_num)+1))
      num02 <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set]^2 %x% t(numeric(mix_num)+1))
      num_list[[1]] <- c(num0,num02)
      
      denom <- sum(re_prob[risk_set,] * ageadj_risk)
      
      #fix age to be first variable
      grad[1] <- grad[1] + (surv_covar[[1]][ind] - num0/denom)
      hess[1,1] <- hess[1,1] + (num02/denom) - (num0/denom)^2
      
      

      for (beta_ind in 2:mix_num){
        num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk))
        grad[beta_ind] <- grad[beta_ind] + (re_prob[ind,beta_ind] - num/denom)
        
        hess[beta_ind,beta_ind] <- hess[beta_ind,beta_ind] + (num/denom) - (num/denom)^2

        num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[1]][risk_set])
        hess[1,beta_ind] <- hess[1,beta_ind] + (denom*num_cross - num*num0)/denom^2
        hess[beta_ind,1] <- hess[beta_ind,1] + (denom*num_cross - num*num0)/denom^2
      }
      
      if (length(surv_covar) > 1){
        for (surv_covar_ind in 2:length(surv_covar)){
          #need to only include people with same covariate status in numerator
          #num0 and num02 will be equal, need to adjust accordingly
          
          num_disc_covar_vec <- c()
          for (surv_covar_indicator in 2:length(surv_coef[[surv_covar_ind]])){
            
            num_disc_covar <- sum(re_prob[risk_set,] * ageadj_risk * 
                                    surv_covar[[surv_covar_ind]][risk_set,surv_covar_indicator] %x% t(numeric(mix_num)+1))
            
            num_disc_covar_vec <- c(num_disc_covar_vec,num_disc_covar)
            num_list[[surv_covar_ind]] <- num_disc_covar_vec
          }
        }
      }
      
      #Covariates
      list_of_lens <- unlist(lapply(num_list,length))
      list_of_lens[1] <- 0
      list_of_cum_lens <- cumsum(list_of_lens)
      for (surv_covar_ind in 2:length(surv_covar)){
        #Get indicies for covariates in gradient/hessian
        starting_index <- mix_num + list_of_cum_lens[surv_covar_ind-1] +1
        ending_index <- mix_num + list_of_cum_lens[surv_covar_ind]
        curr_covar_inds <- c(starting_index:ending_index)
        
        # grad[5] <- grad[5] + (surv_covar[[2]][ind,2] - num_list[[2]][2]/denom)
        grad[curr_covar_inds] <- grad[curr_covar_inds] + surv_covar[[surv_covar_ind]][ind,-1] - num_list[[surv_covar_ind]]/denom
        
        # hess[5,5] <- hess[5,5] + num_list[[2]][2]/denom - (num_list[[2]][2]/denom)^2
        diag(hess)[curr_covar_inds] <- diag(hess)[curr_covar_inds] + num_list[[surv_covar_ind]]/denom - (num_list[[surv_covar_ind]]/denom)^2
        
        #Cross terms for age
        for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
          # num_cross_age2 <- sum(re_prob[risk_set,] * ageadj_risk * age_vec[risk_set]*surv_covar[[2]][risk_set,2]%x% t(numeric(mix_num)+1))
          # hess[1,5] <- hess[1,5] + (denom*num_cross_age2 - num_list[[2]][[2]]*num0)/denom^2
          num_cross_age <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set]*surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1]%x% t(numeric(mix_num)+1))
          hess[1,curr_covar_inds[ind_curr_covar_inds]] <- hess[1,curr_covar_inds[ind_curr_covar_inds]] + (denom*num_cross_age - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num0)/denom^2
          hess[curr_covar_inds[ind_curr_covar_inds],1] <- hess[curr_covar_inds[ind_curr_covar_inds],1] + (denom*num_cross_age - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num0)/denom^2
        }
        
        #Cross terms for betas
        for (beta_ind in 2:mix_num){
          for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
            # num_cross2 <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[2]][risk_set,2])
            # hess[5,beta_ind] <- hess[5,beta_ind] + (denom*num_cross2 - num*num_list[[2]][[2]])/denom^2

            num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1])
            hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] <- hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] + 
              (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2
            
            hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] <- hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] + 
              (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2
          }
        }
        
        #Cross terms for surv covariates
        if(surv_covar_ind != length(surv_covar)){
          for (surv_covar_ind2 in (surv_covar_ind+1):length(surv_covar)){
            
            starting_index2 <- mix_num + list_of_cum_lens[surv_covar_ind2-1] +1
            ending_index2 <- mix_num + list_of_cum_lens[surv_covar_ind2]
            curr_covar_inds2 <- c(starting_index2:ending_index2)
            
            for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
              for (ind_curr_covar_inds2 in 1:length(curr_covar_inds2)){
                
                num_cross_covar <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1]*surv_covar[[surv_covar_ind2]][risk_set,ind_curr_covar_inds2+1]%x% t(numeric(mix_num)+1))
                
                hess[curr_covar_inds2[ind_curr_covar_inds2],curr_covar_inds[ind_curr_covar_inds]] <- hess[curr_covar_inds2[ind_curr_covar_inds2],curr_covar_inds[ind_curr_covar_inds]] + 
                  (denom*num_cross_covar - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]])/denom^2
                
                hess[curr_covar_inds[ind_curr_covar_inds],curr_covar_inds2[ind_curr_covar_inds2]] <- hess[curr_covar_inds[ind_curr_covar_inds],curr_covar_inds2[ind_curr_covar_inds2]] + 
                  (denom*num_cross_covar - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]])/denom^2
              }
            }
            
            
            
          }
        }
        
        
        
      }
      

    }

    step_size <- .01
    step_mat <- matrix(0,dim(hess)[1],dim(hess)[2])
    nr_fact <- numeric(length(grad))-1
    max_val <- 6
    while(slike_diff < 0|all(nr_fact == -1) | max_val > 5){
    # while(all(nr_fact == -1)){
      diag(step_mat) <- (diag(step_mat) + step_size - .01)
      nr_fact <- SolveCatch(hess + step_mat,-grad)
      step_size <- step_size * 10
    
      
      beta_surv_coef_new <- beta_surv_coef-nr_fact
      
      beta_vec_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len)[[1]]
      surv_coef_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len)[[2]]
      max_val <- abs(max(c(beta_vec_new,unlist(surv_coef_new))))
      # surv_coef_new <- list(beta_surv_coef_new[1])
      
      
      new_slike <- SurvLike(beta_vec_new,surv_covar_risk_vec,surv_coef_new)
      slike_diff <- old_slike - new_slike
    }

    l2norm <- sum(sqrt((beta_surv_coef_new - beta_surv_coef)^2))



    beta_surv_coef <- beta_surv_coef_new
    # print(slike_diff)

    # beta_vec_long[beta_vec_long>200] <- 200

  }
  # return(beta_surv_coef)
  return(list(beta_surv_coef,sqrt(diag(solve(hess)))))
}

RemFirCol <- function(x){return(x[,-1])}

CalcBeta <- function(beta_surv_coef, combined_covar_mat,surv_covar_risk_vec ,incl_surv, stop_crit = 1){
  
  if (incl_surv!=0){
    if (incl_surv == 2){stop_crit <- 100}
    beta_surv_coef_new <- CalcBetaManual(beta_surv_coef, surv_covar_risk_vec,stop_crit)
    return(beta_surv_coef_new)
  }
  
  surv_data <- data.frame(time = surv_time,
                          status = surv_event,
                          age = surv_covar[[1]])
  
  surv_data <- cbind(surv_data,re_prob,combined_covar_mat)
  colnames(surv_data)[4] <- "toRem"
  
  fit <- coxph(Surv(time, status) ~ .  - toRem, data = surv_data)
  beta_surv_coef_new <- fit$coefficients 
  se <- sqrt(diag(fit$var))
  
  return(list(beta_surv_coef_new,se))

  
}

Vec2Mat <- function(vect){
  mat <- matrix(0,nrow = length(vect),ncol = max(vect))
  
  for(i in 1:length(vect)){
    mat[i,vect[i]] <- 1
  }
  return (mat)
}

FirstDay2SingleDay <- function(first_day,target_day){
  
  day_to_keep_vec <- numeric(864)
  
  if (first_day == target_day){
    day_to_keep_vec[673:768] <- 1
  } else {
    if (first_day > target_day){
      day_ind <- 7 - (first_day-target_day)
    } else {
      day_ind <- target_day - first_day
    }
    first_day_ind <- (96 * (day_ind)) + 1
    last_day_ind <- first_day_ind + 95
    day_to_keep_vec[first_day_ind:last_day_ind] <- 1
  }
    
  return(day_to_keep_vec)
}


FirstDay2WeekInd <- function(first_day){
  weekday <- numeric(96)
  friday <- c(rep(0,68),rep(1,28))
  saturday <- numeric(96)+1
  sunday <- c(rep(1,68),rep(0,28))
  
  if (first_day == 1){
    covar_vec <- c(sunday,rep(weekday,4),friday,saturday,sunday,weekday)
  } else if (first_day == 2) {
    covar_vec <- c(rep(weekday,4),friday,saturday,sunday,rep(weekday,2))
  } else if (first_day == 3) {
    covar_vec <- c(rep(weekday,3),friday,saturday,sunday,rep(weekday,3))
  } else if (first_day == 4) {
    covar_vec <- c(rep(weekday,2),friday,saturday,sunday,rep(weekday,4))
  } else if (first_day == 5) {
    covar_vec <- c(weekday,friday,saturday,sunday,rep(weekday,4),friday)
  } else if (first_day == 6) {
    covar_vec <- c(friday,saturday,sunday,rep(weekday,4),friday,saturday)
  } else if (first_day == 7) {
    covar_vec <- c(saturday,sunday,rep(weekday,4),friday,saturday,sunday)
  }
  
  return(covar_vec)
}

ParamsArray2DF <- function(params_tran_array, misspecification = 0){
  
  tran_df <- data.frame(prob = c(),
                        type = c(),
                        time = c(),
                        age = c(),
                        weekend = c(),
                        mixture = c())
  
  
  for (re_ind in 1:dim(params_tran_array)[1]){
    for (vcovar_ind in 1:vcovar_num){
      
      params_tran <- params_tran_array[,,vcovar_ind]
      if (mix_num == 1 | misspecification == 1) {params_tran <- matrix(params_tran,nrow = 1)}
      
      tran_mat <- Params2TranVectorTresid(re_ind,96,params_tran)
      tosleep <- tran_mat[,3]
      towake <- tran_mat[,2]
      
      tran_df_working <- data.frame(prob = c(tosleep,towake),
                                    type = rep(c("Falling Asleep", "Waking"),each= 96),
                                    time = rep(c(1:96)/4,2),
                                    weekend = vcovar_ind,
                                    mixture = re_ind)
      tran_df <- rbind(tran_df,tran_df_working)
    }
  }
  
  return(tran_df)
}

CalcPiHelper <- function(nu_mat,nu_covar_vec){
  pi_ind <- exp(colSums(nu_mat * nu_covar_vec))/sum(exp(colSums(nu_mat * nu_covar_vec)))
  if (any(is.na(pi_ind))){
    pi_ind[is.na(pi_ind)] <- 1
    pi_ind <- pi_ind/sum(pi_ind)
  }
  return(pi_ind)
}

CalcPi <- function(nu_mat,nu_covar_mat){
  
  pi_l_new <- matrix(NA,dim(nu_covar_mat)[1],dim(nu_mat)[2])
  for (ind in 1:dim(nu_covar_mat)[1]){
    pi_l_new[ind,] <- CalcPiHelper(nu_mat,nu_covar_mat[ind,])
  }
  
  return(pi_l_new)
}

CalcNu <- function(nu_mat,re_prob,nu_covar_mat){
  if(dim(re_prob)[2] == 1){return(nu_mat)}
  # nu_long <- as.vector(rbind(nu[-1],nu2[-1]))
  # old_mlike <- NuLogLike(nu_long)
  old_mlike <- CalcLikelihood(alpha,CalcPi(nu_mat,nu_covar_mat))
  mlike_diff <- -1
  nnu_covar <- dim(nu_covar_mat)[2]
  gradient_nu <- numeric(mix_num * nnu_covar)
  hess_nu <- matrix(0,mix_num*nnu_covar,mix_num*nnu_covar)
  
  for (ind in 1:num_of_people){
    age_ind_vec <- nu_covar_mat[ind,]
    age_ind_mat <- age_ind_vec %*% t(age_ind_vec)
    
    num <- exp(colSums(nu_mat * nu_covar_mat[ind,]))
    denom <- sum(num)
    p_vec <- num/denom
    
    pvec_grad <- (re_prob[ind,] - p_vec)
    pvec_grad[1] <- 0
    gradient_nu <- gradient_nu + pvec_grad %x% age_ind_vec
    
    # p_vec[1]*(1-p_vec[1]) * age_ind_mat
    
    p_vec[1] <- 0
    p_mat <- p_vec %x% t(p_vec)
    diag(p_mat) <- -p_vec * (1-p_vec)
    hess_nu <- hess_nu + p_mat %x% age_ind_mat
  }
  
  hess_nu <- hess_nu[c(-(1:nnu_covar)),c(-(1:nnu_covar))]
  gradient_nu <- gradient_nu[c(-(1:nnu_covar))]
  
  inf_mat <- numeric(length(gradient_nu))-1
  
  step_size <- .01
  while(mlike_diff < 0 | all(inf_mat == -1)){
  # while(all(inf_mat == -1)){
    
    step_fact <- matrix(0,dim(hess_nu)[1],dim(hess_nu)[2])
    
    diag(step_fact) <- step_size - .01
    
    inf_mat <- SolveCatch(hess_nu+step_fact,gradient_nu)
    step_size <- step_size * 10
    
    
    nu_mat_lm <-cbind(rep(0,nnu_covar),matrix(inf_mat,nrow = nnu_covar,byrow = F))
    nu_mat_new <- nu_mat - nu_mat_lm
    
    pi_l_new <- CalcPi(nu_mat_new,nu_covar_mat)
    new_mlike <- CalcLikelihood(alpha,pi_l_new)
    mlike_diff <- new_mlike - old_mlike
    
  }
  # print(mlike_diff)
  
  
  
  return(nu_mat_new)
}

# NuLogLike <- function(nu_long){
#   nu_mat <- matrix(nu_long,nrow = 2, byrow = F)
#   nu_mat <- cbind(c(0,0),nu_mat)
#   nu <- nu_mat[1,]
#   nu2 <- nu_mat[2,]
#   age_vec <- age_vec/10
#   loglike <- 0
#   for (ind in 1:num_of_people){
#     age_ind_vec <- c(age_vec[ind],age_vec[ind]^2)
#     for (re_ind in 1:mix_num){
#       cont <- (age_ind_vec %*% nu_mat[,re_ind]) - log(sum(exp(age_ind_vec %*% nu_mat)))
#       loglike_ind <- re_prob[ind,re_ind]*cont
#       loglike <- loglike + loglike_ind
#     }
#   }
#   return(loglike)
# }
# 
# CalcNuNum <- function(nu,nu2){
#   nu_long <- c(nu[-1],nu2[-1])
#   grad_nu_num <- grad(NuLogLike,nu_long)  
#   hess_nu_num <- hessian(NuLogLike,nu_long) 
#   
#   inf_mat <- solve(hess_nu_num,grad_nu_num)
#   
#   step_size <- 1.1
#   if (max(abs(inf_mat)) > .25){
#     diag(hess_nu_num) <- diag(hess_nu_num) * step_size
#     inf_mat <- solve(hess_nu_num,grad_nu_num)
#     step_size <- step_size * 10
#   }
#   
#   
#   nu_long_new <- nu_long - inf_mat
#   nu_mat <- matrix(nu_long_new,nrow = 2, byrow = T)
#   # nu_mat[,1] <- 0
#   nu <- c(0,nu_mat[1,])
#   nu2 <- c(0,nu_mat[2,])
#   return(list(nu,nu2))
# }


Viterbi <- function(act,light,vcovar_mat){
  decoded_array <- array(NA, dim = c(day_length,num_of_people,mix_num))
  for (ind in 1:num_of_people){
    for (clust_i in 1:mix_num){
      vit_ind_vec <- ViterbiIndHelper(ind,clust_i,act,light,vcovar_mat)
      decoded_array[,ind,clust_i] <- vit_ind_vec
    }
  }
  return(decoded_array)
}

ViterbiIndHelper <- function(ind,clust_i,act,light,vcovar_mat){
  
  
  tran_list_clust <- tran_list[[clust_i]]
  
  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))
  
  vcovar_vec <- vcovar_mat[,ind]
  
  
  log_class_0_week <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_week[1,1,clust_i], 
                                          emit_act_week[1,2,clust_i], 
                                          emit_light_week[1,1,clust_i],
                                          emit_light_week[1,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,1,1], lintegral_mat[clust_i,1,1],
                                          lambda_act_mat[clust_i,1,1],lambda_light_mat[clust_i,1,1],tobit)
  
  log_class_1_week <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_week[2,1,clust_i], 
                                          emit_act_week[2,2,clust_i], 
                                          emit_light_week[2,1,clust_i],
                                          emit_light_week[2,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,2,1], lintegral_mat[clust_i,2,1],
                                          lambda_act_mat[clust_i,2,1],lambda_light_mat[clust_i,2,1],tobit)
  
  log_class_0_weekend <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_weekend[1,1,clust_i], 
                                          emit_act_weekend[1,2,clust_i], 
                                          emit_light_weekend[1,1,clust_i],
                                          emit_light_weekend[1,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,1,2], lintegral_mat[clust_i,1,2],
                                          lambda_act_mat[clust_i,1,2],lambda_light_mat[clust_i,1,2],tobit)
  
  log_class_1_weekend <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_weekend[2,1,clust_i], 
                                          emit_act_weekend[2,2,clust_i], 
                                          emit_light_weekend[2,1,clust_i],
                                          emit_light_weekend[2,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,2,2], lintegral_mat[clust_i,2,2],
                                          lambda_act_mat[clust_i,2,2],lambda_light_mat[clust_i,2,2],tobit)
  
  
  log_class_0 <- (log_class_0_week * (1-vcovar_vec)) + (log_class_0_weekend * vcovar_vec) 
  log_class_1 <- (log_class_1_week * (1-vcovar_vec)) + (log_class_1_weekend * vcovar_vec)
  
  
  viterbi_mat <- matrix(NA,2,day_length)
  viterbi_mat[1,1] <- log(init[clust_i,1]) + log_class_0[1]
  viterbi_mat[2,1] <- log(init[clust_i,2]) + log_class_1[1]
  
  viterbi_ind_mat <- matrix(NA,2,day_length)
  
  
  for (time in 2:day_length){
    
    tran <- tran_list_clust[[vcovar_vec[time]+1]][[time]]
    
    viterbi_mat[1,time] <- log_class_0[[time]] + 
      max(viterbi_mat[1,time-1] + log(tran[1,1]),
          viterbi_mat[2,time-1] + log(tran[2,1]))
    
    
    viterbi_mat[2,time] <- log_class_1[[time]] + 
      max(viterbi_mat[1,time-1] + log(tran[1,2]),
          viterbi_mat[2,time-1] + log(tran[2,2]))
    
    
    viterbi_ind_mat[1,time] <-  which.max(c(viterbi_mat[1,time-1] + log(tran[1,1]),
                                            viterbi_mat[2,time-1] + log(tran[2,1])))
    
    
    viterbi_ind_mat[2,time] <- which.max(c(viterbi_mat[1,time-1] + log(tran[1,2]),
                                           viterbi_mat[2,time-1] + log(tran[2,2])))
    
    
  }
  
  decoded_mc <- c(which.max(viterbi_mat[,time]))
  for(time in day_length:2){
    decoded_mc <- c(viterbi_ind_mat[decoded_mc[1],time],decoded_mc)
  }
  
  return(decoded_mc-1)
} 


CalcBinom <- function(act,light,weights_array_list){
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])
  
  act_vec <- as.vector(act)
  light_vec <- as.vector(light)
  vcovar_vec <- as.vector(vcovar_mat)
  
  lod_act_weight <- as.numeric(act_vec==lod_act)
  lod_light_weight <- as.numeric(light_vec==lod_light)
  
  
  lambda_act_mat <- array(NA,dim = c(mix_num,2,2))
  lambda_light_mat <- array(NA,dim = c(mix_num,2,2))
  
  for (clust_i in 1:mix_num){
    weights_vec_wake <- as.vector(weights_array_wake[,,clust_i])
    weights_vec_sleep <- as.vector(weights_array_sleep[,,clust_i])
    for (vcovar_ind in c(0,1)){
      
      vcovar_bool_vec <- vcovar_vec == vcovar_ind
      
      lambda_act_mat[clust_i,1,vcovar_ind+1] <- sum(lod_act_weight[vcovar_bool_vec] *weights_vec_wake[vcovar_bool_vec],na.rm=T)/
        sum(weights_vec_wake[vcovar_bool_vec][!is.na(act_vec[vcovar_bool_vec])])
      
      lambda_act_mat[clust_i,2,vcovar_ind+1] <- sum(lod_act_weight[vcovar_bool_vec] *weights_vec_sleep[vcovar_bool_vec],na.rm=T)/
        sum(weights_vec_sleep[vcovar_bool_vec][!is.na(act_vec[vcovar_bool_vec])])
      
      
      
      lambda_light_mat[clust_i,1,vcovar_ind+1] <- sum(lod_light_weight[vcovar_bool_vec] *weights_vec_wake[vcovar_bool_vec],na.rm=T)/
        sum(weights_vec_wake[vcovar_bool_vec][!is.na(light_vec[vcovar_bool_vec])])
      
      lambda_light_mat[clust_i,2,vcovar_ind+1] <- sum(lod_light_weight[vcovar_bool_vec] *weights_vec_sleep[vcovar_bool_vec],na.rm=T)/
        sum(weights_vec_sleep[vcovar_bool_vec][!is.na(light_vec[vcovar_bool_vec])])
    }
  
  }

  
  return(list(lambda_act_mat,lambda_light_mat))
}


CalcMeans <- function(obs_mat,lod_obs,weights_array){
  n <- dim(weights_array)[2]
  mean_vec <- matrix(0,mix_num,2)
  
  leps_indicator <- (obs_mat==lod_obs)
  leps_indicator[is.na(leps_indicator)] <- F
  
  for (clust_i in 1:mix_num){
    for (vcovar_ind in c(0,1)){
      num <- 0
      denom <- 0
      for (ind in 1:n){
        obs_ind0 <- obs_mat[,ind] 
        vcovar_vec <- vcovar_mat[,ind] == vcovar_ind
        
        leps_ind <- leps_indicator[,ind]
        inds_keep <- (!is.na(obs_ind0) & !leps_ind) & vcovar_vec
        
        num <- num + (weights_array[inds_keep,ind,clust_i]) %*% obs_ind0[inds_keep]
        denom <- denom + sum(weights_array[inds_keep,ind,clust_i])
      }
      
      mean_vec[clust_i,vcovar_ind+1] <- num/denom
    }
   
    

  }
  return(mean_vec)
}


CalcSigmas <- function(obs_mat,lod_obs,weights_array,mean_mat){
  n <- dim(weights_array)[2]
  sig_mat <- matrix(0,mix_num,2)
  
  leps_indicator <- (obs_mat==lod_obs)
  leps_indicator[is.na(leps_indicator)] <- F
  
  for (clust_i in 1:mix_num){
    for (vcovar_ind in c(0,1)){
      num <- 0
      denom <- 0
      for (ind in 1:n){
        resid0 <- (obs_mat[,ind]-mean_mat[clust_i,vcovar_ind+1])^2
        
        vcovar_vec <- vcovar_mat[,ind] == vcovar_ind
        leps_ind <- leps_indicator[,ind]
        inds_keep <- (!is.na(resid0) & !leps_ind) & vcovar_vec
        
        num <- num + (weights_array[inds_keep,ind,clust_i]) %*% resid0[inds_keep]
        denom <- denom + sum(weights_array[inds_keep,ind,clust_i])
      }
      
      sig_mat[clust_i,vcovar_ind+1] <- sqrt(num/denom)
    }
    
    
    
  }
  return(sig_mat)
}


Forward <- function(act, light,init,tran_list, 
                    emit_act, emit_light,
                    lod_act, lod_light, corr_mat, beta_vec, surv_coef, surv_covar_risk_vec,
                    event_vec, bline_vec, cbline_vec, lintegral_mat, log_sweights_vec, 
                    surv_covar, vcovar_mat, lambda_act_mat, lambda_light_mat, tobit, incl_surv){
  
  alpha_list <- list()
  day_length <- dim(act)[1]
  
  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))
  
  if (incl_surv == 0 | incl_surv == 1){
    adj_incl_surv <- 0
  } else {
    adj_incl_surv <- 1
  }
  
  # log(bline_vec[ind]) + beta_vec[clust_i]+covar_risk - (cbline_vec[ind] * exp(beta_vec[clust_i]+(covar_risk)))
  
  for (ind in 1:dim(act)[2]){
    alpha_array <- array(NA,dim = c(day_length,2,mix_num))
    ind_log_sweight <- log_sweights_vec[ind]
    covar_risk <- surv_covar_risk_vec[ind]
    
    for (clust_i in 1:mix_num){
      alpha_array[,,clust_i] <- ForwardIndC(act[,ind], light[,ind], init[clust_i,], tran_list, emit_act_week, emit_light_week, 
                                            emit_act_weekend, emit_light_weekend,clust_i-1, lod_act, lod_light, corr_mat, 
                                            beta_vec,covar_risk,surv_event[ind], bline_vec[ind], cbline_vec[ind], lintegral_mat, 
                                            ind_log_sweight, vcovar_mat[,ind],lambda_act_mat, lambda_light_mat, tobit, adj_incl_surv*beta_bool)
      alpha_list[[ind]] <- alpha_array
    }
  }
  return(alpha_list)
}


Backward <- function(act, light, tran_list, 
                     emit_act, emit_light, 
                     lod_act, lod_light, corr_mat, lintegral_mat, vcovar_mat,
                     lambda_act_mat, lambda_light_mat, tobit){
  
  beta_list <- list()
  
  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))
  
  for (ind in 1:dim(act)[2]){
    beta_array <- array(NA,dim = c(day_length,2,mix_num))
    
    
    for (clust_i in 1:mix_num){
      #NEED TO ADD WEIGHT
      beta_array[,,clust_i] <- BackwardIndC(act[,ind], light[,ind], tran_list, emit_act_week, emit_light_week, 
                                            emit_act_weekend, emit_light_weekend,clust_i-1, lod_act, lod_light, corr_mat, 
                                            lintegral_mat, 
                                            vcovar_mat[,ind],lambda_act_mat, lambda_light_mat, tobit)
      beta_list[[ind]] <- beta_array
    }
  }
  return(beta_list)
  
}


CalcS <- function(event_time,cbline_vec_new,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment,surv_covar_risk_vec){
  surv_mat_ind <- matrix(NA,length(event_time),dim(re_prob)[1])
  surv_vec <- numeric(length(event_time))
  
  for (ind in 1:dim(re_prob)[1]){
    for (t in 1:length(event_time)){
      rs <- 0
      
      
      for (re_ind in 1:mix_num){
        rs <- rs + exp(-cbline_vec_new[t]* exp(beta_vec[re_ind]+surv_covar_risk_vec[ind])) * re_prob[ind,re_ind]
      }
      
      surv_vec[t] <- rs
    }
    surv_mat_ind[,ind] <- surv_vec
  }
  
  return(surv_mat_ind)
}

Vec2StepPlot <- function(cens_dist,cens_dist_time,t_i){
  for (i in 1:(length(cens_dist_time)-1)){
    if ((t_i >= cens_dist_time[i]) & (t_i < cens_dist_time[i+1])){
      return(cens_dist[i])
    }
  }
  return(cens_dist[length(cens_dist_time)])
  
}

BrierScore <- function(bs_t,surv_event,surv_time,cens_dist,cens_dist_time,surv_mat_ind,disc){
  bs_sum <- 0
  
  if (disc){
    to_add <- 1
  } else {
    to_add <- .00000000001
  }
  
  for (ind in 1:num_of_people){
    
    sprob <- Vec2StepPlot(surv_mat_ind[,ind],event_time,bs_t)
    
    num1 <- sprob^2 * (surv_time[ind] <= bs_t) * surv_event[ind]
    num2 <- (1-sprob)^2 * (surv_time[ind] > bs_t) 
    
    denom1 <- Vec2StepPlot(cens_dist,cens_dist_time,surv_time[ind]-to_add)
    denom2 <-  Vec2StepPlot(cens_dist,cens_dist_time,bs_t)
    
    bs_sum <- bs_sum + (num1/denom1) + (num2/denom2)
    if (is.na(bs_sum)){break}
  }
  return(bs_sum/num_of_people)
}

vBrierScore <- Vectorize(BrierScore,vectorize.args = "bs_t")

IntBS <- function(vbs,surv_event,surv_time,cens_dist,cens_dist_time,surv_mat_ind,disc){
  # vbs <- vBrierScore(unique(sort(surv_time)),surv_event,surv_time,cens_dist,cens_dist_time,surv_mat_ind,disc)
  sorted_surv_time <- unique(sort(surv_time))
  
  riem_sum <- 0
  for (b_t in 2:length(sorted_surv_time)){
    riem_sum <- riem_sum + (vbs[b_t-1] * (sorted_surv_time[b_t]-sorted_surv_time[b_t-1]))
  }
  IBS <- riem_sum/max(surv_time)
  
  return(IBS)
}

IndBSold <- function(){
  bs_vec <- numeric(num_of_people)
  for (bs_ind in 1:num_of_people){
    bs_t <- unique(sort(surv_time))[bs_ind]
    ind <- 1
    bssum <- 0
    
    
    
    for (ind in 1:num_of_people){
      f1 <- pred_surv[ind,bs_ind]^2 * (surv_time[ind] <= bs_t) * surv_event[ind]/cens_dist[which(cens_dist_time >= surv_time[ind])[1]-1]
      f2 <- (1-pred_surv[ind,bs_ind])^2 * (surv_time[ind] > bs_t) /cens_dist[which(cens_dist_time > bs_t)[1]-1]
      
      f1_vec[ind] <- f1
      f2_vec[ind] <- f2
      
      if (is.na(f1 + f2)){break}
      bssum <- bssum + f1 + f2
    }
    
    bs_vec[bs_ind] <- bssum/num_of_people
  }
  return(bs_vec)
}


CalcCindex <- function(surv_time,surv_event,beta_vec,surv_coef,re_prob,surv_covar,surv_covar_risk_vec){
  denom <- 0
  num <- 0
  for (j in 1:length(surv_time)){
    
    
    for (i in 1:length(surv_time)){
      if (surv_event[j] == 1){
        if (i != j){
          
          risk_i <- (re_prob[i,] %*% beta_vec)+surv_covar_risk_vec[i]
          risk_j <- (re_prob[j,] %*% beta_vec)+surv_covar_risk_vec[j]
          
          if (surv_time[i] > surv_time[j]){
            denom <- denom + 1
            
            if (risk_j > risk_i){
              num <- num + 1
            }
            
          }
          
          
        }
      }
      
    }
  }
  
  return(num/denom)
}


CalcIBS <- function(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment){
  event_time <- unique(sort(surv_time))
  cbline_vec_new <- unique(sort(cbline_vec))
  # cbline_vec_new <- basehaz(cox1,centered = F)[,1]
  
  stime_vec <- c()
  for (stime in event_time){
    stime_vec <- c(stime_vec,which(surv_time == stime)[1])
  }
  
  cbline_vec_new <- cbline_vec[stime_vec]
  
  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment,surv_covar_risk_vec_new)
  
  
  
  IBS_times <- sort(unique(surv_time))
  # IBS_times <- sort(surv_time)
  # Calculate the IBS, again including the observed event times, censoring
  # variables, and the prediction timepoints corresponding to each row of the
  # survival prob matrix
  
  ibs <- integrated_brier_score(y_true = Surv(surv_time, surv_event), 
                                surv = t(surv_mat_ind), 
                                times = IBS_times)
  
  return(ibs)
}

SurvCovar2Coef <- function(covar_mat,max_val = .1){
  return(seq(0,max_val,length.out = dim(covar_mat)[2]))
}

readCpp( "Scripting/cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )
################## EM Setup ################## 
###### True Settings ###### 
if (sim_size == 0){
  day_length <- 96 * 2
  num_of_people <- 1000
  missing_perc <- .2
} else if (sim_size== 1) {
  day_length <- 96 * 4
  num_of_people <- 4000
  missing_perc <- 0
} else if (sim_size== 2) {
  day_length <- 96 * 7
  num_of_people <- 7000
  missing_perc <- 0
}
  


if (misspecification){
  model_name_loadin <- "JMHMM"
  if (!tobit){model_name_loadin <- paste0(model_name_loadin,"NonTob")}
  if (incl_surv==0){model_name_loadin <- paste0(model_name_loadin,"NoSurv")}
  if (incl_surv==1){model_name_loadin <- paste0(model_name_loadin,"HalfSurv")}
  if (!incl_light){model_name_loadin <- paste0(model_name_loadin,"NoLight")}
  if (!incl_act){model_name_loadin <- paste0(model_name_loadin,"NoAct")}
  model_name_loadin <- paste0(model_name_loadin,"Mix",misspecification,"Seed",".rda")
  
  print(paste("Loading",model_name_loadin))
  setwd("Data")
  load(model_name_loadin)
  setwd("..")
  
  
  init_true <- to_save[[2]][[1]]
  params_tran_array_true <- to_save[[2]][[2]]
  emit_act_true <- to_save[[2]][[3]]
  emit_light_true <- to_save[[2]][[4]]
  corr_mat_true <- to_save[[2]][[5]]
  nu_mat_true <- to_save[[2]][[6]]
  beta_vec_true <- to_save[[2]][[7]]
  surv_coef_true <- to_save[[2]][[8]]
  re_prob_true <- to_save[[2]][[10]]
  #NEED TO CHANGE THIS TOO
  #FIRST 13 SHOULD BE 12
  lambda_act_mat_true <- to_save[[2]][[13]]
  lambda_light_mat_true <- to_save[[2]][[13]]
  
  if (misspecification == mix_num){
    init <- init_true
    params_tran_array <- params_tran_array_true
    emit_act <- emit_act_true
    emit_light <- emit_light_true
    corr_mat <- corr_mat_true
    beta_vec <- beta_vec_true
    beta_age <- beta_age_true
    nu_mat <- nu_mat_true
    lambda_act_mat <- lambda_act_mat_true
    lambda_light_mat <- lambda_light_mat_true
  } else {
    
    params_tran_week <- matrix(rep(c(0,0,0,0,0,0),mix_num),ncol = 6,byrow = T)
    params_tran_week[1:misspecification,] <- params_tran_array_true[,,1]
    params_tran_week[(misspecification+1):mix_num,] <- params_tran_array_true[1,,1]
    
    params_tran_weekend <- matrix(rep(c(0,0,0,0,0,0),mix_num),ncol = 6,byrow = T)
    params_tran_weekend[1:misspecification,] <- params_tran_array_true[,,2]
    params_tran_weekend[(misspecification+1):mix_num,] <- params_tran_array_true[1,,2]
    
    params_tran_array_dim <- c(mix_num,6,vcovar_num)
    params_tran_array <- array(NA,dim = params_tran_array_dim)
    params_tran_array[,,1] <- params_tran_week 
    params_tran_array[,,2] <- params_tran_weekend
    params_tran_array[(misspecification+1):mix_num,,] <- params_tran_array[(misspecification+1):mix_num,,] + 
      runif(length(unlist(params_tran_array[(misspecification+1):mix_num,,])),-.5,.5)
    
    emit_act_week <- array(emit_act_true[,,,1], c(2,2,mix_num))
    emit_act_weekend <- array(emit_act_true[,,1,2], c(2,2,mix_num))
    emit_act <- array(NA, c(2,2,mix_num,2))
    emit_act[,,,1] <- emit_act_week
    emit_act[,,,2] <- emit_act_weekend
    emit_act[,,(misspecification+1):mix_num,] <- emit_act[,,(misspecification+1):mix_num,] + runif(4*(mix_num-misspecification),-.5,.5)
        
    emit_light_week <- array(emit_light_true[,,1,1], c(2,2,mix_num))
    emit_light_weekend <- array(emit_light_true[,,1,2], c(2,2,mix_num))
    emit_light <- array(NA, c(2,2,mix_num,2))
    emit_light[,,,1] <- emit_light_week
    emit_light[,,,2] <- emit_light_weekend
    emit_light[,,(misspecification+1):mix_num,] <- emit_light[,,(misspecification+1):mix_num,] + runif(4*(mix_num-misspecification),-.5,.5)
    
    corr_mat <- array(NA, c(mix_num,2,2))
    corr_mat[,1,1] <- corr_mat_true[1,1,1]
    corr_mat[,2,1] <- corr_mat_true[1,2,1]
    corr_mat[,1,2] <- corr_mat_true[1,1,2]
    corr_mat[,2,2] <- corr_mat_true[1,2,2]
    corr_mat[(misspecification+1):mix_num,,] <- corr_mat[(misspecification+1):mix_num,,] + runif(2*(mix_num-misspecification),-.15,.15)
    
    beta_vec <- seq(0,3,length.out = mix_num)
    beta_age <- 0.07 + runif(1,-.01,.01)
    
    nu <- c(0,seq(-.04,.05,length.out = (mix_num-1)))
    nu2 <- c(0,seq(.0005,-.0015,length.out = (mix_num-1)))
    nu_stat <- c(0,seq(.01,-.025,length.out = (mix_num-1)))
    nu2_stat <- c(0,seq(-.001,.002,length.out = (mix_num-1)))
    
    nu_mat <- matrix(NA,nrow = 4, ncol = mix_num)
    nu_mat[1,] <- nu
    nu_mat[2,] <- nu2
    nu_mat[3,] <- nu_stat
    nu_mat[4,] <- nu2_stat
    
    ###NEED TO UPDATE IF USING THIS
    lambda_act_mat <- array(NA,dim = c(mix_num,2,2))
    lambda_act_mat[,1,1] <- seq(.01,.05,length.out = mix_num)
    lambda_act_mat[,2,1] <- seq(.3,.6,length.out = mix_num)
    lambda_act_mat[,1,2] <- seq(.01,.1,length.out = mix_num)
    lambda_act_mat[,2,2] <- seq(.3,.7,length.out = mix_num)
    
    lambda_light_mat <- array(NA,dim = c(mix_num,2,2))
    lambda_light_mat[,1,1] <- seq(.01,.15,length.out = mix_num)
    lambda_light_mat[,2,1] <- seq(.2,.6,length.out = mix_num)
    lambda_light_mat[,1,2] <- seq(.01,.2,length.out = mix_num)
    lambda_light_mat[,2,2] <- seq(.3,.8,length.out = mix_num)
  }
  
} else {
  
  init_true <- matrix(NA,ncol = 2,nrow = mix_num)
  init_true[,1] <- seq(.1,.9,length.out = mix_num)
  init_true[,2] <- 1 - init_true[,1]
  
  params_tran_week_true <- matrix(rep(c(0,0,0,0,0,0),mix_num),ncol = 6,byrow = T) 
  params_tran_week_true[,1] <- seq(-3.2,-2,length.out = mix_num)
  params_tran_week_true[,2] <- seq(1.6,.8,length.out = mix_num)
  params_tran_week_true[,3] <- seq(.1,.9,length.out = mix_num)
  params_tran_week_true[,4] <- seq(-1.6,-2.4,length.out = mix_num)
  params_tran_week_true[,5] <- seq(-1.7,-.8,length.out = mix_num)
  params_tran_week_true[,6] <- seq(-.3,-.9,length.out = mix_num)
  
  params_tran_weekend_true <- matrix(rep(c(0,0,0,0,0,0),mix_num),ncol = 6,byrow = T) 
  params_tran_weekend_true[,1] <- seq(-3.2,-2.1,length.out = mix_num)
  params_tran_weekend_true[,2] <- seq(1,.5,length.out = mix_num)
  params_tran_weekend_true[,3] <- seq(.1,.9,length.out = mix_num)
  params_tran_weekend_true[,4] <- seq(-1.2,-2.4,length.out = mix_num)
  params_tran_weekend_true[,5] <- seq(-1.6,-.8,length.out = mix_num)
  params_tran_weekend_true[,6] <- seq(-.4,-1.1,length.out = mix_num)
  
  
  params_tran_array_dim <- c(mix_num,6,vcovar_num)
  params_tran_array_true <- array(NA,dim = params_tran_array_dim)
  params_tran_array_true[,,1] <- params_tran_week_true 
  params_tran_array_true[,,2] <- params_tran_weekend_true 
  
  emit_act_week_true <- array(NA, c(2,2,mix_num))
  emit_act_week_true[1,1,] <- seq(5,7,length.out = mix_num)
  emit_act_week_true[1,2,] <- seq(2,8,length.out = mix_num)
  emit_act_week_true[2,1,] <- seq(3,1,length.out = mix_num)
  emit_act_week_true[2,2,] <- seq(3,1,length.out = mix_num)
  
  emit_act_weekend_true <- array(NA, c(2,2,mix_num))
  emit_act_weekend_true[1,1,] <- seq(5,7,length.out = mix_num)
  emit_act_weekend_true[1,2,] <- seq(1,4,length.out = mix_num)
  emit_act_weekend_true[2,1,] <- seq(3,0,length.out = mix_num)
  emit_act_weekend_true[2,2,] <- seq(5,2,length.out = mix_num)
  
  emit_act_true <- array(NA, c(2,2,mix_num,2))
  emit_act_true[,,,1] <- emit_act_week_true
  emit_act_true[,,,2] <- emit_act_weekend_true
  
  
  emit_light_week_true <- array(NA, c(2,2,mix_num))
  emit_light_week_true[1,1,] <- seq(-5,-1,length.out = mix_num)
  emit_light_week_true[1,2,] <- seq(8,3,length.out = mix_num)
  emit_light_week_true[2,1,] <- seq(-15,-19,length.out = mix_num)
  emit_light_week_true[2,2,] <- seq(13,19,length.out = mix_num)
  
  emit_light_weekend_true <- array(NA, c(2,2,mix_num))
  emit_light_weekend_true[1,1,] <- seq(-5,-2,length.out = mix_num)
  emit_light_weekend_true[1,2,] <- seq(12,7,length.out = mix_num)
  emit_light_weekend_true[2,1,] <- seq(-15,-17,length.out = mix_num)
  emit_light_weekend_true[2,2,] <- seq(13,19,length.out = mix_num)
  
  emit_light_true <- array(NA, c(2,2,mix_num,2))
  emit_light_true[,,,1] <- emit_light_week_true
  emit_light_true[,,,2] <- emit_light_weekend_true
  
  corr_mat_true <- array(NA, c(mix_num,2,2))
  corr_mat_true[,1,1] <- seq(0,.5,length.out = mix_num)
  corr_mat_true[,2,1] <- seq(-.2,.3,length.out = mix_num)
  corr_mat_true[,1,2] <- seq(-.5,0,length.out = mix_num)
  corr_mat_true[,2,2] <- seq(-.5,.5,length.out = mix_num)
  
  beta_vec_true <- seq(0,3,length.out = mix_num)
  beta_age_true <- 0.07
  
  nu_true <- c(0,seq(-.04,.05,length.out = (mix_num-1)))
  nu2_true <- c(0,seq(.0005,-.0015,length.out = (mix_num-1)))
  nu_stat_true <- c(0,seq(.01,-.025,length.out = (mix_num-1)))
  nu2_stat_true <- c(0,seq(-.001,.002,length.out = (mix_num-1)))
  
  nu_mat_true <- matrix(NA,nrow = 4, ncol = mix_num)
  nu_mat_true[1,] <- nu_true
  nu_mat_true[2,] <- nu2_true
  nu_mat_true[3,] <- nu_stat_true
  nu_mat_true[4,] <- nu2_stat_true
  
  lambda_act_mat_true <- array(NA,dim = c(mix_num,2,2))
  lambda_act_mat_true[,1,1] <- seq(.01,.05,length.out = mix_num)
  lambda_act_mat_true[,2,1] <- seq(.3,.6,length.out = mix_num)
  lambda_act_mat_true[,1,2] <- seq(.01,.1,length.out = mix_num)
  lambda_act_mat_true[,2,2] <- seq(.3,.7,length.out = mix_num)
  
  lambda_light_mat_true <- array(NA,dim = c(mix_num,2,2))
  lambda_light_mat_true[,1,1] <- seq(.01,.15,length.out = mix_num)
  lambda_light_mat_true[,2,1] <- seq(.2,.6,length.out = mix_num)
  lambda_light_mat_true[,1,2] <- seq(.01,.2,length.out = mix_num)
  lambda_light_mat_true[,2,2] <- seq(.3,.8,length.out = mix_num)
  
}

  
if (load_data){
  # setwd("Data")
  # load("JMHMMMix4Seed.rda")
  # setwd("..")
  
  load(paste0("Inter",model_name))
  
  init_true <- to_save[[2]][[1]]
  params_tran_array_true <- to_save[[2]][[2]]
  emit_act_true <- to_save[[2]][[3]]
  emit_light_true <- to_save[[2]][[4]]
  corr_mat_true <- to_save[[2]][[5]]
  nu_mat_true <- to_save[[2]][[6]]
  beta_vec_true <- to_save[[2]][[7]]
  surv_coef_true <- to_save[[2]][[8]]
  #check this
  # re_prob_true <- to_save[[2]][[10]]
  #NEED TO CHANGE THIS TOO
  #FIRST 13 SHOULD BE 12
  # lambda_act_mat_true <- to_save[[2]][[13]]
  # lambda_light_mat_true <- to_save[[2]][[13]]
}

###### Simulate Data ###### 
if (!real_data){
  lod_act_true <- -5.809153
  lod_light_true <- -1.560658
  
  lod_act <- lod_act_true
  lod_light <- lod_light_true
  
  beta_covar_sim <- c(0,.6,-.5)
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                               init=init_true,params_tran_array = params_tran_array_true,
                               emit_act = emit_act_true,emit_light = emit_light_true,
                               corr_mat = corr_mat_true,
                               lod_act = lod_act_true,lod_light = lod_light_true,
                               nu_mat = nu_mat_true,
                               beta_age_true = beta_age_true,beta_covar_sim = beta_covar_sim,
                               missing_perc = missing_perc, beta_vec_true = beta_vec_true,
                               lambda_act_mat = lambda_act_mat_true,lambda_light_mat = lambda_light_mat_true,
                               misspecification = misspecification)
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  light <- simulated_hmm[[3]]
  mixture_mat <- simulated_hmm[[4]]
  age_vec <- simulated_hmm[[5]]
  nu_covar_mat <- simulated_hmm[[6]]
  vcovar_mat <-  simulated_hmm[[7]]
  # emp_params <- simulated_hmm[[7]]
  surv_list <- simulated_hmm[[8]]
  surv_covar_sim <- simulated_hmm[[9]]
  
  id_sim <- cbind(age_vec,surv_covar_sim-1)
  surv_covar <- list(age_vec,Vec2Mat(surv_covar_sim))
  surv_coef <- list(beta_age_true,beta_covar_sim)
  combined_covar_mat <- matrix(surv_covar_sim-1,nrow = num_of_people)
  combined_covar_mat <- as.factor(combined_covar_mat)
  
  surv_time <- surv_list[[1]]
  surv_event <- surv_list[[2]]
  
  # init_emp <- emp_params[[1]]
  # params_tran_array_emp <- emp_params[[2]]
  # params_tran_array_fcovar_emp <- emp_params[[3]]
  # emit_act_emp <- emp_params[[4]]
  # emit_light_emp <- emp_params[[5]]
  # corr_mat_emp <- emp_params[[6]]
  # pi_l_emp <- emp_params[[7]]
  # beta_mat_emp <- emp_params[[8]]
  
  log_sweights_vec <- numeric(dim(act)[2])
  
  
} 
###### Read in Data ###### 
if (real_data) {
  setwd("Data/")
  load("NHANES_2011_2012_2013_2014.rda")
  nhanes1 <- NHANES_mort_list[[1]] %>% filter(eligstat == 1)
  nhanes2 <- NHANES_mort_list[[2]] %>% filter(eligstat == 1)
  lmf_data <- rbind(nhanes1,nhanes2)
  
  
  load("Wavedata_G.rda")
  load("Wavedata_H.rda")
  setwd("..")
  
  # sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  # sim_num <- 1
  
  act_G <- wave_data_G[[1]]
  act_H <- wave_data_H[[1]]
  act <- rbind(act_G,act_H)
  act <- t(act[,-1])
  act0 <- act == 0
  act <- log(act)
  lod_act <- min(act[act!=-Inf],na.rm = T) - 1e-5
  act[act0] <- lod_act
  
  light_G <- wave_data_G[[2]]
  light_H <- wave_data_H[[2]]
  light <- rbind(light_G,light_H)
  light <- t(light[,-1])
  light0 <- light == 0
  light <- log(light)
  lod_light <- min(light[light!=-Inf],na.rm = T) - 1e-5
  light[light0] <- lod_light
  
  id_G <- wave_data_G[[3]]
  id_H <- wave_data_H[[3]]
  id <- rbind(id_G,id_H)
  
  seqn_com_id <- id$SEQN %in% lmf_data$seqn
  seqn_com_lmf <- lmf_data$seqn %in% id$SEQN
  
  id <- id[seqn_com_id,]
  act <- act[,seqn_com_id]
  light <- light[,seqn_com_id]
  
  lmf_data <- lmf_data[seqn_com_lmf,]
  
  if (sum(id$SEQN - lmf_data$seqn) != 0){print("LMF NOT LINKED CORRECTLY")}
  
  log_sweights_vec <- log(id$sweights/2)
  
  id <- id %>% mutate(age_disc = case_when(age <=30 ~ 1,
                                           age <=50 & age > 30 ~ 2,
                                           age <=65 & age > 50 ~ 3,
                                           age > 65 ~ 4))
  
  # id <- id %>% mutate(age_disc = case_when(age <=30 ~ 1,
  #                                          age <=40 & age > 30 ~ 2,
  #                                          age <=50 & age > 40 ~ 3,
  #                                          age <=60 & age > 50 ~ 4,
  #                                          age <=70 & age > 60 ~ 5,
  #                                          age > 70 ~ 6))
  
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  
  id$modact <- id$modact - 1
  
  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm

  if (leave_out){
    
    setwd("Data")
    load("LeaveOutMat.rda")
    setwd("..")
    leave_out_inds <- leave_out_mat[sim_num,]
    leave_out_inds <- leave_out_inds[!is.na(leave_out_inds)]
    
    act_old <- act
    light_old <- light
    id_old <- id
    surv_event_old <- surv_event
    surv_time_old <- surv_time
    log_sweights_vec_old <- log_sweights_vec
    
    first_day_vec_old <- as.numeric(id_old$PAXDAYWM)
    vcovar_mat_old <- sapply(first_day_vec_old,FirstDay2WeekInd)
    
    surv_covar_old <- list(id_old$age,
                       Vec2Mat(id_old$gender+1),
                       Vec2Mat(id_old$race+1),
                       Vec2Mat(id_old$overall_health+1),
                       Vec2Mat(id_old$education+1),
                       Vec2Mat(id_old$bmi_disc+1),
                       Vec2Mat(id_old$diabetes+1),
                       Vec2Mat(id_old$CHD+1),
                       Vec2Mat(id_old$CHF+1),
                       Vec2Mat(id_old$heart_attack+1),
                       Vec2Mat(id_old$stroke+1),
                       Vec2Mat(id_old$alcohol+1),
                       Vec2Mat(id_old$smoking+1),
                       Vec2Mat(id_old$phyfunc+1))
    
    age_vec_old <-id_old$age
    statact_vec_old <- id_old$statact
    nu_covar_mat_old <- cbind(age_vec_old/10,(age_vec_old/10)^2,statact_vec_old,statact_vec_old^2)
    
    act <- act[,-c(leave_out_inds)]
    light <- light[,-c(leave_out_inds)]
    id <- id[-c(leave_out_inds),]
    surv_event <- surv_event[-c(leave_out_inds)]
    surv_time <- surv_time[-c(leave_out_inds)]
    
    log_sweights_vec <- log(id$sweights/2)
    
  }
  
  first_day_vec <- as.numeric(id$PAXDAYWM)
  vcovar_mat <- sapply(first_day_vec,FirstDay2WeekInd)
  
  if (single_day != 0){
    single_day_mat <- sapply(first_day_vec,FirstDay2SingleDay,target_day = single_day)
    new_act <- matrix(NA,96,dim(act)[2])
    new_light <- matrix(NA,96,dim(light)[2])
    vcovar_mat <- matrix(0,96,dim(light)[2])
    
    for (i in 1:dim(act)[2]){
      new_act[,i] <- act[,i][single_day_mat[,i]==1]
      new_light[,i] <- light[,i][single_day_mat[,i]==1]
    }
    
    act <- new_act
    light <- new_light
  }
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
  
  age_vec <-id$age
  modact_vec <- id$modact
  statact_vec <- id$statact
  # nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2,modact_vec)
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)

  surv_covar <- list(id$age,
                     Vec2Mat(id$gender+1),
                     Vec2Mat(id$race+1),
                     Vec2Mat(id$overall_health+1),
                     Vec2Mat(id$education+1),
                     Vec2Mat(id$bmi_disc+1),
                     Vec2Mat(id$diabetes+1),
                     Vec2Mat(id$CHD+1),
                     Vec2Mat(id$CHF+1),
                     Vec2Mat(id$heart_attack+1),
                     Vec2Mat(id$stroke+1),
                     Vec2Mat(id$alcohol+1),
                     Vec2Mat(id$smoking+1),
                     Vec2Mat(id$phyfunc+1))
  
  surv_coef <- lapply(surv_covar[-1],SurvCovar2Coef)
  surv_coef <- append(list(.05),surv_coef)
  surv_coef_len <- unlist(lapply(surv_coef,length))
  
  combined_covar_mat <- id %>% dplyr::select(gender,race,overall_health,education,bmi_disc,diabetes,
                                             race,CHD,CHF,heart_attack,stroke,alcohol,smoking,phyfunc)
  combined_covar_mat <- lapply(combined_covar_mat, factor)
  
  if (weekend_only){
    #Only weekend data
    
    act[vcovar_mat == 0] <- NA
    light[vcovar_mat == 0] <- NA
  }
    
}

###### Initial Settings ###### 

##########
# # load("/Users/aronjr/Documents/Predoc/JointHMM/InterJMHMM2SeedNA.rda")
#

if (leave_out){
  # # Correct way to do it
  # model_name_loadin <- "JMHMM"
  # # if (!real_data){model_name <- paste0(model_name,"SimSize",sim_size)}
  # if (!tobit){model_name_loadin <- paste0(model_name_loadin,"NonTob")}
  # if (incl_surv==0){model_name_loadin <- paste0(model_name_loadin,"NoSurv")}
  # if (incl_surv==1){model_name_loadin <- paste0(model_name_loadin,"HalfSurv")}
  # if (!incl_light){model_name_loadin <- paste0(model_name_loadin,"NoLight")}
  # if (!incl_act){model_name_loadin <- paste0(model_name_loadin,"NoAct")}
  # if (misspecification){model_name_loadin <- paste0(model_name_loadin,"Misspecified",misspecification)}
  # model_name_loadin <- paste0(model_name_loadin,"Mix",mix_num,"Seed",".rda")
  # 
  # print(paste("Loading",model_name_loadin))
  # setwd("Data")
  # load(model_name_loadin)
  # setwd("..")
  
  if (incl_surv == 2){
    model_name_loadin <- paste0("JMHMMMix",mix_num,"Seed",sim_num,".rda")
    foldername <- paste0(mix_num)
  } else {
    model_name_loadin <- paste0("JMHMMNoSurvMix",mix_num,"Seed",sim_num,".rda")
    foldername <- paste0("NS",mix_num)
  }
  print(paste("Loading",model_name_loadin))
  setwd(foldername)
  load(model_name_loadin)
  setwd("..")
    
    
  init_true <- to_save[[2]][[1]]
  params_tran_array_true <- to_save[[2]][[2]]
  emit_act_true <- to_save[[2]][[3]]
  emit_light_true <- to_save[[2]][[4]]
  corr_mat_true <- to_save[[2]][[5]]
  nu_mat_true <- to_save[[2]][[6]]
  pi_l_true <- CalcPi(nu_mat_true,nu_covar_mat)
  beta_vec_true <- to_save[[2]][[7]]
  beta_age_true <- to_save[[2]][[8]]
  # re_prob <- to_save[[2]][[10]]
  re_prob <- pi_l_true
  
  if (leave_out){
    setwd("Data")
    if (incl_surv == 2){
      load(paste0("JMHMMMix",mix_num,"Seed.rda"))
    } else {
      load(paste0("JMHMMNoSurvMix",mix_num,"Seed.rda"))
    }
      
    setwd("..")
    mix_assignment_true <- apply(to_save[[2]][[10]],1,which.max)
    re_prob <- to_save[[2]][[10]][-leave_out_inds,]
  }

  #ISSUE HERE IF NON TOBIT
  #FIRST 13 should be a 12
  lambda_act_mat <- to_save[[2]][[13]]
  lambda_light_mat <- to_save[[2]][[13]]
  
}
  
#########

# keep_inds <- c(1,2,3,5,7)
# init <- init[keep_inds,]
# params_tran_array <- params_tran_array[keep_inds,,]
# emit_act <- emit_act[,,keep_inds,]
# emit_light <- emit_light[,,keep_inds,]
# corr_mat <- corr_mat[keep_inds,,]
# nu <- nu[keep_inds]
# nu2 <- nu2[keep_inds]
# beta_vec <- beta_vec[keep_inds]
# re_prob <- re_prob[,keep_inds]
# re_prob <- re_prob/rowSums(re_prob)
################## EM ##################

runif_tol <- 1
if (!randomize_init){runif_tol <- 0}

# init <- matrix(rep(.5,mix_num*2),ncol = 2)
init <- init_true

params_tran_array <- params_tran_array_true + runif(unlist(length(params_tran_array_true)),-runif_tol*2,runif_tol*2)


emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-runif_tol,runif_tol)
emit_act[,2,,] <- abs(emit_act[,2,,])

emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-runif_tol*2,runif_tol*2)
emit_light[,2,,] <- abs(emit_light[,2,,])

corr_mat <- corr_mat_true + runif(length(unlist(corr_mat_true)),-runif_tol/5,runif_tol/5)
corr_mat[corr_mat>.99] <- .99
corr_mat[corr_mat<-.99] <- -.99

if(misspecification == 1){beta_vec <- numeric(length(beta_vec))}
beta_vec <- beta_vec_true + runif(mix_num,-runif_tol,runif_tol)
beta_vec[1] <- 0

for (i in 1:length(surv_coef)){
  surv_coef[[i]] <-surv_coef[[i]]  +  runif(length(surv_coef[[i]]),-runif_tol/10,runif_tol/10)
  if (length(surv_coef[[i]]) != 1){surv_coef[[i]][1] <- 0} 
}
surv_coef[[1]] <- .065

nu_mat <- nu_mat_true
lambda_act_mat <- lambda_act_mat_true
lambda_light_mat <- lambda_light_mat_true

time_vec <- c()
pi_l <- CalcPi(nu_mat,nu_covar_mat)

#maybe comment out this
re_prob <- pi_l

surv_coef_len <- unlist(lapply(surv_coef,length))
surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

if(!incl_light){
  light <- matrix(NA,dim(light)[1],dim(light)[2])
  light_old <- matrix(NA,dim(light_old)[1],dim(light_old)[2])
}
if(!incl_act){
  act <- matrix(NA,dim(act)[1],dim(act)[2])
  act_old <- matrix(NA,dim(act_old)[1],dim(act_old)[2])
}

bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
print("Pre Alpha")
alpha <- Forward(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         lod_act = lod_act, lod_light = lod_light, 
         corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
         surv_covar = surv_covar, vcovar_mat = vcovar_mat,
         lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit,incl_surv = incl_surv)

beta <- Backward(act = act,light = light, tran_list = tran_list,
                 emit_act = emit_act,emit_light = emit_light,
                  lod_act = lod_act, lod_light =  lod_light, 
                  corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                  lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
         
print("Post Beta")
new_likelihood <- CalcLikelihood(alpha,pi_l)
# if ((incl_surv*beta_bool) == 1){new_likelihood <- new_likelihood - SurvLike(beta_vec_long)}
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood
# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
iter_count <- 1
stop_crit <- 1e-3
# if(!real_data){stop_crit <- stop_crit/10}
if(mix_num > 8){stop_crit <- stop_crit * 10}
if(mix_num > 12){stop_crit <- stop_crit * 10}
if(mix_num > 15){stop_crit <- stop_crit * 5}

while(abs(like_diff) > stop_crit){
# for (i in 1:1){
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####
  
  #### Mixing Proportion  #####
  ##### Survival ####
  re_prob <- CalcProbRE(alpha,pi_l)
  
  if(beta_bool){

    nu_mat  <- CalcNu(nu_mat,re_prob,nu_covar_mat)
    pi_l <- CalcPi(nu_mat,nu_covar_mat)
    re_prob <- CalcProbRE(alpha,pi_l)
    # 
    if (incl_surv == 2){
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,surv_covar_risk_vec,incl_surv)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],surv_coef_len)
      beta_vec <- beta_surv_coef_temp_list[[1]]
      surv_coef <- beta_surv_coef_temp_list[[2]]
      beta_se <- beta_surv_coef_se[[2]]
      
      surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]
    }
      
  }
  
  # #### Weights  #####
  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])

  
  ##### Bivariate Normal Est  #####
  
  ##### Mixture Normal Param
  if (tobit){
    
    if(incl_light){

      emit_light[1,1,,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_light[2,1,,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

      emit_light[1,2,,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_light[2,2,,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    }

    if (incl_act){
      emit_act[1,2,,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_act[2,2,,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      # #
      # #
      emit_act[1,1,,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_act[2,1,,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    }

    if (incl_act & incl_light){
      corr_mat[,1,] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      corr_mat[,2,] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

    }

  } else {
    lambda_mat_list <- CalcBinom(act,light,weights_array_list)
    lambda_act_mat <- lambda_mat_list[[1]]
    lambda_light_mat <- lambda_mat_list[[2]]
    
    emit_act[1,1,,] <- CalcMeans(act,lod_act,weights_array_wake)
    emit_act[2,1,,] <- CalcMeans(act,lod_act,weights_array_sleep)
    
    emit_act[1,2,,] <- CalcSigmas(act,lod_act,weights_array_wake,emit_act[1,1,,])
    emit_act[2,2,,] <- CalcSigmas(act,lod_act,weights_array_sleep,emit_act[2,1,,])
    
    emit_light[1,1,,] <- CalcMeans(light,lod_light,weights_array_wake)
    emit_light[2,1,,] <- CalcMeans(light,lod_light,weights_array_sleep)
    
    emit_light[1,2,,] <- CalcSigmas(light,lod_light,weights_array_wake,emit_light[1,1,,])
    emit_light[2,2,,] <- CalcSigmas(light,lod_light,weights_array_sleep,emit_light[2,1,,])
  }
  
  ###
  params_tran_array_old <- params_tran_array
  tran_gradhess_list <- CalcTranCHelper(alpha,beta,act,light,params_tran_array,
                                        emit_act,emit_light,corr_mat,
                                        pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat,
                                        lambda_act_mat, lambda_light_mat, tobit, check_tran,likelihood)
  
  params_tran_array <- LM(tran_gradhess_list[[1]],tran_gradhess_list[[2]],params_tran_array,check_tran,likelihood,pi_l)
  
  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  ####
  
  # #avoids label swapping within a group
  # #Needs work
  # for (re_ind in 1:mix_num){
  #   for (week_ind in 1:2){
  #     if (emit_act[2,1,re_ind,week_ind] > emit_act[1,1,re_ind,week_ind]){
  #       print(re_ind)
  #       print(week_ind)
  # 
  #       if (week_ind == 1){
  #         temp <- init[re_ind,1]
  #         #NO WEEKEND INIT
  #         #SWAPPING MAY SLIGHTLEY DECREASE LIKE?
  #         init[re_ind,1] <- init[re_ind,2]
  #         init[re_ind,2] <- temp
  #       }
  # 
  # 
  #       temp <- emit_act[1,,re_ind,week_ind]
  #       emit_act[1,,re_ind,week_ind] <- emit_act[2,,re_ind,week_ind]
  #       emit_act[2,,re_ind,week_ind] <- temp
  # 
  #       temp <- emit_light[1,,re_ind,week_ind]
  #       emit_light[1,,re_ind,week_ind] <- emit_light[2,,re_ind,week_ind]
  #       emit_light[2,,re_ind,week_ind] <- temp
  # 
  #       temp <- params_tran_array[re_ind,1:3,week_ind]
  #       params_tran_array[re_ind,1:3,week_ind] <- params_tran_array[re_ind,4:6,week_ind]
  #       params_tran_array[re_ind,4:6,week_ind] <- temp
  # 
  #       temp <- corr_mat[re_ind,1,week_ind]
  #       corr_mat[re_ind,1,week_ind] <- corr_mat[re_ind,2,week_ind]
  #       corr_mat[re_ind,2,week_ind] <- temp
  #     }
  #   }
  # }
  
  ##### #####
  


  # Calculate tran_list after reordering 
  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  
  alpha <- Forward(act = act,light = light,
                   init = init,tran_list = tran_list,
                   emit_act= emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light = lod_light, 
                   corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
                   event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                   lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                   surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat, tobit = tobit,incl_surv = incl_surv*beta_bool)
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  
  # if ((incl_surv*beta_bool) == 1){new_likelihood <- new_likelihood - SurvLike(beta_vec_long)}
  
  like_diff <- new_likelihood - likelihood
  
  if (like_diff < 0){
    print("Transition Likelihood Decrease")

    params_tran_array <- params_tran_array_old
    tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
    
    alpha <- Forward(act = act,light = light,
                     init = init,tran_list = tran_list,
                     emit_act = emit_act,emit_light= emit_light,
                     lod_act = lod_act, lod_light = lod_light, 
                     corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef, surv_covar_risk_vec = surv_covar_risk_vec,
                     event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                     lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                     surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                     lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat, tobit = tobit,incl_surv = incl_surv)
    
    new_likelihood <- CalcLikelihood(alpha,pi_l)
    # if ((incl_surv*beta_bool) == 1){new_likelihood <- new_likelihood - SurvLike(beta_vec_long)}
    like_diff <- new_likelihood - likelihood
  }
  
  beta <- Backward(act = act,light = light, tran_list = tran_list,
                   emit_act = emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light =  lod_light, 
                   corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
  
  print(paste("RE num:",mix_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  # if(iter_count > 5 & min(surv_event %*% re_prob) > 1e-10){beta_bool <- T}
  if (iter_count == 2){
    beta_bool <- T
    print("Starting to Est Survival and Age Mixing Effect")
  }
  
  iter_count <- iter_count + 1
  
  if (iter_count %% 10 == 0){
    
    tran_df <- ParamsArray2DF(params_tran_array)
    if (!real_data & !misspecification){
      tran_df_true <- ParamsArray2DF(params_tran_array_true,misspecification)
      tran_df_truth <- tran_df_true[,1]
      
      tran_df <- tran_df %>% mutate(truth = tran_df_truth)
      tran_df <- tran_df %>% mutate(resid = prob - truth)
    }
    
    
    true_params <- list(init_true,params_tran_array_true,
                        emit_act_true,emit_light_true,
                        corr_mat_true,nu_mat_true,
                        beta_vec_true,beta_age_true)
                        #lambda_act_mat_true,lambda_light_mat_true)
    
    est_params <- list(init,params_tran_array,
                       emit_act,emit_light,
                       corr_mat,nu_mat,
                       beta_vec,surv_coef)
                       # tran_df,
                       # re_prob,
                       # new_likelihood,
                       #lambda_act_mat,lambda_light_mat,
                       # beta_se)
    
    # bic <- CalcBIC(new_likelihood,mix_num,act,light)
    to_save <- list(true_params,est_params)#,bic)
    
    save(to_save,file = paste0("Inter",model_name))
      
  }
}

if (incl_surv != 2){
  beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef)
  beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,surv_covar_risk_vec,incl_surv)
  beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],surv_coef_len)
  beta_vec <- beta_surv_coef_temp_list[[1]]
  surv_coef <- beta_surv_coef_temp_list[[2]]
  beta_se <- beta_surv_coef_se[[2]]
  
  surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  
  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]
}

  
##### Reorder #####
#Reorder to avoid label switching
#Cluster means go from small to large by activity
reord_inds <- order(beta_vec)
# reord_inds <- c(0,rev(order(beta_vec[-1])))+1
if (!all(reord_inds == c(1:mix_num))){
  print("Swapping Labels")
  emit_act <- emit_act[,,reord_inds,]
  emit_light <- emit_light[,,reord_inds,]
  nu_mat <- nu_mat[,reord_inds]
  nu_mat <- nu_mat - nu_mat[,1]
  params_tran_array <- params_tran_array[reord_inds,,]
  corr_mat <- corr_mat[reord_inds,,]
  init <- init[reord_inds,]
  beta_vec <- beta_vec[reord_inds]
  beta_vec <- beta_vec-min(beta_vec)

  pi_l <- CalcPi(nu_mat,nu_covar_mat)
  re_prob <- CalcProbRE(alpha,pi_l)
  
  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]
  
}

for (re_ind in 1:mix_num){
  for (week_ind in 1:2){
    if (emit_act[2,1,re_ind,week_ind] > emit_act[1,1,re_ind,week_ind]){
      
      if (week_ind == 1){
        temp <- init[re_ind,1]
        #NO WEEKEND INIT
        #SWAPPING MAY SLIGHTLEY DECREASE LIKE?
        init[re_ind,1] <- init[re_ind,2]
        init[re_ind,2] <- temp
      }
      
      temp <- emit_act[1,,re_ind,week_ind]
      emit_act[1,,re_ind,week_ind] <- emit_act[2,,re_ind,week_ind]
      emit_act[2,,re_ind,week_ind] <- temp
      
      temp <- emit_light[1,,re_ind,week_ind]
      emit_light[1,,re_ind,week_ind] <- emit_light[2,,re_ind,week_ind]
      emit_light[2,,re_ind,week_ind] <- temp
      
      temp <- params_tran_array[re_ind,1:3,week_ind]
      params_tran_array[re_ind,1:3,week_ind] <- params_tran_array[re_ind,4:6,week_ind]
      params_tran_array[re_ind,4:6,week_ind] <- temp
      
      temp <- corr_mat[re_ind,1,week_ind]
      corr_mat[re_ind,1,week_ind] <- corr_mat[re_ind,2,week_ind]
      corr_mat[re_ind,2,week_ind] <- temp
    }
  }
}




if(!leave_out){
  decoded_mat <- Viterbi(act,light,vcovar_mat)
} else {
  decoded_mat <- matrix(NA,2,2)
}

tran_df <- ParamsArray2DF(params_tran_array)
if (!real_data & !misspecification){
  tran_df_true <- ParamsArray2DF(params_tran_array_true,misspecification)
  tran_df_truth <- tran_df_true[,1]
  
  tran_df <- tran_df %>% mutate(truth = tran_df_truth)
  tran_df <- tran_df %>% mutate(resid = prob - truth)
}

true_params <- list(init_true,params_tran_array_true,
                    emit_act_true,emit_light_true,
                    corr_mat_true,nu_mat_true,
                    beta_vec_true,beta_age_true,
                    lambda_act_mat_true,lambda_light_mat_true)

est_params <- list(init,params_tran_array,
                   emit_act,emit_light,
                   corr_mat,nu_mat,
                   beta_vec,surv_coef,
                   tran_df,
                   re_prob,
                   new_likelihood,
                   decoded_mat,
                   lambda_act_mat,lambda_light_mat,
                   bline_vec,cbline_vec,
                   beta_se)

# init <- to_save[[2]][[1]]
# params_tran_array <- to_save[[2]][[2]]
# emit_act <- to_save[[2]][[3]]
# emit_light <- to_save[[2]][[4]]
# corr_mat <- to_save[[2]][[5]]
# nu_mat <- to_save[[2]][[6]]
# beta_vec <- to_save[[2]][[7]]
# beta_age <- to_save[[2]][[8]]
# re_prob <- to_save[[2]][[10]]


SubsetSurvCovar <- function(surv_covar,leave_out_inds){
  for (i in 1:length(surv_covar)){
    if (is.null(dim(surv_covar[[i]]))){
      surv_covar[[i]] <- surv_covar[[i]][leave_out_inds]
    } else {
      surv_covar[[i]] <- surv_covar[[i]][leave_out_inds,]
    }
  }
  return(surv_covar)
}

if (leave_out){
  new_act <- act_old[,leave_out_inds]
  new_light <- light_old[,leave_out_inds]
  new_vcovar_mat <- vcovar_mat_old[,leave_out_inds]
  len <- dim(new_act)[1]
  num_of_people <- dim(new_act)[2]
  new_surv_covar <- SubsetSurvCovar(surv_covar_old,leave_out_inds)
  new_pi_l <- CalcPi(nu_mat,nu_covar_mat_old[leave_out_inds,])
  
  surv_covar_risk_vec_new <- SurvCovarRiskVec(new_surv_covar,surv_coef)

  surv_event_new <- surv_event_old[leave_out_inds]
  surv_time_new <- surv_time_old[leave_out_inds]
  
  log_sweights_vec_new <- log_sweights_vec_old[leave_out_inds]


  
  empty_list <- vector(mode = "list", length = 9)
  for (i in 1:9){
    # empty_list[[i]] <- vector(mode = "list", length = 9+1-i)
    empty_list[[i]] <- list()
  }
  
  empty_mat_list <- vector(mode = "list", length = 9)
  for (i in 1:9){
    # empty_list[[i]] <- vector(mode = "list", length = 9+1-i)
    empty_mat_list[[i]] <- matrix(0,mix_num,mix_num)
  }
  
  empty_vec_list <- vector(mode = "list", length = 9)
  
  
  re_prob_new_list <- empty_list
  mix_assignment_pred_list <- empty_list
  conf_mat_list <- empty_mat_list
  cindex_new_list <- empty_vec_list
  ibs_new_list <- empty_vec_list
  
  inclexcl_mat <- as.matrix(expand.grid(replicate(9, 0:1, simplify = FALSE),stringsAsFactors = FALSE))
  inclexcl_mat <- inclexcl_mat[-1,]
  
  days_incl_vec <- rowSums(inclexcl_mat)

  for (inclexcl_ind in 1:dim(inclexcl_mat)[1]){
    
    days_incl <- days_incl_vec[inclexcl_ind]
    
    inclexcl_vec <- rep(inclexcl_mat[inclexcl_ind,],each=96)
    inclexcl_ind_mat <- replicate(num_of_people,inclexcl_vec)
    
    
    new_act_working <- new_act
    new_light_working <- new_light
    # new_vcovar_mat_working <- new_vcovar_mat
    
    new_act_working[inclexcl_ind_mat == 0] <- NA
    new_light_working[inclexcl_ind_mat == 0] <- NA
    # new_vcovar_mat_working[inclexcl_ind_mat == 0] <- NA
    
    
    
    if (wake_sleep){
      
      decoded_mat <- Viterbi(new_act_working,new_light_working,new_vcovar_mat)
      alpha <- ForwardAltC(decoded_mat,init,tran_list,new_vcovar_mat)
      
    } else {
      
      alpha <- Forward(act = new_act_working,light = new_light_working,
                       init = init,tran_list = tran_list,
                       emit_act = emit_act,emit_light = emit_light,
                       lod_act = lod_act, lod_light = lod_light,
                       corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef, surv_covar_risk_vec = surv_covar_risk_vec_new,
                       event_vec = numeric(100), bline_vec = numeric(100), cbline_vec = numeric(100),
                       lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec_new,
                       surv_covar = new_surv_covar, vcovar_mat = new_vcovar_mat,
                       lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,
                       tobit = T,incl_surv = 0)
    }
    
    
    re_prob_new <- CalcProbRE(alpha,new_pi_l)
    mix_assignment_pred <- apply(re_prob_new,1,which.max)
    
    mix_assignment_pred <- c(mix_assignment_pred,c(1:mix_num))
    mix_assignment_true_ind <- c(mix_assignment_true[leave_out_inds],c(1:mix_num))
    conf_mat_ind <- table(mix_assignment_pred,mix_assignment_true_ind)
    diag(conf_mat_ind) <- diag(conf_mat_ind) - 1
    
    cindex <- CalcCindex(surv_time_new,surv_event_new,beta_vec,surv_coef,re_prob_new,new_surv_covar,surv_covar_risk_vec_new)
    ibs <- CalcIBS(surv_time_new,surv_event_new,cbline_vec,beta_vec,surv_coef,new_surv_covar,re_prob_new,incl_surv,mix_assignment_pred)
    
    
    re_prob_new_list[[days_incl]] <- append(re_prob_new_list[[days_incl]],list(re_prob_new))
    mix_assignment_pred_list[[days_incl]] <- append(mix_assignment_pred_list[[days_incl]],list(mix_assignment_pred))
    conf_mat_list[[days_incl]] <- conf_mat_list[[days_incl]] + conf_mat_ind
    cindex_new_list[[days_incl]] <- c(cindex_new_list[[days_incl]],cindex)
    ibs_new_list[[days_incl]] <- c(ibs_new_list[[days_incl]],ibs)

    
  }
  
  # leave_out_list <- list(leave_out_inds,re_prob_new_list,mix_assignment_pred_list,conf_mat_list,cindex_new_list,ibs_new_list)
  leave_out_list <- list(leave_out_inds,conf_mat_list,cindex_new_list,ibs_new_list)
      
  #Only weekend data
  new_act_weekend <- new_act
  new_light_weekend <- new_light
  
  new_act_weekend[new_vcovar_mat == 0] <- NA
  new_light_weekend[new_vcovar_mat == 0] <- NA
  
  
  alpha <- Forward(act = new_act_weekend,light = new_light_weekend,
                   init = init,tran_list = tran_list,
                   emit_act = emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light = lod_light,
                   corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec_new,
                   event_vec = numeric(100), bline_vec = numeric(100), cbline_vec = numeric(100),
                   lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec_new,
                   surv_covar = new_surv_covar, vcovar_mat = new_vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,
                   tobit = T,incl_surv = 0)
  
  re_prob_new <- CalcProbRE(alpha,new_pi_l)
  mix_assignment_pred <- apply(re_prob_new,1,which.max)
  
  conf_mat <- table(c(mix_assignment_pred,c(1:mix_num)),
                        c(mix_assignment_true[leave_out_inds],c(1:mix_num)))
  diag(conf_mat) <- diag(conf_mat) - 1
  
  
  cindex <- CalcCindex(surv_time_new,surv_event_new,beta_vec,surv_coef,re_prob_new,new_surv_covar,surv_covar_risk_vec_new)
  ibs <- CalcIBS(surv_time_new,surv_event_new,cbline_vec,beta_vec,surv_coef,new_surv_covar,re_prob_new,incl_surv,mix_assignment_pred)
   
  
  # leave_out_weekend_list <- list(leave_out_inds,re_prob_new,mix_assignment_pred,conf_mat,cindex,ibs)
  leave_out_weekend_list <- list(leave_out_inds,conf_mat,cindex,ibs)
  leave_out_to_save <- list(leave_out_list,leave_out_weekend_list)
} else {
  leave_out_to_save <- list()
}

if (real_data | save_space){
  simulated_hmm <- list()
}



mix_assignment <- apply(re_prob,1,which.max)
if (!real_data){
  tab <- table(c(mixture_mat,c(1:mix_num)),c(mix_assignment,c(1:mix_num)))
  diag(tab) <- diag(tab) - 1
} else {
  tab <- table(c(mix_assignment,c(1:mix_num)),c(mix_assignment,c(1:mix_num)))
  diag(tab) <- diag(tab) - 1
}
  
if (save_space){
  # est_params[[10]] <- 0
  est_params[[12]] <- 0
  est_params[[15]] <- 0
  est_params[[16]] <- 0
}

ibs <- CalcIBS(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment)
cindex <- CalcCindex(surv_time,surv_event,beta_vec,surv_coef,re_prob,surv_covar,surv_covar_risk_vec)
diagnostics <- list(cindex,ibs,tab)

bic <- CalcBIC(new_likelihood,mix_num,act,light)
to_save <- list(true_params,est_params,bic,leave_out_to_save,simulated_hmm,diagnostics)
# model_name <- paste0("2",model_name)
setwd("/gpfs/gsfs12/users/aronjr/JM/Routputs")
#"~/JM/Routputs"
save(to_save,file = model_name)




# surv_data <- data.frame(time = surv_time,
#                         status = surv_event,
#                         age = age_vec)
# 
# 
# fit <- coxph(Surv(time, status) ~ age, data = surv_data)
# beta_age <- fit$coefficients
# ibs <- CalcIBS(surv_time,surv_event,cbline_vec,numeric(mix_num),beta_age,age_vec,re_prob,incl_surv,mix_assignment)
# cindex <- CalcCindex(surv_time,surv_event,numeric(mix_num),beta_age,re_prob,incl_surv)
