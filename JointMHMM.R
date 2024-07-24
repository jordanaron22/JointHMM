################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)

RE_type <- "norm"

real_data <- F
set_seed <- T

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 999}
if (set_seed){set.seed(sim_num)}

RE_num <- as.numeric(commandArgs(TRUE)[1])
# sim_size <- as.numeric(commandArgs(TRUE)[2])
# RE_type <- as.character(commandArgs(TRUE)[3])
# print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",RE_num))
print(paste("Sim Seed:",sim_num,"HMM Num:",RE_num))

if(is.na(RE_num)){RE_num <- 3}
# if(is.na(sim_size)){sim_size <- 0}
# if(is.na(RE_type)){RE_type <- "norm"}

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

TranByTimeVec <- function(re_ind, params_tran, time_vec){
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
  re_num <- dim(emit_act)[3]
  if (is.na(re_num)){re_num <- 1}
  
  lintegral_mat <- matrix(NA,nrow = re_num, ncol = 2)
  for (i in 1:re_num){
    lintegral_mat[i,1] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                        emit_act[1,1,i],emit_act[1,2,i],emit_light[1,1,i],emit_light[1,2,i],corr_mat[i,1],lod_light)[[1]])
    
    lintegral_mat[i,2] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                        emit_act[2,1,i],emit_act[2,2,i],emit_light[2,1,i],emit_light[2,2,i],corr_mat[i,2],lod_light)[[1]])
  }
  
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

ChooseCovar <- function(covar_vec){
  return(which(covar_vec == 1))
}

SimulateMC <- function(day_length,init,tran_list_ind,covar_ind){
  hidden_states <- numeric(day_length)
  
  for (i in 1:day_length){
    tran <- tran_list_ind[[(i-1)%%96+1]]
    
    if (i == 1) {
      hidden_states[1] <- rbinom(1,1,init[covar_ind,2])
    } else {
      hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
    }
  }
  
  
  return(hidden_states)
}

finv <- function(lam,time, randu, xb_ind){
  return((1-pexp(time,rate = lam))^(exp(xb_ind)) - randu)
}

SimSurvival <- function(covar_mat,beta_vec,lam = 1/20){
  
  num_of_people <- dim(covar_mat)[1]
  xb <- covar_mat %*% beta_vec
  
  failure_times <- numeric(num_of_people)
  censor_times <- numeric(num_of_people)
  
  for(i in 1:num_of_people){
    evnt<-runif(1)
    cens<-runif(1)
    failure_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=evnt, xb_ind=xb[i], extendInt = "yes")$root
    censor_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=cens, xb_ind=xb[i], extendInt = "yes")$root
    
  }
  
  time <- pmin(failure_times, censor_times)
  event <- as.integer(failure_times<censor_times)
  
  # time <- failure_times
  # event <- rep(1,num_of_people)
  
  return(list(time,event))
}

SimulateHMM <- function(day_length,num_of_people,init,params_tran,emit_act,
                        emit_light,corr_mat,lod_act,lod_light, pi_l,beta_vec_true,missing_perc){
  mixture_num <- dim(emit_act)[3]
  
  covar_mat <- t(rmultinom(num_of_people,1,pi_l_true))
  
  tran_list <- lapply(c(1:mixture_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:day_length))
  covar_vec <- apply(covar_mat,1,ChooseCovar)
  
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)
    
    covar_ind <- covar_vec[ind]
    
    hidden_states <- SimulateMC(day_length,init,tran_list[[covar_ind]],covar_ind)
    
    for (i in 1:day_length){
      
      mu_act <- emit_act[hidden_states[i] + 1,1,covar_ind] 
      sig_act <- emit_act[hidden_states[i] + 1,2,covar_ind]
      
      mu_light <- emit_light[hidden_states[i] + 1,1,covar_ind] 
      sig_light <- emit_light[hidden_states[i] + 1,2,covar_ind]
      
      bivar_corr <- corr_mat[covar_ind,hidden_states[i] + 1]
      
      sigma_mat <- matrix(c(sig_act^2,bivar_corr * sig_act* sig_light,
                            bivar_corr * sig_act* sig_light,sig_light^2),2,2,byrow = T)
      
      act_light <- mvrnorm(n = 1,
                           mu = c(mu_act,mu_light),
                           Sigma = sigma_mat)
    
      activity[i] <-act_light[1]
      light[i] <-act_light[2]
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
  
  init_emp <- init
  params_tran_emp <- params_tran
  emit_act_emp <- emit_act
  emit_light_emp <- emit_light
  corr_mat_emp <- corr_mat
  pi_l_emp <- pi_l
  
  for (re_ind in 1:mixture_num){
    covar_indicator <- (covar_vec == re_ind)
    act_covar <- activity_matrix[,covar_indicator]
    light_covar <- light_matrix[,covar_indicator]
    mc_wake_covar <- hidden_states_matrix[,covar_indicator]==0
    mc_sleep_covar <- hidden_states_matrix[,covar_indicator]==1
    
    init_emp[re_ind,1] <- sum(mc_wake_covar[1,])/length(mc_wake_covar[1,])
    init_emp[re_ind,2] <- 1 - init_emp[re_ind,1] 
    
    emit_act_emp[1,1,re_ind] <- mean(act_covar[mc_wake_covar])
    emit_act_emp[1,2,re_ind]  <- sqrt(var(act_covar[mc_wake_covar]))
    emit_act_emp[2,1,re_ind]  <- mean(act_covar[mc_sleep_covar])
    emit_act_emp[2,2,re_ind]  <- sqrt(var(act_covar[mc_sleep_covar]))
    
    emit_light_emp[1,1,re_ind] <- mean(light_covar[mc_wake_covar])
    emit_light_emp[1,2,re_ind] <- sqrt(var(light_covar[mc_wake_covar]))
    emit_light_emp[2,1,re_ind] <- mean(light_covar[mc_sleep_covar])
    emit_light_emp[2,2,re_ind] <- sqrt(var(light_covar[mc_sleep_covar]))
    
    corr_mat_emp[re_ind,1] <- cor(act_covar[mc_wake_covar],light_covar[mc_wake_covar])
    corr_mat_emp[re_ind,2] <- cor(act_covar[mc_sleep_covar],light_covar[mc_sleep_covar])
    
    pi_l_emp[re_ind] <- sum(covar_indicator)/num_of_people
  }
  
  surv_list <- SimSurvival(covar_mat,beta_vec_true)
  cph_fit <- coxph(Surv(surv_list[[1]], surv_list[[2]]) ~ covar_mat[,2] + covar_mat[,3] )
  beta_vec_emp <- c(0,summary(cph_fit)$coefficients[,1])
  
  emp_params <- list(init_emp,params_tran_emp,emit_act_emp,
                     emit_light_emp,corr_mat_emp,pi_l_emp,beta_vec_emp)
  

  light_matrix[light_matrix<lod_light] <- lod_light
  activity_matrix[activity_matrix<lod_act] <- lod_act
  
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)
  
  activity_matrix[act_missing==1] <- NA
  light_missing[light_missing==1] <- NA
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,covar_mat,emp_params,surv_list))
}

CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[re_ind])
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
  init_0_vec <- matrix(0,length(alpha),length(pi_l))
  init_1_vec <- matrix(0,length(alpha),length(pi_l))
  init_mat <- matrix(0,length(pi_l),2)
  
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)
  
  
  for(ind in 1:length(alpha)){ 
    ind_like <- ind_like_vec[ind]
    
    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l) - ind_like + log_sweights_vec[ind]
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l) - ind_like + log_sweights_vec[ind]

  }
  
  for (re_ind in 1:length(pi_l)){
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
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[re_ind])
    }
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
    
  }
  
  return(re_weights)
  
}

CalcTranC <- function(alpha,beta,act,light,params_tran,emit_act,emit_light,corr_mat,covar_mat_tran,pi_l,lod_act,lod_light,lintegral_mat){
  
  len <- dim(act)[1]
  re_num <- dim(emit_act)[3]
  params_tran_working <- params_tran
  
  gradient <- array(0,c(2,3,re_num))
  hessian_vec <- array(0,c(2,3,re_num))
  cos_part_vec <- matrix(0,2,re_num)
  sin_part_vec <- matrix(0,2,re_num)
  cos_sin_part <- matrix(0,2,re_num)
  
  tran_list_mat <- lapply(c(1:re_num),Params2TranVectorT, len = len, params_tran = params_tran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  tran_vals_re_00 <- CalcTranHelperC(init_state = 0,new_state = 0,act = act,
                                     light = light,tran_list_mat = tran_list_mat,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_01 <- CalcTranHelperC(init_state = 0,new_state = 1,act = act,
                                     light = light,tran_list_mat = tran_list_mat,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_10 <- CalcTranHelperC(init_state = 1,new_state = 0,act = act,
                                     light = light,tran_list_mat = tran_list_mat,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_11 <- CalcTranHelperC(init_state = 1,new_state = 1,act = act,
                                     light = light,tran_list_mat = tran_list_mat,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals <- tran_vals_re_11}
      
      for (ind in 1:length(alpha)){
        for(re_ind in 1:re_num){
          
          tran_vec <- tran_list_mat[[re_ind]]
          
          if(init_state == 1 & new_state == 1){
            # tran_prime <- -tran[1,2]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- -tran_vec[,3]
            tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
            
          } else if(init_state == 1 & new_state == 2){ 
            # tran_prime <- tran[1,1]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- tran_vec[,1]
            tran_prime_prime <- -tran_vec[,3]*tran_vec[,1]
            
          } else if(init_state == 2 & new_state == 2){ 
            # tran_prime <- -tran[2,1]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- -tran_vec[,2]
            tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
            
          } else if(init_state == 2 & new_state == 1){ 
            # tran_prime <- tran[2,2]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- tran_vec[,4]
            tran_prime_prime <- -tran_vec[,2] * tran_vec[,4]
          }
          
          cos_vec <- cos(2*pi*c(2:(len))/96)
          sin_vec <- sin(2*pi*c(2:(len))/96)
          
          gradient[init_state,1,re_ind] <- gradient[init_state,1,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime)
          gradient[init_state,2,re_ind] <- gradient[init_state,2,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*cos_vec)
          gradient[init_state,3,re_ind] <- gradient[init_state,3,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*sin_vec)
          
          hessian_vec[init_state,1,re_ind] <- hessian_vec[init_state,1,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime)
          hessian_vec[init_state,2,re_ind] <- hessian_vec[init_state,2,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2)
          hessian_vec[init_state,3,re_ind] <- hessian_vec[init_state,3,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2)
          
          cos_part_vec[init_state,re_ind] <- cos_part_vec[init_state,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec)
          sin_part_vec[init_state,re_ind] <- sin_part_vec[init_state,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec)
          
          cos_sin_part[init_state,re_ind] <- cos_sin_part[init_state,re_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec)

        }
      }
    }
  }
  
  
  for(re_ind in 1:re_num){
  
    gradient_re <- as.vector(t(gradient[,,re_ind]))
    hessian_vec_re <- hessian_vec[,,re_ind]
    hessian_re <- matrix(0,6,6)
    hessian_re1 <- matrix(0,3,3)
    hessian_re2 <- matrix(0,3,3)
    
    #### HESS 1
    
    diag(hessian_re1) <- c(hessian_vec_re[1,])
    hessian_re1[1,2] <- cos_part_vec[1,re_ind]
    hessian_re1[1,3] <- sin_part_vec[1,re_ind]
    hessian_re1[2,3] <- cos_sin_part[1,re_ind]
    
    hessian_re1 <- hessian_re1+t(hessian_re1)
    diag(hessian_re1) <- diag(hessian_re1)/2
    
    #### HESS 2
    
    diag(hessian_re2) <- c(hessian_vec_re[2,])
    hessian_re2[1,2] <- cos_part_vec[2,re_ind]
    hessian_re2[1,3] <- sin_part_vec[2,re_ind]
    hessian_re2[2,3] <- cos_sin_part[2,re_ind]
    
    hessian_re2 <- hessian_re2+t(hessian_re2)
    diag(hessian_re2) <- diag(hessian_re2)/2
    
    #######
    
    
    hessian_re[1:3,1:3] <- hessian_re1
    hessian_re[4:6,4:6] <- hessian_re2
    
  
    params_tran_working[re_ind,] <- params_tran_working[re_ind,] - solve(-hessian_re,-gradient_re)
    
  }

  
  return(params_tran_working)
}

EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,weights_vec){
  
  lintegral <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                             mu_act,sig_act,mu_light,sig_light,
                             bivar_corr,lod_light)[[1]])
  
  log_like <- logClassificationC(act,light,mu_act,sig_act,mu_light,sig_light,
                                 lod_act,lod_light,bivar_corr,lintegral)
  
  return(-sum(log_like * weights_vec))
}

CalcActMean <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind){
  
  mu_act <- optimize(EmitLogLike, c(-10,10), act = act, light = light, 
                     sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind], sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state],
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcActSig <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind){
  
  mu_act <- optimize(EmitLogLike, c(0,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind], sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state],
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcLightMean <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind){
  
  mu_act <- optimize(EmitLogLike, c(-10,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state],
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcLightSig <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind){
  
  mu_act <- optimize(EmitLogLike, c(0,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state],
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcBivarCorr <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind){
  
  mu_act <- optimize(EmitLogLike, c(-1,1), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind],sig_light = emit_light[mc_state,2,re_ind],
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

UpdateNorm <- function(FUN,mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep){
  opt_param_vec <- c()
  
  if(mc_state == 1){weights_array <- weights_array_wake}
  if(mc_state == 2){weights_array <- weights_array_sleep}
  
  for (re_ind in 1:dim(emit_act)[3]){
    opt_param <- FUN(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind)
    opt_param_vec <- c(opt_param_vec,opt_param)
  }
  
  return(opt_param_vec)
  
}

SolveCatch <- function(block_ind_hess,block_ind_grad) {
  tryCatch(
    {
      solve(block_ind_hess,block_ind_grad)
    },
    error = function(cond) {
      # message("Non-Invertible Matrix")
      numeric(length(block_ind_grad))
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
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(pi_l)))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(pi_l)))
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

CalcBIC <- function(new_likelihood,RE_num,act,light){
  
  #init tran emit pi
  num_of_param <- RE_num + RE_num*6 + RE_num*4*2 + RE_num*2 + (RE_num-1)
  
  bic <- num_of_param * log(sum(!is.na(act)) + sum(!is.na(light))) - (2 * new_likelihood)
  return(bic)
}

CalcBLHaz <- function(beta_vec, re_prob,surv_event){

  bline_vec <- numeric(num_of_people)
  for (time_ind in 1:num_of_people){
    risk_set <- c(time_ind:num_of_people)
    bline_vec[time_ind] <- surv_event[time_ind]/sum(re_prob[risk_set,] %*% exp(beta_vec))
  }
  cbline_vec <- cumsum(bline_vec)

  
  return(list(bline_vec,cbline_vec))
}

CalcBeta <- function(beta_vec,re_prob,surv_event){
  
  re_num <- length(beta_vec)
  l2norm <- 1
  
  # while (l2norm > 1e-10){
    
  grad <- numeric(re_num)
  hess <- numeric(re_num)
  
  for (beta_ind in 1:re_num){
    for (ind in 1:num_of_people){
      risk_set <- c(ind:num_of_people)
      
      num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind]))
      denom <- sum(re_prob[risk_set,] %*% exp(beta_vec))
      
      grad[beta_ind] <- grad[beta_ind] + surv_event[ind] * (re_prob[ind,beta_ind] - num/denom)
      hess[beta_ind] <- hess[beta_ind] + surv_event[ind] * ((num/denom)^2 - (num/denom))
      
    }
  }
  
  
  gradient <- grad[2:re_num]
  hessian <- diag(hess[2:re_num])
  
  beta_vec_new <- beta_vec[2:re_num] - solve(-hessian,-gradient)
  beta_vec_new <- c(0,beta_vec_new)
  
  l2norm <- sum(sqrt((beta_vec_new - beta_vec)^2))
  
  beta_vec <- beta_vec_new
  # }
  return(beta_vec)
}

readCpp( "cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )

################## EM Setup ################## 

###### True Settings ###### 
day_length <- 96 * 2
num_of_people <- 2000
missing_perc <- .1

init_true <- matrix(c(.15,.85,.3,.7,.4,.6),ncol = 2,byrow = T)
# params_tran_true <- matrix(c(-3,.9,.1,-2.2,-.75,.5,
#                              -2.5,.1,.3,-1.9,.1,.5,
#                              -2,-1,1,-3,.5,1),ncol = 6,byrow = T)
params_tran_true <- matrix(c(-3,.9,.1,-2.2,-.75,.5,
                             -2.5,.1,.3,-1.9,.1,.5,
                             -2.5,1,-1,-2.6,.5,1),ncol = 6,byrow = T)


emit_act_true <- array(NA, c(2,2,3))
emit_act1 <- matrix(c(1.5,.8,0,.5),2,2,byrow = T)
emit_act2 <- matrix(c(2.5,1.2,.25,.75),2,2,byrow = T)
emit_act3 <- matrix(c(3,1.5,.1,.6),2,2,byrow = T)
emit_act_true[,,1] <- emit_act1
emit_act_true[,,2] <- emit_act2
emit_act_true[,,3] <- emit_act3
num_re <- dim(emit_act_true)[3]

emit_light_true <- array(NA, c(2,2,3))
emit_light1 <- matrix(c(7,2.5,2,1.5),2,2,byrow = T)
emit_light2 <- matrix(c(8,3,3,2),2,2,byrow = T)
emit_light3 <- matrix(c(7.5,2,1,2),2,2,byrow = T)
emit_light_true[,,1] <- emit_light1
emit_light_true[,,2] <- emit_light2
emit_light_true[,,3] <- emit_light3

corr_mat_true <- matrix(c(.6,.8,
                          -.5,-.7,
                          -.4,.55),ncol = 2,byrow = T)

beta_vec_true <- c(0,-.5,.75)

lod_act_true <- 0
lod_light_true <- 0
pi_l_true <- c(.25,.4,.35)
simulated_hmm <- SimulateHMM(day_length,num_of_people,
                              init=init_true,params_tran = params_tran_true,
                              emit_act = emit_act_true,emit_light = emit_light_true,corr_mat = corr_mat_true,
                              lod_act = lod_act_true,lod_light = lod_light_true,pi_l = pi_l_true,
                              missing_perc = missing_perc, beta_vec_true = beta_vec_true)

mc <- simulated_hmm[[1]]
act <- simulated_hmm[[2]]
light <- simulated_hmm[[3]]
covar_mat <- simulated_hmm[[4]]

emp_params <- simulated_hmm[[5]]
surv_list <- simulated_hmm[[6]]

surv_time <- surv_list[[1]]
surv_event <- surv_list[[2]]

init_emp <- emp_params[[1]]
params_tran_emp <- emp_params[[2]]
emit_act_emp <- emp_params[[3]]
emit_light_emp <- emp_params[[4]]
corr_mat_emp <- emp_params[[5]]
pi_l_emp <- emp_params[[6]]
beta_vec_emp <- emp_params[[7]]

surv_order <- order(surv_time)
surv_time <- surv_time[surv_order]
surv_event <- surv_event[surv_order]

mc <- mc[,surv_order]
act <- act[,surv_order]
light <- light[,surv_order]
covar_mat <- covar_mat[surv_order,]


###### Initial Settings ###### 

init <- matrix(rep(.5,RE_num*2),ncol = 2)
params_tran <- matrix(rep(c(-3,.9,.1,-2.2,-.75,.5),RE_num),ncol = 6,byrow = T)
# params_tran[,2:6] <- 0

runif_tol <- .2

emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-runif_tol,runif_tol)
emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-runif_tol,runif_tol)
corr_mat <- corr_mat_true + runif(length(unlist(corr_mat_true)),-runif_tol,runif_tol)
beta_vec <- beta_vec_true[2:RE_num] + runif((RE_num-1),-runif_tol,runif_tol)
beta_vec <- c(0,beta_vec)

lod_act <- lod_act_true
lod_light <- lod_light_true

log_sweights_vec <- numeric(dim(act)[2])
time_vec <- c()
pi_l <- rep(1/RE_num,RE_num)
re_prob <- matrix(rep(pi_l,num_of_people),ncol = length(pi_l),byrow = T)

# beta_vec <- c(0, seq(-1,1,length.out = (RE_num-1)))


# init <- init_emp
# params_tran <- params_tran_emp
# emit_act <- emit_act_emp
# emit_light <- emit_light_emp
# corr_mat <- corr_mat_emp
# pi_l <- pi_l_emp
# beta_vec <- beta_vec_true


################## EM ################## 

bhaz_vec <- CalcBLHaz(beta_vec,re_prob,surv_event)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)

tran_list <- lapply(c(1:RE_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:96))

alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         lod_act = lod_act, lod_light =  lod_light, 
         corr_mat = corr_mat, beta_vec = beta_vec,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]))

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
              emit_act = emit_act,emit_light = emit_light,
              lod_act = lod_act, lod_light =  lod_light, 
              corr_mat = corr_mat,lintegral_mat = lintegral_mat)
         

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
# break
while(abs(like_diff) > 1e-4*1e-10){
# for (i in 1:2){
  
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####

  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  params_tran <- CalcTranC(alpha,beta,act,light,params_tran,emit_act,emit_light,corr_mat,covar_mat_tran,pi_l,lod_act,lod_light,lintegral_mat)
  
  ##### Weights  #####
  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])
  
  # # Uncomment to use real mc states as weights
  # weights_array_wake <- 1-mc
  # weights_array_sleep <- mc
  # dim(weights_array_wake)[3] <- 1
  # dim(weights_array_sleep)[3] <- 1
  
  ##### Mixing Proportion  #####
  
  re_prob <- CalcProbRE(alpha,pi_l)
  pi_l <- colSums(re_prob)/dim(act)[2]
  pi_l[pi_l<1e-100] <- 1e-100
  if (RE_num <= 1){pi_l <- c(1)}
  
  ##### Bivariate Normal Est  #####
  
  corr_mat[,1] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  corr_mat[,2] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)

  emit_act[1,1,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  emit_act[2,1,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)

  emit_light[1,1,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  emit_light[2,1,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)

  emit_act[1,2,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  emit_act[2,2,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)

  emit_light[1,2,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  emit_light[2,2,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  
  ##### Survival ####
  
  beta_vec <- CalcBeta(beta_vec,re_prob,surv_event)

  bhaz_vec <- CalcBLHaz(beta_vec,re_prob,surv_event)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]
  
  
  ##### Reorder #####
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  if (RE_num > 1){
    reord_inds <- c(1,order(beta_vec[2:length(beta_vec)])+1)
    emit_act <- emit_act[,,reord_inds]
    emit_light <- emit_light[,,reord_inds]
    pi_l <- pi_l[reord_inds]
    re_prob <- re_prob[,reord_inds]
    params_tran <- params_tran[reord_inds,]
    corr_mat <- corr_mat[reord_inds,]
    init <- init[reord_inds,]
    beta_vec <- beta_vec[reord_inds]
    covar_mat <- covar_mat[,reord_inds]
  }
  
  ##### #####

  # Calculate tran_list after reordering 
  tran_list <- lapply(c(1:RE_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:day_length))
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat, beta_vec = beta_vec,
                    event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                    lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]))
  
  beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat)
  
  
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",RE_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
}
true_params <- list(init_true,params_tran_true,emit_act_true,
                   emit_light_true,corr_mat_true,pi_l_true)

emp_params <- list(init_emp,params_tran_emp,emit_act_emp,
                   emit_light_emp,corr_mat_emp,pi_l_emp)

est_params <- list(init,params_tran,emit_act,
                   emit_light,corr_mat,pi_l)

bic <- CalcBIC(new_likelihood,RE_num,act,light)
to_save <- list(true_params,emp_params,est_params,bic)

# save(to_save,file = paste0("MHMM",RE_num,"Seed",sim_num,".rda"))
