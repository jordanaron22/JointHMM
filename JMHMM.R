################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(dplyr)

RE_type <- "norm"

real_data <- F
set_seed <- T

epsilon <- 1e-5
lepsilon <- log(epsilon)

fcovar_num <- 4
vcovar_num <- 2

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 99}
if (set_seed){set.seed(sim_num)}

mix_num <- as.numeric(commandArgs(TRUE)[1])
# sim_size <- as.numeric(commandArgs(TRUE)[2])
# RE_type <- as.character(commandArgs(TRUE)[3])
# print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",mix_num))
print(paste("Sim Seed:",sim_num,"HMM Num:",mix_num))

if(is.na(mix_num)){mix_num <- 3}
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
  mix_num <- dim(emit_act)[3]
  if (is.na(mix_num)){mix_num <- 1}
  
  lintegral_mat <- matrix(NA,nrow = mix_num, ncol = 2)
  for (i in 1:mix_num){
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

SimSurvival <- function(mixture_mat,beta_vec,lam = 1/20){
  
  num_of_people <- dim(mixture_mat)[1]
  xb <- mixture_mat %*% beta_vec
  
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
  
  return(list(time,event))
}

# rev_order <- c(dim(tran_expand_df)[2]:1
# for (i in 1:dim(tran_expand_df)[1]){
#   list_index <- paste0("[[",tran_expand_df[2,rev_order],"]]",collapse="")
# }


GenTranColVecList <- function(params_tran_array_array,mix_num,fcovar_num,vcovar_num){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  
  fcovar_mixture_vcovar_tran_list <- list()
  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (fcovar_ind in 1:fcovar_num){
    for (mixture_ind in 1:mix_num){
      for (vcovar_ind in 1:vcovar_num){
        vcovar_tran_list[[vcovar_ind]] <- Params2TranVectorT(mixture_ind,len,params_tran_array[,,fcovar_ind,vcovar_ind])
        
      }
      mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
    }
    fcovar_mixture_vcovar_tran_list[[fcovar_ind]] <- mixture_vcovar_tran_list
  }
  return(fcovar_mixture_vcovar_tran_list)
}

GenTranList <- function(params_tran_array,time_vec,mix_num,fcovar_num,vcovar_num){
  
  
  mixture_vcovar_tran_list <- list()
  fcovar_mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (fcovar_ind in 1:fcovar_num){
    for (mixture_ind in 1:mix_num){
      for (vcovar_ind in 1:vcovar_num){
        
        vcovar_tran_list[[vcovar_ind]] <- TranByTimeVec(re_ind = mixture_ind,
                                                          params_tran = params_tran_array[,,fcovar_ind,vcovar_ind],
                                                          time_vec = time_vec)
        
      }
      mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
    }
    fcovar_mixture_vcovar_tran_list[[fcovar_ind]] <- mixture_vcovar_tran_list
  }
  return(fcovar_mixture_vcovar_tran_list)
}


SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,emit_act,
                        emit_light,corr_mat,lod_act,lod_light, pi_l,beta_vec_true,missing_perc){
  mix_num <- dim(emit_act)[3]
  
  mixture_mat <- t(rmultinom(num_of_people,1,pi_l_true))
  fcovar_mat <- t(rmultinom(num_of_people,1,rep(1/fcovar_num,fcovar_num)))
  
  vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48))
  if (vcovar_num == 1){vcovar_vec <- rep(0,192)}
  vcovar_mat <- replicate(num_of_people, vcovar_vec)
  
  # tran_list <- lapply(c(1:mix_num),TranByTimeVec, params_tran = params_tran_array, time_vec = c(1:day_length))
  mixture_vec <- apply(mixture_mat,1,ChooseMixture)
  fcovar_vec <- apply(fcovar_mat,1,ChooseMixture)

  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,fcovar_num,vcovar_num)
  
  
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)
    
    mixture_ind <- mixture_vec[ind]
    fcovar_ind <- fcovar_vec[ind]
    
    tran_vcovar_list <- tran_list[[fcovar_ind]][[mixture_ind]]
    vcovar_vec_ind <- vcovar_mat[,ind]+1
    
    hidden_states <- SimulateMC(day_length,init,tran_vcovar_list,mixture_ind,vcovar_vec_ind)
    
    for (i in 1:day_length){
      
      mu_act <- emit_act[hidden_states[i] + 1,1,mixture_ind] 
      sig_act <- emit_act[hidden_states[i] + 1,2,mixture_ind]
      
      mu_light <- emit_light[hidden_states[i] + 1,1,mixture_ind] 
      sig_light <- emit_light[hidden_states[i] + 1,2,mixture_ind]
      
      bivar_corr <- corr_mat[mixture_ind,hidden_states[i] + 1]
      
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
  params_tran_array_emp <- params_tran_array
  emit_act_emp <- emit_act
  emit_light_emp <- emit_light
  corr_mat_emp <- corr_mat
  pi_l_emp <- pi_l
  
  for (re_ind in 1:mix_num){
    mixture_indicator <- (mixture_vec == re_ind)
    act_mixture <- activity_matrix[,mixture_indicator]
    light_mixture <- light_matrix[,mixture_indicator]
    mc_wake_mixture <- hidden_states_matrix[,mixture_indicator]==0
    mc_sleep_mixture <- hidden_states_matrix[,mixture_indicator]==1
    
    init_emp[re_ind,1] <- sum(mc_wake_mixture[1,])/length(mc_wake_mixture[1,])
    init_emp[re_ind,2] <- 1 - init_emp[re_ind,1] 
    
    emit_act_emp[1,1,re_ind] <- mean(act_mixture[mc_wake_mixture])
    emit_act_emp[1,2,re_ind]  <- sqrt(var(act_mixture[mc_wake_mixture]))
    emit_act_emp[2,1,re_ind]  <- mean(act_mixture[mc_sleep_mixture])
    emit_act_emp[2,2,re_ind]  <- sqrt(var(act_mixture[mc_sleep_mixture]))
    
    emit_light_emp[1,1,re_ind] <- mean(light_mixture[mc_wake_mixture])
    emit_light_emp[1,2,re_ind] <- sqrt(var(light_mixture[mc_wake_mixture]))
    emit_light_emp[2,1,re_ind] <- mean(light_mixture[mc_sleep_mixture])
    emit_light_emp[2,2,re_ind] <- sqrt(var(light_mixture[mc_sleep_mixture]))
    
    corr_mat_emp[re_ind,1] <- cor(act_mixture[mc_wake_mixture],light_mixture[mc_wake_mixture])
    corr_mat_emp[re_ind,2] <- cor(act_mixture[mc_sleep_mixture],light_mixture[mc_sleep_mixture])
    
    pi_l_emp[re_ind] <- sum(mixture_indicator)/num_of_people
  }
  
  surv_list <- SimSurvival(mixture_mat,beta_vec_true)
  
  factors <- paste0("mixture_mat[,", 2:mix_num,"]")
  cphform <- as.formula(paste("Surv(surv_list[[1]], surv_list[[2]]) ~", paste(factors, collapse="+")))
  cph_fit <- coxph(cphform)
  
  beta_vec_emp <- c(0,summary(cph_fit)$coefficients[,1])
  
  emp_params <- list(init_emp,params_tran_array_emp,emit_act_emp,
                     emit_light_emp,corr_mat_emp,pi_l_emp,beta_vec_emp)
  

  light_matrix[light_matrix<lod_light] <- lod_light
  activity_matrix[activity_matrix<lod_act] <- lod_act
  
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)
  
  activity_matrix[act_missing==1] <- NA
  light_missing[light_missing==1] <- NA
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,mixture_mat,fcovar_vec,vcovar_mat,emp_params,surv_list))
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

CalcTranHelperR <- function(init_state, new_state, act, light, tran_list_mat, emit_act, 
                           emit_light, ind_like_vec, alpha, beta, lod_act, lod_light, 
                           corr_mat, lintegral_mat, pi_l,
                           fcovar_vec, vcovar_mat){
  
  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  
  tran_vals_re_cube <- array(NA,c(len-1,num_people,num_re))
  
  for (clust_i in 1:num_re){
    
    tran_vals_re_mat <- matrix(NA, len-1,num_of_people)
    corr_vec <- corr_mat[clust_i,]
    
    for (ind in 1:num_of_people){
      tran_mat_week <- tran_list_mat[[fcovar_vec[ind]]][[clust_i]][[1]]
      tran_mat_weekend <- tran_list_mat[[fcovar_vec[ind]]][[clust_i]][[2]]
      vcovar_vec <- vcovar_mat[-1,ind]
      
      tran_mat <- tran_mat_week * (1-vcovar_vec) + tran_mat_weekend * vcovar_vec
      
      #0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
      tran_vec_ind <- (init_state-1) + ((new_state-1) * 2) + 1
      
      alpha_ind <- alpha[[ind]]
      beta_ind <- beta[[ind]]
      likelihood <- ind_like_vec[ind]
      
      
      act_ind <- act[2:len,ind]
      light_ind <- light[2:len,ind]
      
      class_vec <- logClassificationC( act_ind, light_ind, emit_act[new_state,1,clust_i], 
                                       emit_act[new_state,2,clust_i], emit_light[new_state,1,clust_i], 
                                       emit_light[new_state,2,clust_i], lod_act, lod_light, 
                                       corr_vec[new_state], lintegral_mat[clust_i,new_state] )
      
      alpha_ind_slice <- alpha_ind[1:(len-1),init_state,clust_i]
      beta_ind_slice <- beta_ind[2:(len),new_state,clust_i]
      tran_vec_slice <- tran_mat[,tran_vec_ind]
      
      tran_vals_re_ind = exp(alpha_ind_slice + beta_ind_slice + log(tran_vec_slice) + log(pi_l[clust_i]) + class_vec - likelihood)
      
      tran_vals_re_mat[,ind] = tran_vals_re_ind
      
    }
    
    tran_vals_re_cube[,,clust_i] = tran_vals_re_mat
    
  }
  
  return(tran_vals_re_cube)
}


CalcTranHelper <- function(act, light, tran_list_mat, emit_act, 
                 emit_light, ind_like_vec, alpha, beta, lod_act, lod_light, 
                 corr_mat, lintegral_mat, pi_l,fcovar_vec,vcovar_mat){
  
  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  
  tran_vals_re_array <- array(NA,c(2,2,len - 1,num_people,num_re))
  
  for(init_state in 1:2){
    for(new_state in 1:2){
      for (clust_i in 1:num_re){
        tran_vals_re_array[init_state,new_state,,,clust_i] <- CalcTranHelperC(init_state = init_state-1, new_state = new_state-1,act = act,
                                                          light = light,tran_list_mat = tran_list_mat,
                                                          emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                                          alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                                          corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                                          clust_i = clust_i-1, fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)
      }
    }
  }
  
  return(tran_vals_re_array)
}

  
CalcTranC <- function(alpha,beta,act,light,params_tran_array,emit_act,emit_light,corr_mat,pi_l,lod_act,lod_light,lintegral_mat, fcovar_vec, fcovar_mat){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  params_tran_array_working <- params_tran_array
  
  gradient <- array(0,c(2,3,mix_num,fcovar_num,vcovar_num))
  hessian_vec <- array(0,c(2,3,mix_num,fcovar_num,vcovar_num))
  cos_part_vec <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  sin_part_vec <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  cos_sin_part <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  
  # tran_list_mat <- lapply(c(1:mix_num),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_list_mat <- GenTranColVecList(params_tran_array,mix_num,fcovar_num,vcovar_num)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  tran_vals_re_array <- CalcTranHelper(act = act,
                                       light = light,tran_list_mat = tran_list_mat,
                                       emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                       alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                       corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                       fcovar_vec = fcovar_vec,vcovar_mat = vcovar_mat[-1,])
  
  
  cos_vec <- cos(2*pi*c(2:(len))/96)
  sin_vec <- sin(2*pi*c(2:(len))/96)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals <- tran_vals_re_array[init_state,new_state,,,]
      
      for (ind in 1:length(alpha)){
        
        fcovar_ind <- fcovar_vec[ind]
        vcovar_vec <- vcovar_mat[-1,ind]
        vcovar_vecR <- vcovar_vec + 1
        
        for(re_ind in 1:mix_num){
          
          if (vcovar_num == 1){
            tran_mat <- tran_list_mat[[fcovar_ind]][[re_ind]][[1]]
          } else if (vcovar_num == 2){
            tran_mat_week <- tran_list_mat[[fcovar_ind]][[re_ind]][[1]]
            tran_mat_weekend <- tran_list_mat[[fcovar_ind]][[re_ind]][[2]]
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
          
            gradient[init_state,1,re_ind,fcovar_ind,vcovar_ind] <- gradient[init_state,1,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*(vcovar_vecR==vcovar_ind))
            gradient[init_state,2,re_ind,fcovar_ind,vcovar_ind] <- gradient[init_state,2,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            gradient[init_state,3,re_ind,fcovar_ind,vcovar_ind] <- gradient[init_state,3,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            hessian_vec[init_state,1,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec[init_state,1,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,2,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec[init_state,2,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,3,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec[init_state,3,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2*(vcovar_vecR==vcovar_ind))
            
            cos_part_vec[init_state,re_ind,fcovar_ind,vcovar_ind] <- cos_part_vec[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            sin_part_vec[init_state,re_ind,fcovar_ind,vcovar_ind] <- sin_part_vec[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            cos_sin_part[init_state,re_ind,fcovar_ind,vcovar_ind] <- cos_sin_part[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec*(vcovar_vecR==vcovar_ind))
          }
          
        }
      }
    }
  }
  
  
  for(re_ind in 1:mix_num){
    for(fcovar_ind in 1:fcovar_num){
      for(vcovar_ind in 1:vcovar_num){
    
        gradient_re <- as.vector(t(gradient[,,re_ind,fcovar_ind,vcovar_ind]))
        hessian_vec_re <- hessian_vec[,,re_ind,fcovar_ind,vcovar_ind]
        hessian_re <- matrix(0,6,6)
        hessian_re1 <- matrix(0,3,3)
        hessian_re2 <- matrix(0,3,3)
        
        #### HESS 1
        
        diag(hessian_re1) <- c(hessian_vec_re[1,])
        hessian_re1[1,2] <- cos_part_vec[1,re_ind,fcovar_ind,vcovar_ind]
        hessian_re1[1,3] <- sin_part_vec[1,re_ind,fcovar_ind,vcovar_ind]
        hessian_re1[2,3] <- cos_sin_part[1,re_ind,fcovar_ind,vcovar_ind]
        
        hessian_re1 <- hessian_re1+t(hessian_re1)
        diag(hessian_re1) <- diag(hessian_re1)/2
        
        #### HESS 2
        
        diag(hessian_re2) <- c(hessian_vec_re[2,])
        hessian_re2[1,2] <- cos_part_vec[2,re_ind,fcovar_ind,vcovar_ind]
        hessian_re2[1,3] <- sin_part_vec[2,re_ind,fcovar_ind,vcovar_ind]
        hessian_re2[2,3] <- cos_sin_part[2,re_ind,fcovar_ind,vcovar_ind]
        
        hessian_re2 <- hessian_re2+t(hessian_re2)
        diag(hessian_re2) <- diag(hessian_re2)/2
        
        #######
        
        
        hessian_re[1:3,1:3] <- hessian_re1
        hessian_re[4:6,4:6] <- hessian_re2
        
        
        params_tran_array_working[re_ind,,fcovar_ind,vcovar_ind] <- params_tran_array_working[re_ind,,fcovar_ind,vcovar_ind] - solve(-hessian_re,-gradient_re)
      }
    }
    
  }
  
  
  return(params_tran_array_working)
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
  
  mu_act <- optimize(EmitLogLike, c(-20,10), act = act, light = light, 
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
  
  mu_act <- optimize(EmitLogLike, c(-20,10), act = act, light = light, 
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

CalcBIC <- function(new_likelihood,mix_num,act,light){
  
  #init tran emit pi
  num_of_param <- mix_num + mix_num*6 + mix_num*4*2 + mix_num*2 + (mix_num-1)
  
  bic <- num_of_param * log(sum(!is.na(act)) + sum(!is.na(light))) - (2 * new_likelihood)
  return(bic)
}

CalcBLHaz <- function(beta_vec, re_prob,surv_event,surv_time){
  n <- length(surv_event)
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)
  
  for (time_ind in 1:n){
    
    risk_set <- surv_time >= surv_time[time_ind]
    # risk_set <- c(time_ind:n)
    bline_vec[time_ind] <- surv_event[time_ind]/sum(re_prob[risk_set,] %*% exp(beta_vec))
  }
  
  for(time_ind in 1:n){
    anti_risk_set <- surv_time <= surv_time[time_ind]
    cbline_vec[time_ind] <- sum(bline_vec[anti_risk_set])
  }
  
  # cbline_vec <- cumsum(bline_vec)
  
  
  return(list(bline_vec,cbline_vec))
}

SurvLike <- function(beta_vec,re_prob,surv_event,bline_vec,cbline_vec){
  return(sum(log(bline_vec^surv_event) + (re_prob %*% beta_vec) * surv_event - cbline_vec*exp(re_prob %*% beta_vec)))
}

SurvLikeOpt <- function(x){
  beta_vec <- c(0,x)
  return(-SurvLike(beta_vec,re_prob,surv_event,bline_vec,cbline_vec))
}

BetaGrad <- function(x){
  beta_vec <- c(0,x)
  grad <- numeric(mix_num)
  
  for (beta_ind in 2:mix_num){
    for (ind in which(surv_event == 1)){
      
      risk_set <- surv_time >= surv_time[ind]
      
      num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind]))
      denom <- sum(re_prob[risk_set,] %*% exp(beta_vec))
      
      grad[beta_ind] <- grad[beta_ind] + surv_event[ind] * (re_prob[ind,beta_ind] - num/denom)
      
    }
  }
  return(grad[2:mix_num])
  
}

CalcBeta <- function(beta_vec,re_prob,surv_event,surv_time){
  
  mix_num <- length(beta_vec)
  l2norm <- 1
  
  while (l2norm > 1e-1){
    
    grad <- numeric(mix_num)
    hess <- numeric(mix_num)
    
    for (beta_ind in 2:mix_num){
      for (ind in which(surv_event == 1)){
        
        risk_set <- surv_time >= surv_time[ind]
        # risk_set <- c(ind:num_of_people)
        
        num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind]))
        denom <- sum(re_prob[risk_set,] %*% exp(beta_vec))
        
        grad[beta_ind] <- grad[beta_ind] + surv_event[ind] * (re_prob[ind,beta_ind] - num/denom)
        hess[beta_ind] <- hess[beta_ind] + surv_event[ind] * ((num/denom)^2 - (num/denom))
        
      }
    }
    
    
    
    
    # N-R
    gradient <- grad[2:mix_num]
    if(mix_num == 2){hessian = hess[2:mix_num]
    } else {hessian <- diag(hess[2:mix_num])}

    inf_fact <- solve(-hessian,-gradient)
    if(any(abs(inf_fact)> 2) ){inf_fact <- inf_fact/max(abs(inf_fact))}
    beta_vec_new <- beta_vec[2:mix_num] - inf_fact
    beta_vec_new <- c(0,beta_vec_new)

    # step_size <- 1
    # old_like <- SurvLike(beta_vec,re_prob,surv_event,bline_vec,cbline_vec)
    # new_like <- -Inf
    # counter <- 1
    # while (new_like < old_like){
    #   beta_vec_new <- beta_vec - grad*step_size 
    #   new_like <- SurvLike(beta_vec_new,re_prob,surv_event,bline_vec,cbline_vec)
    #   step_size <- step_size/2
    #   counter <- counter + 1
    #   if(counter %% 100 == 0){print(counter)}
    # }
    
    
    
    
    
    l2norm <- sum(sqrt((beta_vec_new - beta_vec)^2))
    
    beta_vec <- beta_vec_new
  }
  return(beta_vec)
}

readCpp( "cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )
################## EM Setup ################## 

###### True Settings ###### 
day_length <- 96 * 2
num_of_people <- 2000
missing_perc <- .1

init_true <- matrix(NA,ncol = 2,nrow = mix_num)
init_true[,1] <- seq(.1,.9,length.out = mix_num)
init_true[,2] <- 1 - init_true[,1]

params_tran_true <- matrix(rep(c(-3,-1.5,-.6,-1.9,.6,.5),mix_num),ncol = 6,byrow = T) 
params_tran_array_dim <- c(mix_num,6,fcovar_num,vcovar_num)
params_tran_array_true <- array(NA,dim = params_tran_array_dim)
for (fcovar_num_ind in 1:fcovar_num){
  for (vcovar_num_ind in 1:vcovar_num){
    params_tran_array_true[,,fcovar_num_ind,vcovar_num_ind] <- params_tran_true - fcovar_num_ind*.2
  }
}


params_tran_array_true <- params_tran_array_true + runif(unlist(length(params_tran_array_true)),-.3,.3)

emit_act_true <- array(NA, c(2,2,mix_num))
emit_act_true[2,1,] <- seq(6,7,length.out = mix_num)
emit_act_true[1,2,] <- seq(1,2,length.out = mix_num)
emit_act_true[1,1,] <- seq(1,0,length.out = mix_num)
emit_act_true[2,2,] <- seq(2,1,length.out = mix_num)

emit_light_true <- array(NA, c(2,2,mix_num))
emit_light_true[1,1,] <- seq(8,9,length.out = mix_num)
emit_light_true[1,2,] <- seq(2,3,length.out = mix_num)
emit_light_true[2,1,] <- seq(2,1,length.out = mix_num)
emit_light_true[2,2,] <- seq(1,2,length.out = mix_num)

corr_mat_true <- matrix(0,ncol = 2,nrow = mix_num)

beta_vec_true <- c(0,seq(-.75,1,length.out = mix_num-1))
# beta_vec_true <- rep(0,mix_num)

lod_act_true <- 0
lod_light_true <- 0
pi_l_true <- rep(1/mix_num,mix_num)

###### Simulate Data ###### 
if (!real_data){
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                                init=init_true,params_tran_array = params_tran_array_true,
                                emit_act = emit_act_true,emit_light = emit_light_true,corr_mat = corr_mat_true,
                                lod_act = lod_act_true,lod_light = lod_light_true,pi_l = pi_l_true,
                                missing_perc = missing_perc, beta_vec_true = beta_vec_true)
  
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  light <- simulated_hmm[[3]]
  mixture_mat <- simulated_hmm[[4]]
  fcovar_vec <- simulated_hmm[[5]]
  vcovar_mat <- simulated_hmm[[6]]
  
  emp_params <- simulated_hmm[[7]]
  surv_list <- simulated_hmm[[8]]
  
  surv_time <- surv_list[[1]]
  surv_event <- surv_list[[2]]
  
  init_emp <- emp_params[[1]]
  params_tran_array_emp <- emp_params[[2]]
  emit_act_emp <- emp_params[[3]]
  emit_light_emp <- emp_params[[4]]
  corr_mat_emp <- emp_params[[5]]
  pi_l_emp <- emp_params[[6]]
  beta_vec_emp <- emp_params[[7]]
  
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
  
  sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  
  act_G <- wave_data_G[[1]]
  act_H <- wave_data_H[[1]]
  act <- rbind(act_G,act_H)
  act <- t(act[,-1])
  act <- log(act + epsilon)
  
  light_G <- wave_data_G[[2]]
  light_H <- wave_data_H[[2]]
  light <- rbind(light_G,light_H)
  light <- t(light[,-1])
  light <- log(light + epsilon)
  
  id_G <- wave_data_G[[4]]
  id_H <- wave_data_H[[4]]
  id <- rbind(id_G,id_H)
  
  seqn_com_id <- id$SEQN %in% lmf_data$seqn
  seqn_com_lmf <- lmf_data$seqn %in% id$SEQN
  
  id <- id[seqn_com_id,]
  act <- act[,seqn_com_id]
  light <- light[,seqn_com_id]
  
  lmf_data <- lmf_data[seqn_com_lmf,]
  
  if (sum(id$SEQN - lmf_data$seqn) != 0){print("LMF NOT LINKED CORRECTLY")}
  
  log_sweights_vec <- log(id$sweights/2)
  
  id <- id %>% mutate(age_disc = case_when(age <= 10 ~ 1,
                                           age <=20 & age > 10 ~ 2,
                                           age <=35 & age > 20 ~ 3,
                                           age <=50 & age > 35 ~ 4,
                                           age <=65 & age > 50 ~ 5,
                                           age > 65 ~ 6))
  
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  
  id <- id %>% mutate(bmi_disc = case_when(BMI <= 18.5 ~ 1,
                                           BMI <=25 & BMI > 18.5 ~ 2,
                                           BMI <=30 & BMI > 25 ~ 3,
                                           BMI <=35 & BMI > 30 ~ 4,
                                           BMI <=40 & BMI > 35 ~ 5,
                                           BMI > 40 ~ 6))
  
  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
  

  #need to initialize starting values
  
  # num_of_people <- 2000
  # act <- act[,1:num_of_people]
  # light <- light[,1:num_of_people]
  # surv_event <- surv_event[1:num_of_people]
  # surv_time <- surv_time[1:num_of_people]
  
  
}

###### Initial Settings ###### 

init <- matrix(rep(.5,mix_num*2),ncol = 2)
params_tran <- matrix(rep(c(-3,-1.5,-.6,-1.9,.6,.5),mix_num),ncol = 6,byrow = T)
params_tran_array <- params_tran_array_true
# params_tran[,2:6] <- 0

runif_tol <- .1

emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-runif_tol,runif_tol)
emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-runif_tol,runif_tol)
corr_mat <- corr_mat_true + runif(length(unlist(corr_mat_true)),-runif_tol,runif_tol)
beta_vec <- beta_vec_true[2:mix_num] + runif((mix_num-1),-runif_tol,runif_tol)
beta_vec <- c(0,beta_vec)

lod_act <- lod_act_true
lod_light <- lod_light_true

log_sweights_vec <- numeric(dim(act)[2])
time_vec <- c()
pi_l <- rep(1/mix_num,mix_num)
re_prob <- matrix(rep(pi_l,num_of_people),ncol = length(pi_l),byrow = T)

# beta_vec <- c(0, seq(-1,1,length.out = (mix_num-1)))


# init <- init_emp
# params_tran <- params_tran_emp
# emit_act <- emit_act_emp
# emit_light <- emit_light_emp
# corr_mat <- corr_mat_emp
# pi_l <- pi_l_emp
# beta_vec <- beta_vec_true


################## EM ################## 
bhaz_vec <- CalcBLHaz(beta_vec,re_prob,surv_event,surv_time)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)

# tran_list <- lapply(c(1:mix_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:96))
tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,fcovar_num,vcovar_num)



alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         lod_act = lod_act, lod_light =  lod_light, 
         corr_mat = corr_mat, beta_vec = beta_vec,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
         fcovar_vec = fcovar_vec, vcovar_mat = vcovar_mat)

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
              emit_act = emit_act,emit_light = emit_light,
              lod_act = lod_act, lod_light =  lod_light, 
              corr_mat = corr_mat,lintegral_mat = lintegral_mat,
              fcovar_vec = fcovar_vec, vcovar_mat = vcovar_mat)
         

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
# break
beta_counter <- 0
beta_bool <- F

while(abs(like_diff) > 1e-4 * 1e-10){
# for (i in 1:15){
  
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####

  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  params_tran_array <- CalcTranC(alpha,beta,act,light,params_tran_array,emit_act,emit_light,corr_mat,pi_l,lod_act,lod_light,lintegral_mat,fcovar_vec,vcovar_mat)
  
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
  if (mix_num <= 1){pi_l <- c(1)}
  
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
  if(beta_bool){
    beta_vec <- CalcBeta(beta_vec,re_prob,surv_event,surv_time)
  
  
    bhaz_vec <- CalcBLHaz(beta_vec,re_prob,surv_event,surv_time)
    bline_vec <- bhaz_vec[[1]]
    cbline_vec <- bhaz_vec[[2]]
  }
  
  
  ##### Reorder #####
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  if (mix_num > 1){
    reord_inds <- c(1,order(beta_vec[2:length(beta_vec)])+1)
    emit_act <- emit_act[,,reord_inds]
    emit_light <- emit_light[,,reord_inds]
    pi_l <- pi_l[reord_inds]
    re_prob <- re_prob[,reord_inds]
    params_tran_array <- params_tran_array[reord_inds,,,]
    dim(params_tran_array) <- params_tran_array_dim
    corr_mat <- corr_mat[reord_inds,]
    init <- init[reord_inds,]
    beta_vec <- beta_vec[reord_inds]
    if(!real_data){mixture_mat <- mixture_mat[,reord_inds]}
  }
  
  for (re_ind in 1:mix_num){
    if (emit_act[2,1,re_ind] > emit_act[1,1,re_ind]){
      temp <- init[re_ind,1]
      init[re_ind,1] <- init[re_ind,2]
      init[re_ind,2] <- temp
      
      temp <- emit_act[1,,re_ind]
      emit_act[1,,re_ind] <- emit_act[2,,re_ind]
      emit_act[2,,re_ind] <- temp
      
      temp <- emit_light[1,,re_ind]
      emit_light[1,,re_ind] <- emit_light[2,,re_ind]
      emit_light[2,,re_ind] <- temp
      
      temp <- params_tran_array[re_ind,1:3,,]
      params_tran_array[re_ind,1:3,,] <- params_tran_array[re_ind,4:6,,]
      params_tran_array[re_ind,4:6,,] <- temp
      
      temp <- corr_mat[re_ind,1]
      corr_mat[re_ind,1] <- corr_mat[re_ind,2]
      corr_mat[re_ind,2] <- temp
    }
  }
  
  ##### #####

  # Calculate tran_list after reordering 
  # tran_list <- lapply(c(1:mix_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:day_length))
  # x <- lapply(c(1:mix_num),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,fcovar_num,vcovar_num)
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat, beta_vec = beta_vec,
                    event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                    lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
                    fcovar_vec = fcovar_vec, vcovar_mat = vcovar_mat)
  
  beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat,
                    fcovar_vec = fcovar_vec, vcovar_mat = vcovar_mat)
  
  
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",mix_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  beta_counter <- beta_counter + 1
  if(beta_counter > 3 & min(surv_event %*% re_prob) > 1e-10){beta_bool <- T}
  
}
true_params <- list(init_true,params_tran_array_true,emit_act_true,
                   emit_light_true,corr_mat_true,pi_l_true,
                   beta_vec_true)

if(!real_data){
  emp_params <- list(init_emp,params_tran_array_emp,emit_act_emp,
                     emit_light_emp,corr_mat_emp,pi_l_emp,
                     beta_vec_emp)
}

est_params <- list(init,params_tran_array,emit_act,
                   emit_light,corr_mat,pi_l,
                   beta_vec,bhaz_vec)

bic <- CalcBIC(new_likelihood,mix_num,act,light)

if(!real_data){to_save <- list(true_params,emp_params,est_params,bic)
} else {to_save <- list(true_params,est_params,bic)}


# if (like_diff > -1){save(to_save,file = paste0("JMHMM",mix_num,"Seed",sim_num,".rda"))}

readCpp( "cFunctions.cpp" )

