################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(dplyr)
library(numDeriv)
library(Matrix)

RE_type <- "norm"

real_data <- F
set_seed <- T

epsilon <- 1e-5
lepsilon <- log(epsilon)

fcovar_num <- 4
vcovar_num <- 2

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 999}
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

Params2TranVectorT <- function(re_ind,fcovar_ind,len,params_tran,params_tran_fcovar){
  return(t(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,re_ind=re_ind,params_tran_fcovar=params_tran_fcovar,fcovar_ind=fcovar_ind)))
}

TranByTimeVec <- function(re_ind, fcovar_ind, params_tran, params_tran_fcovar,time_vec){
  return(lapply(time_vec, Params2Tran, params_tran = params_tran,re_ind=re_ind,params_tran_fcovar=params_tran_fcovar,fcovar_ind=fcovar_ind))
}

Param2TranHelper <- function(p12,p21){
  tran <- matrix(0,2,2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1- tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  return(tran)
}

Params2Tran <- function(params_tran,time,re_ind,params_tran_fcovar,fcovar_ind){
  
  param_matrix <- matrix(params_tran[re_ind,],ncol=3,nrow=2, byrow = T)
  param_fcovar_matrix <- matrix(params_tran_fcovar[fcovar_ind,],ncol=3,nrow=2, byrow = T)
  
  param_matrix_comb <- param_matrix + param_fcovar_matrix

  tran <- Param2TranHelper(param_matrix_comb[1,1]+param_matrix_comb[1,2]*cos(2*pi*time/96)+param_matrix_comb[1,3]*sin(2*pi*time/96),
                           param_matrix_comb[2,1]+param_matrix_comb[2,2]*cos(2*pi*time/96)+param_matrix_comb[2,3]*sin(2*pi*time/96))


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

CalcLintegralMat <- function(emit_act,emit_light,emit_act_fcovar,emit_light_fcovar,corr_mat,lod_act,lod_light){
  mix_num <- dim(emit_act)[3]
  if (is.na(mix_num)){mix_num <- 1}
  
  lintegral_mat <- array(NA,dim = c(mix_num,fcovar_num,2))
  for (i in 1:mix_num){
    for (j in 1:fcovar_num){
      lintegral_mat[i,j,1] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[1,1,i]+emit_act_fcovar[1,1,j],sqrt(emit_act[1,2,i]^2+emit_act_fcovar[1,2,j]^2),
                                          emit_light[1,1,i]+emit_light_fcovar[1,1,j],sqrt(emit_light[1,2,i]^2+emit_light_fcovar[1,2,j]),
                                          corr_mat[i,1],lod_light)[[1]])
      
      lintegral_mat[i,j,2] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                          emit_act[2,1,i]+emit_act_fcovar[2,1,j],sqrt(emit_act[2,2,i]^2+emit_act_fcovar[2,2,j]^2),
                                          emit_light[2,1,i]+emit_light_fcovar[2,1,j],sqrt(emit_light[2,2,i]^2+emit_light_fcovar[2,2,j]^2),
                                          corr_mat[i,2],lod_light)[[1]])
    }
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

SimSurvival <- function(mixture_mat,beta_mat,fcovar_vec,lam = 1/20){
  
  num_of_people <- dim(mixture_mat)[1]
  
  failure_times <- numeric(num_of_people)
  censor_times <- numeric(num_of_people)
  
  for(i in 1:num_of_people){
    beta_vec <- beta_mat[,fcovar_vec[i]]
    xb <- mixture_mat[i,] %*% beta_vec
    
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


GenTranColVecList <- function(params_tran_array,params_tran_array_fcovar,mix_num,fcovar_num,vcovar_num){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  
  fcovar_mixture_vcovar_tran_list <- list()
  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (fcovar_ind in 1:fcovar_num){
    for (mixture_ind in 1:mix_num){
      for (vcovar_ind in 1:vcovar_num){
        vcovar_tran_list[[vcovar_ind]] <- Params2TranVectorT(mixture_ind,fcovar_ind,len,params_tran_array[,,vcovar_ind],params_tran_array_fcovar[,,vcovar_ind])
        
      }
      mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
    }
    fcovar_mixture_vcovar_tran_list[[fcovar_ind]] <- mixture_vcovar_tran_list
  }
  return(fcovar_mixture_vcovar_tran_list)
}

GenTranList <- function(params_tran_array,params_tran_array_fcovar,time_vec,mix_num,fcovar_num,vcovar_num){
  
  
  mixture_vcovar_tran_list <- list()
  fcovar_mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()
  
  for (fcovar_ind in 1:fcovar_num){
    for (mixture_ind in 1:mix_num){
      for (vcovar_ind in 1:vcovar_num){
        
        vcovar_tran_list[[vcovar_ind]] <- TranByTimeVec(re_ind = mixture_ind,
                                                        fcovar_ind = fcovar_ind,
                                                        params_tran = params_tran_array[,,vcovar_ind],
                                                        params_tran_fcovar = params_tran_array_fcovar[,,vcovar_ind],
                                                        time_vec = time_vec)
        
      }
      mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
    }
    fcovar_mixture_vcovar_tran_list[[fcovar_ind]] <- mixture_vcovar_tran_list
  }
  return(fcovar_mixture_vcovar_tran_list)
}


SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,params_tran_array_fcovar,
                        emit_act,emit_light,emit_act_fcovar,emit_light_fcovar,corr_mat,
                        lod_act,lod_light, pi_l,beta_mat_true,missing_perc){
  mix_num <- dim(emit_act)[3]
  
  mixture_mat <- t(rmultinom(num_of_people,1,pi_l_true))
  fcovar_mat <- t(rmultinom(num_of_people,1,rep(1/fcovar_num,fcovar_num)))
  fcovar_vec <- apply(fcovar_mat,1,ChooseMixture)
  
  vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48),rep(0,day_length-192))
  if (vcovar_num == 1){vcovar_vec <- rep(0,192)}
  vcovar_mat <- replicate(num_of_people, vcovar_vec)
  
  # tran_list <- lapply(c(1:mix_num),TranByTimeVec, params_tran = params_tran_array, time_vec = c(1:day_length))
  mixture_vec <- apply(mixture_mat,1,ChooseMixture)

  tran_list <- GenTranList(params_tran_array,params_tran_array_fcovar,c(1:day_length),mix_num,fcovar_num,vcovar_num)
  
  
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)
    
    mixture_ind <- mixture_vec[ind]
    fcovar_ind <- fcovar_vec[ind]
    
    tran_vcovar_list <- tran_list[[fcovar_ind]][[mixture_ind]]
    vcovar_vec_ind <- vcovar_mat[,ind]+1
    
    hidden_states <- SimulateMC(day_length,init,tran_vcovar_list,mixture_ind,vcovar_vec_ind)
    
    for (i in 1:day_length){
      
      mu_act <- emit_act[hidden_states[i] + 1,1,mixture_ind] + emit_act_fcovar[hidden_states[i] + 1,1,fcovar_ind]
      
      sig_act <- sqrt(emit_act[hidden_states[i] + 1,2,mixture_ind]^2 + 
                        emit_act_fcovar[hidden_states[i] + 1,2,fcovar_ind]^2)
      
      mu_light <- emit_light[hidden_states[i] + 1,1,mixture_ind] + emit_light_fcovar[hidden_states[i] + 1,1,fcovar_ind]
      
      sig_light <- sqrt(emit_light[hidden_states[i] + 1,2,mixture_ind]^2 +
                          emit_light_fcovar[hidden_states[i] + 1,2,fcovar_ind]^2)
      
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
  
  surv_list <- SimSurvival(mixture_mat,beta_mat_true,fcovar_vec)
  
  # beta_vec_long_emp <- as.vector(rev(coxph(Surv(surv_list[[1]], surv_list[[2]]) ~fcovar_mat[,c(fcovar_num:1)] + mixture_mat[,c(mix_num:1)])$coefficients))
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
  

  light_matrix[light_matrix<lod_light] <- lod_light
  activity_matrix[activity_matrix<lod_act] <- lod_act
  
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)
  
  activity_matrix[act_missing==1] <- NA
  light_missing[light_missing==1] <- NA
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,mixture_mat,fcovar_mat,vcovar_mat,surv_list))
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

CalcTranHelperRDEP <- function(init_state, new_state, act, light, tran_list_mat, emit_act, 
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
                 emit_light,emit_act_fcovar,emit_light_fcovar, ind_like_vec, alpha, beta, lod_act, lod_light, 
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
                                                          emit_act = emit_act,emit_light = emit_light,
                                                          emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
                                                          ind_like_vec = ind_like_vec,
                                                          alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                                          corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                                          clust_i = clust_i-1, fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)
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

  
CalcTranC <- function(alpha,beta,act,light,params_tran_array,emit_act,emit_light,emit_act_fcovar,emit_light_fcovar,corr_mat,pi_l,lod_act,lod_light,lintegral_mat, fcovar_vec, vcovar_mat){
  
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  params_tran_array_working <- params_tran_array
  
  gradient <- array(0,c(2,3,mix_num,vcovar_num))
  hessian_vec <- array(0,c(2,3,mix_num,vcovar_num))
  cos_part_vec <- array(0,c(2,mix_num,vcovar_num))
  sin_part_vec <- array(0,c(2,mix_num,vcovar_num))
  cos_sin_part <- array(0,c(2,mix_num,vcovar_num))
  
  gradient_fcovar <- array(0,c(2,3,fcovar_num,vcovar_num))
  hessian_vec_fcovar <- array(0,c(2,3,fcovar_num,vcovar_num))
  cos_part_vec_fcovar <- array(0,c(2,fcovar_num,vcovar_num))
  sin_part_vec_fcovar <- array(0,c(2,fcovar_num,vcovar_num))
  cos_sin_part_fcovar <- array(0,c(2,fcovar_num,vcovar_num))
  
  hessian_vec_both <- array(0,c(2,3,mix_num,fcovar_num,vcovar_num))
  cos_part_vec_both <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  sin_part_vec_both <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  cos_sin_part_both <- array(0,c(2,mix_num,fcovar_num,vcovar_num))
  
  # tran_list_mat <- lapply(c(1:mix_num),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_list_mat <- GenTranColVecList(params_tran_array,params_tran_array_fcovar,mix_num,fcovar_num,vcovar_num)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  tran_vals_re_array <- CalcTranHelper(act = act,
                                       light = light,tran_list_mat = tran_list_mat,
                                       emit_act = emit_act,emit_light = emit_light,
                                       emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
                                       ind_like_vec = ind_like_vec,
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
        # fovar_vec_ind <- fcovar_vec == fcovar_ind
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
            
            
            # Double derivatives of a's
            gradient_fcovar[init_state,1,fcovar_ind,vcovar_ind] <- gradient_fcovar[init_state,1,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*(vcovar_vecR==vcovar_ind))
            gradient_fcovar[init_state,2,fcovar_ind,vcovar_ind] <- gradient_fcovar[init_state,2,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            gradient_fcovar[init_state,3,fcovar_ind,vcovar_ind] <- gradient_fcovar[init_state,3,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            hessian_vec_fcovar[init_state,1,fcovar_ind,vcovar_ind] <- hessian_vec_fcovar[init_state,1,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*(vcovar_vecR==vcovar_ind))
            hessian_vec_fcovar[init_state,2,fcovar_ind,vcovar_ind] <- hessian_vec_fcovar[init_state,2,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2*(vcovar_vecR==vcovar_ind))
            hessian_vec_fcovar[init_state,3,fcovar_ind,vcovar_ind] <- hessian_vec_fcovar[init_state,3,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2*(vcovar_vecR==vcovar_ind))
            
            cos_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind] <- cos_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            sin_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind] <- sin_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            cos_sin_part_fcovar[init_state,fcovar_ind,vcovar_ind] <- cos_sin_part_fcovar[init_state,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec*(vcovar_vecR==vcovar_ind))
            
            # Partial derivatives of both b's and a's
            hessian_vec_both[init_state,1,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec_both[init_state,1,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*(vcovar_vecR==vcovar_ind))
            hessian_vec_both[init_state,2,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec_both[init_state,2,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2*(vcovar_vecR==vcovar_ind))
            hessian_vec_both[init_state,3,re_ind,fcovar_ind,vcovar_ind] <- hessian_vec_both[init_state,3,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2*(vcovar_vecR==vcovar_ind))
            
            cos_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind] <- cos_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            sin_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind] <- sin_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec*(vcovar_vecR==vcovar_ind))
            
            cos_sin_part_both[init_state,re_ind,fcovar_ind,vcovar_ind] <- cos_sin_part_both[init_state,re_ind,fcovar_ind,vcovar_ind] + sum(tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec*(vcovar_vecR==vcovar_ind))
            
          }
    
        }
      }
    }
  }
  
  params_tran_array_new <- params_tran_array
  params_tran_array_fcovar_new <- params_tran_array_fcovar
  
  for (init_state in 1:2){
    for(vcovar_ind in 1:vcovar_num){
      
      bdiag_mat_list <- list()
      bdiag_mat_fcovar_list <- list()
      
      for (re_ind in 1:mix_num){
        bdiag_mat <- matrix(0,3,3)
        diag(bdiag_mat) <- hessian_vec[init_state,,re_ind,vcovar_ind]
        bdiag_mat[1,2] <- cos_part_vec[init_state,re_ind,vcovar_ind]
        bdiag_mat[1,3] <- sin_part_vec[init_state,re_ind,vcovar_ind]
        bdiag_mat[2,3] <- cos_sin_part[init_state,re_ind,vcovar_ind]
        
        bdiag_mat_list <- append(bdiag_mat_list, list(bdiag_mat))
      }
      
      for (fcovar_ind in 1:fcovar_num){
        bdiag_mat_fcovar <- matrix(0,3,3)
        diag(bdiag_mat_fcovar) <- hessian_vec_fcovar[init_state,,fcovar_ind,vcovar_ind]
        bdiag_mat_fcovar[1,2] <- cos_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind]
        bdiag_mat_fcovar[1,3] <- sin_part_vec_fcovar[init_state,fcovar_ind,vcovar_ind]
        bdiag_mat_fcovar[2,3] <- cos_sin_part_fcovar[init_state,fcovar_ind,vcovar_ind]
        
        bdiag_mat_fcovar_list <- append(bdiag_mat_fcovar_list, list(bdiag_mat_fcovar))
      }
      
      hessian <- as.matrix(bdiag(list(bdiag(bdiag_mat_list),bdiag(bdiag_mat_fcovar_list))))
      
      for (re_ind in 1:mix_num){
        bdiag_mat_both <- matrix(0,3,0)
        for (fcovar_ind in 1:fcovar_num){
          bdiag_mat_both_working <- matrix(0,3,3)
          diag(bdiag_mat_both_working) <- hessian_vec_both[init_state,,re_ind,fcovar_ind,vcovar_num]
          bdiag_mat_both_working[1,2] <- cos_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind]
          bdiag_mat_both_working[1,3] <- sin_part_vec_both[init_state,re_ind,fcovar_ind,vcovar_ind]
          bdiag_mat_both_working[2,3] <- cos_sin_part_both[init_state,re_ind,fcovar_ind,vcovar_ind]
          bdiag_mat_both_working <- Symmetricize(bdiag_mat_both_working)
          bdiag_mat_both <- cbind(bdiag_mat_both,bdiag_mat_both_working)
        }
        start_ind_r <- (re_ind*3)-2
        stop_ind_r <- re_ind*3
        
        start_ind_c <- (mix_num*3)+1
        stop_ind_c <- dim(hessian)[2]
        
        hessian[c(start_ind_r:stop_ind_r),c(start_ind_c:stop_ind_c)] <- bdiag_mat_both
      }
      
      hess <- Symmetricize(hessian)
      grad <- c(as.vector(gradient[init_state,,,vcovar_ind]),as.vector(gradient_fcovar[init_state,,,vcovar_ind]))
      
      
      if (init_state == 1){
        array_inds <- c(1:3)
      } else {array_inds <- c(4:6)}
      
      params_tran_vec <- c(as.vector(t(params_tran_array[,array_inds,vcovar_ind])),
                           as.vector(t(params_tran_array_fcovar[,array_inds,vcovar_ind])))
      
      inf_mat <- solve(-hess,-grad, tol = 1e-50)
      if(any(abs(inf_mat)> 2) ){inf_mat <- inf_mat/max(abs(inf_mat))}
      
      params_tran_vec_new <- params_tran_vec - inf_mat
      
      
      params_tran_array_new[,array_inds,vcovar_ind] <- matrix(params_tran_vec_new[1:(3*mix_num)],ncol = 3,byrow = T)
      params_tran_array_fcovar_new[,array_inds,vcovar_ind] <- matrix(params_tran_vec_new[((3*mix_num)+1):length(params_tran_vec_new)],ncol = 3,byrow = T)
    }
  }
  
  return(list(params_tran_array_new,params_tran_array_fcovar_new))
}


EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,fcovar_ind,weights_vec){
  
  lintegral <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                             mu_act,sig_act,mu_light,sig_light,
                             bivar_corr,lod_light)[[1]])
  
  log_like <- logClassificationC(act,light,mu_act,sig_act,mu_light,sig_light,
                                 lod_act,lod_light,bivar_corr,lintegral)
  return(-sum(log_like * weights_vec * (fcovar_vec == fcovar_ind)))
}

CalcActMean <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind){
  
  mu_act <- optimize(EmitLogLike, c(-20,10), act = act, light = light, 
                     sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind], sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state],
                     lod_act = lod_act, lod_light = lod_light, fcovar_ind == fcovar_ind,
                     weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcActSig <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind){
  
  mu_act <- optimize(EmitLogLike, c(0,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind], sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state], fcovar_ind == fcovar_ind,
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcLightMean <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind){
  
  mu_act <- optimize(EmitLogLike, c(-20,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state], fcovar_ind == fcovar_ind,
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcLightSig <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind){
  
  mu_act <- optimize(EmitLogLike, c(0,10), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state], fcovar_ind == fcovar_ind,
                     lod_act = lod_act, lod_light = lod_light, weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

CalcBivarCorr <- function(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind){
  
  mu_act <- optimize(EmitLogLike, c(-1,1), act = act, light = light, 
                     mu_act = emit_act[mc_state,1,re_ind],sig_act = emit_act[mc_state,2,re_ind],
                     mu_light = emit_light[mc_state,1,re_ind],sig_light = emit_light[mc_state,2,re_ind],
                     lod_act = lod_act, lod_light = lod_light,  fcovar_ind == fcovar_ind,
                     weights_vec = as.vector(weights_array[,,re_ind]))$minimum
  return(mu_act)
}

UpdateNorm <- function(FUN,mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep){
  opt_param_mat <- matrix(NA,nrow = mix_num,ncol = fcovar_num)
  
  if(mc_state == 1){weights_array <- weights_array_wake}
  if(mc_state == 2){weights_array <- weights_array_sleep}
  
  for(fcovar_ind in 1:fcovar_num){
    for (re_ind in 1:dim(emit_act)[3]){
      opt_param <- FUN(mc_state,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,fcovar_ind)
      opt_param_mat[re_ind,fcovar_ind] <- opt_param
    }
  }
  
  return(opt_param_mat)
  
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

CalcBLHaz <- function(beta_mat, re_prob,surv_event,surv_time,fcovar_vec){
  n <- length(surv_event)
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)
  
  beta_vec <- beta_mat[,1]
  beta_vec_fcovar <- beta_mat[1,]
  
  for (time_ind in 1:n){
    
    risk_set <- surv_time >= surv_time[time_ind]
    # risk_set <- c(time_ind:n)
    
    denom <- sum((re_prob[risk_set,] %*% exp(beta_vec)) * (fcovar_mat[risk_set,] %*% exp(beta_vec_fcovar)))
    
    bline_vec[time_ind] <- surv_event[time_ind]/denom
  }
  
  for(time_ind in 1:n){
    anti_risk_set <- surv_time <= surv_time[time_ind]
    cbline_vec[time_ind] <- sum(bline_vec[anti_risk_set])
  }
  
  
  
  return(list(bline_vec,cbline_vec))
}

SurvLikeInt <- function(beta_mat_vec){
  beta_mat_vec <- c(0,beta_mat_vec)
  beta_mat <- matrix(beta_mat_vec,ncol = fcovar_num,nrow = mix_num)
  
  rl <- 0
  for (fcovar_ind in 1:fcovar_num){
    fcovar_index <- fcovar_vec==fcovar_ind
    beta_vec <- beta_mat[,fcovar_ind]
    rl <- rl + sum(log(bline_vec[fcovar_index]^surv_event[fcovar_index]) + (re_prob[fcovar_index,] %*% beta_vec) * surv_event[fcovar_index] - cbline_vec[fcovar_index]*exp(re_prob[fcovar_index,] %*% beta_vec))
  }
  return(-rl)
}

CalcBetaHelperInt <- function(beta_mat){
  gradient <- grad(SurvLikeInt, as.vector(beta_mat)[-1])
  invhessian <- diag(1/diag(hessian(SurvLikeInt, as.vector(beta_mat)[-1])))
  nr_fact <- as.vector(-gradient %*% -invhessian)
  nr_fact <- c(0,nr_fact)
  
  beta_mat_vec <- as.vector(beta_mat) - nr_fact
  
  beta_mat_new <- matrix(beta_mat_vec,ncol = fcovar_num,nrow = mix_num)
  return(beta_mat_new)
}

SurvLike <- function(beta_vec_long){
  beta_vec <- beta_vec_long[1:(mix_num-1)]
  beta_vec_fcovar <- beta_vec_long[mix_num:length(beta_vec_long)]
  
  beta_vec <- c(0,beta_vec)
  beta_vec_fcovar <- c(0,beta_vec_fcovar)
  
  loglike <- sum(log(bline_vec^surv_event) + 
                   ((re_prob %*% beta_vec)+(fcovar_mat %*% beta_vec_fcovar)) * surv_event - 
                   cbline_vec*exp((re_prob %*% beta_vec)+(fcovar_mat %*% beta_vec_fcovar)))

  return(-loglike)
}


CalcBetaHelper <- function(beta_vec_long){
  
  gradient <- grad(SurvLike, beta_vec_long[-c(mix_num+1,1)])
  hessian <- hessian(SurvLike, beta_vec_long[-c(mix_num+1,1)])
  
  nr_fact <- solve(-hessian,-gradient)
  nr_fact <- c(0,nr_fact)
  nr_fact <- c(nr_fact[1:mix_num],0,nr_fact[(mix_num+1):length(nr_fact)])
  
  beta_vec_long <- beta_vec_long - nr_fact
  return(beta_vec_long)
}

CalcBeta <- function(beta_mat){
  beta_vec_long <- c(beta_mat[,1],beta_mat[1,])
  l2norm <- 2
  while (l2norm > 1e-4){
    beta_vec_long_new <- CalcBetaHelper(beta_vec_long)
    l2norm <- sum(sqrt((beta_vec_long_new - beta_vec_long)^2))
    beta_vec_long <- beta_vec_long_new
  }
  
  beta_vec <- beta_vec_long[1:mix_num]
  beta_vec_fcovar <- beta_vec_long[(mix_num+1):length(beta_vec_long)]
  
  if (fcovar_num == 1){
    beta_mat <- replicate(fcovar_num, beta_vec)
  } else {
    beta_mat <- replicate(fcovar_num, beta_vec) + t(replicate(mix_num,beta_vec_fcovar))
  }
  
  return(beta_mat)
}

CalcBetaManual <- function(beta_mat,re_prob,surv_event,surv_time){
  
  mix_num <- dim(beta_mat)[1]
  l2norm <- 2
  beta_mat_new <- matrix(0,dim(beta_mat)[1],dim(beta_mat)[2])
  
  
  
  while (l2norm > 1e0){
    
    grad <- matrix(0,mix_num,fcovar_num)
    hess <- matrix(0,mix_num,fcovar_num)
    
    for (beta_ind in 1:mix_num){
      for (ind in which(surv_event == 1)){
        
        fcovar_ind <- fcovar_vec[ind]
        beta_vec <- beta_mat[,fcovar_ind]
        
        fcovar_index <- fcovar_vec==fcovar_ind
        
        if (beta_ind != 1 | fcovar_ind != 1){
          risk_set <- surv_time >= surv_time[ind]
          
          # risk_set <- c(ind:num_of_people)
          
          num <- sum(re_prob[risk_set & fcovar_index,beta_ind] * exp(beta_vec[beta_ind]))
          denom <- sum(re_prob[risk_set & fcovar_index,] %*% exp(beta_vec))
          
          grad[beta_ind,fcovar_ind] <- grad[beta_ind,fcovar_ind] + surv_event[ind] * (re_prob[ind,beta_ind] - num/denom)
          hess[beta_ind,fcovar_ind] <- hess[beta_ind,fcovar_ind] + surv_event[ind] * ((num/denom)^2 - (num/denom))
        }
      }
    }
    
    for (fcovar_ind in 1:fcovar_num){
      if (fcovar_ind == 1){ref_ind <- c(2:mix_num)
      } else {ref_ind <- c(1:mix_num)}
      
      
      gradient <- grad[ref_ind,fcovar_ind]
      if(mix_num == 2){hessian = hess[ref_ind,fcovar_ind]
      } else {hessian <- diag(hess[ref_ind,fcovar_ind])}
      
      inf_fact <- solve(-hessian,-gradient)
      if(any(abs(inf_fact)> 2) ){inf_fact <- inf_fact/max(abs(inf_fact))}
      
      beta_mat_new[ref_ind,fcovar_ind] <- beta_mat_new[ref_ind,fcovar_ind] - inf_fact

    }

    
    
    l2norm <- sum(sqrt((beta_mat_new - beta_mat)^2))
    
    beta_mat <- beta_mat_new
  }
  return(beta_mat)
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
params_tran_array_dim <- c(mix_num,6,vcovar_num)
params_tran_array_true <- array(NA,dim = params_tran_array_dim)
for (vcovar_num_ind in 1:vcovar_num){
  params_tran_array_true[,,vcovar_num_ind] <- params_tran_true 
}


params_tran_fcovar_true <- matrix(rep(c(-1,-.2,-.1,-.7,.2,.1),fcovar_num),ncol = 6,byrow = T) 
params_tran_array_fcovar_dim <- c(fcovar_num,6,vcovar_num)
params_tran_array_fcovar_true <- array(NA,dim = params_tran_array_fcovar_dim)
for (vcovar_num_ind in 1:vcovar_num){
  params_tran_array_fcovar_true[,,vcovar_num_ind] <- params_tran_fcovar_true 
}

emit_act_true <- array(NA, c(2,2,mix_num))
emit_act_true[1,1,] <- seq(4,5,length.out = mix_num)
emit_act_true[1,2,] <- seq(1,2,length.out = mix_num)
emit_act_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_act_true[2,2,] <- seq(2,1,length.out = mix_num)

emit_light_true <- array(NA, c(2,2,mix_num))
emit_light_true[1,1,] <- seq(6,7,length.out = mix_num)
emit_light_true[1,2,] <- seq(1,2,length.out = mix_num)
emit_light_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_light_true[2,2,] <- seq(1,2,length.out = mix_num)

emit_act_fcovar_true <- array(NA, c(2,2,fcovar_num))
emit_act_fcovar_true[1,1,] <- seq(2,0,length.out = fcovar_num)
emit_act_fcovar_true[1,2,] <- seq(2,0,length.out = fcovar_num)
emit_act_fcovar_true[2,1,] <- seq(0,1,length.out = fcovar_num)
emit_act_fcovar_true[2,2,] <- seq(0,1,length.out = fcovar_num)
emit_act_fcovar_true <- array(0, c(2,2,fcovar_num))

emit_light_fcovar_true <- array(NA, c(2,2,fcovar_num))
emit_light_fcovar_true[1,1,] <- seq(3,0,length.out = fcovar_num)
emit_light_fcovar_true[1,2,] <- seq(2,1,length.out = fcovar_num)
emit_light_fcovar_true[2,1,] <- seq(0,1,length.out = fcovar_num)
emit_light_fcovar_true[2,2,] <- seq(0,1,length.out = fcovar_num)
emit_light_fcovar_true <- array(0, c(2,2,fcovar_num))

corr_mat_true <- matrix(0,ncol = 2,nrow = mix_num)

beta_vec_true <- c(0,seq(-.75,1,length.out = mix_num-1))
beta_vec_fcovar_true <- seq(0,1,length.out = fcovar_num)
beta_vec_long_true <- c(beta_vec_true,beta_vec_fcovar_true)

if (fcovar_num == 1){
  beta_mat_true <- replicate(fcovar_num, beta_vec_true)
} else {
  beta_mat_true <- replicate(fcovar_num, beta_vec_true) + t(replicate(mix_num,seq(0,1,length.out = fcovar_num)))
}

# beta_vec_true <- rep(0,mix_num)

lod_act_true <- 0
lod_light_true <- 0
pi_l_true <- rep(1/mix_num,mix_num)
###### Simulate Data ###### 
if (!real_data){
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                               init=init_true,params_tran_array = params_tran_array_true, params_tran_array_fcovar = params_tran_array_fcovar_true,
                               emit_act = emit_act_true,emit_light = emit_light_true,
                               emit_act_fcovar = emit_act_fcovar_true,emit_light_fcovar = emit_light_fcovar_true,corr_mat = corr_mat_true,
                               lod_act = lod_act_true,lod_light = lod_light_true,pi_l = pi_l_true,
                               missing_perc = missing_perc, beta_mat_true = beta_mat_true)
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  light <- simulated_hmm[[3]]
  mixture_mat <- simulated_hmm[[4]]
  fcovar_mat <- simulated_hmm[[5]]
  fcovar_vec <- apply(fcovar_mat,1,ChooseMixture)
  vcovar_mat <- simulated_hmm[[6]]
  
  # emp_params <- simulated_hmm[[7]]
  surv_list <- simulated_hmm[[7]]
  
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
                                           age <=50 & age > 35 ~ 2,
                                           age <=65 & age > 50 ~ 3,
                                           age > 65 ~ 4))
  
  # id <- id %>% mutate(age_disc = case_when(age <=30 ~ 1,
  #                                          age <=40 & age > 30 ~ 2,
  #                                          age <=50 & age > 40 ~ 3,
  #                                          age <=60 & age > 50 ~ 4,
  #                                          age <=70 & age > 60 ~ 5,
  #                                          age > 70 ~ 6))
  
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  
  id <- id %>% mutate(bmi_disc = case_when(BMI <= 18.5 ~ 1,
                                           BMI <=25 & BMI > 18.5 ~ 2,
                                           BMI <=30 & BMI > 25 ~ 3,
                                           BMI <=35 & BMI > 30 ~ 4,
                                           BMI <=40 & BMI > 35 ~ 5,
                                           BMI > 40 ~ 6))
  
  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm
  
  

  #need to initialize starting values
  
  # num_of_people <- 2000
  # act <- act[,1:num_of_people]
  # light <- light[,1:num_of_people]
  # surv_event <- surv_event[1:num_of_people]
  # surv_time <- surv_time[1:num_of_people]
  
  first_day_vec <- as.numeric(id$PAXDAYWM)
  
  vcovar_mat <- sapply(first_day_vec,FirstDay2WeekInd)
  fcovar_vec <- id$age_disc
  
  to_keep_inds <- !is.na(fcovar_vec)
  
  act <- act[,to_keep_inds]
  light <- light[,to_keep_inds]
  vcovar_mat <- vcovar_mat[,to_keep_inds]
  fcovar_vec <- fcovar_vec[to_keep_inds]
  surv_event <- surv_event[to_keep_inds]
  surv_time <- surv_time[to_keep_inds]
  
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
}

###### Initial Settings ###### 

init <- matrix(rep(.5,mix_num*2),ncol = 2)
params_tran_array <- params_tran_array_true + runif(unlist(length(params_tran_array_true)),-.3,.3)
params_tran_array_fcovar <- params_tran_array_fcovar_true + runif(unlist(length(params_tran_array_fcovar_true)),-.3,.3)
# params_tran[,2:6] <- 0

runif_tol <- .2

emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-runif_tol,runif_tol)
emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-runif_tol,runif_tol)

emit_act_fcovar <- emit_act_fcovar_true + runif(length(unlist(emit_act_fcovar_true)),-runif_tol,runif_tol)
emit_light_fcovar <- emit_light_fcovar_true + runif(length(unlist(emit_light_fcovar_true)),-runif_tol,runif_tol)
emit_act_fcovar[emit_act_fcovar <0]<-0
emit_light_fcovar[emit_light_fcovar <0]<-0

# emit_act_fcovar[emit_act_fcovar < 9999] <- 0
# emit_light_fcovar[emit_light_fcovar < 9999] <- 0

corr_mat <- corr_mat_true + runif(length(unlist(corr_mat_true)),-runif_tol,runif_tol)
# beta_vec <- beta_vec_true[2:mix_num] + runif((mix_num-1),-runif_tol,runif_tol)
# beta_vec <- c(0,beta_vec)
beta_mat <- beta_mat_true + runif((mix_num*fcovar_num),-runif_tol,runif_tol)
beta_mat[1,1] <- 0

lod_act <- lod_act_true
lod_light <- lod_light_true

log_sweights_vec <- numeric(dim(act)[2])
time_vec <- c()
pi_l <- rep(1/mix_num,mix_num)
re_prob <- matrix(rep(pi_l,num_of_people),ncol = length(pi_l),byrow = T)

# beta_vec <- c(0, seq(-1,1,length.out = (mix_num-1)))


# init <- init_emp
# params_tran <- params_tran_array_emp
# emit_act <- emit_act_emp
# emit_light <- emit_light_emp
# corr_mat <- corr_mat_emp
# pi_l <- pi_l_emp
# beta_mat <- beta_mat_emp
# beta_mat[1,1] <- 0

################## EM ################## 

bhaz_vec <- CalcBLHaz(beta_mat,re_prob,surv_event,surv_time,fcovar_vec)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,emit_act_fcovar,emit_light_fcovar,corr_mat,lod_act,lod_light)

# tran_list <- lapply(c(1:mix_num),TranByTimeVec, params_tran = params_tran, time_vec = c(1:96))
tran_list <- GenTranList(params_tran_array,params_tran_array_fcovar,c(1:day_length),mix_num,fcovar_num,vcovar_num)



alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
         lod_act = lod_act, lod_light =  lod_light, 
         corr_mat = corr_mat, beta_mat = beta_mat,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
         fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
              emit_act = emit_act,emit_light = emit_light,
              emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
              lod_act = lod_act, lod_light =  lod_light, 
              corr_mat = corr_mat,lintegral_mat = lintegral_mat,
              fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)
         

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
beta_counter <- 0
beta_bool <- F
while(abs(like_diff) > 1e-4 * 1e-10){
# for (i in 1:15){
  
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####

  # init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  # params_tran_array_list <- CalcTranC(alpha,beta,act,light,params_tran_array,emit_act,emit_light,
  #                                     emit_act_fcovar,emit_light_fcovar,corr_mat,pi_l,lod_act,
  #                                     lod_light,lintegral_mat,fcovar_vec,vcovar_mat)
  params_tran_array <- params_tran_array_list[[1]]
  params_tran_array_fcovar <- params_tran_array_list[[2]]
  
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
  
  # re_prob <- CalcProbRE(alpha,pi_l)
  # pi_l <- colSums(re_prob)/dim(act)[2]
  # pi_l[pi_l<1e-100] <- 1e-100
  # if (mix_num <= 1){pi_l <- c(1)}
  
  ##### Bivariate Normal Est  #####
  
  emit_act[1,1,] <- TestMixMu()
  
  # corr_mat[,1] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # corr_mat[,2] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # 
  # emit_act[1,1,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # emit_act[2,1,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # 
  # emit_light[1,1,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # emit_light[2,1,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # 
  # emit_act[1,2,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # emit_act[2,2,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # 
  # emit_light[1,2,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  # emit_light[2,2,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep)
  
  ##### Survival ####
  # if(beta_bool){
  #   beta_mat <- CalcBeta(beta_mat)
  #   bhaz_vec <- CalcBLHaz(beta_mat,re_prob,surv_event,surv_time,fcovar_vec)
  #   bline_vec <- bhaz_vec[[1]]
  #   cbline_vec <- bhaz_vec[[2]]
  # }
  
  
  ##### Reorder #####
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  # if (mix_num > 1){
  #   reord_inds <- c(1,order(beta_mat[2:mix_num,1])+1)
  #   emit_act <- emit_act[,,reord_inds]
  #   emit_light <- emit_light[,,reord_inds]
  #   pi_l <- pi_l[reord_inds]
  #   re_prob <- re_prob[,reord_inds]
  #   params_tran_array <- params_tran_array[reord_inds,,]
  #   dim(params_tran_array) <- params_tran_array_dim
  #   corr_mat <- corr_mat[reord_inds,]
  #   init <- init[reord_inds,]
  #   beta_mat <- beta_mat[reord_inds,]
  #   if(!real_data){mixture_mat <- mixture_mat[,reord_inds]}
  # }
  
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
  tran_list <- GenTranList(params_tran_array,params_tran_array_fcovar,c(1:day_length),mix_num,fcovar_num,vcovar_num)
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,emit_act_fcovar,emit_light_fcovar,corr_mat,lod_act,lod_light)
  
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat, beta_mat = beta_mat,
                    event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                    lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
                    fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)
  
  beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat,
                    fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)
  
  
  
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




TranGrad <- function(params_tran_vec){
  len <- prod(params_tran_array_dim)
  params_tran_vec1 <- params_tran_vec[1:len]
  params_tran_vec2 <- params_tran_vec[(len+1):length(params_tran_vec)]

  params_tran_array <- array(params_tran_vec1, dim = params_tran_array_dim)
  params_tran_array_fcovar <- array(params_tran_vec2, dim = params_tran_array_fcovar_dim)

  tran_list <- GenTranList(params_tran_array,params_tran_array_fcovar,c(1:day_length),mix_num,fcovar_num,vcovar_num)
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    emit_act_fcovar = emit_act_fcovar,emit_light_fcovar = emit_light_fcovar,
                    lod_act = lod_act, lod_light =  lod_light,
                    corr_mat = corr_mat, beta_mat = beta_mat,
                    event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                    lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
                    fcovar_vec = fcovar_vec-1, vcovar_mat = vcovar_mat)

  return(CalcLikelihood(alpha,pi_l))
}

params_tran_vec <- c(as.vector(params_tran_array),as.vector(params_tran_array_fcovar))
TranGrad(params_tran_vec)
tgrad <- grad(TranGrad, params_tran_vec)

TestMixMu <- function(){
  mu_vec <- numeric(mix_num)
  for (re_ind in 1:mix_num){
    num <- 0
    denom <- sum(weights_array_wake[,,re_ind])
    for (i in 1:num_of_people){
      fcovar_ind <- fcovar_vec[i]
      num <- num + sum(weights_array_wake[,i,re_ind] * (act[,i] - emit_act_fcovar[1,1,fcovar_ind]),na.rm = T)
    }
    mu_vec[re_ind] <- num/denom
  }
  return(mu_vec)
}


