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

RE_type <- "norm"

beta_bool <- F
real_data <- T
set_seed <- T

epsilon <- 1e-3
lepsilon <- log(epsilon)

vcovar_num <- 2
fcovar_num <- 4

NR_count <- 2
lm_step_size <- 100

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 10}
if (set_seed){set.seed(sim_num)}

mix_num <- as.numeric(commandArgs(TRUE)[1])
# sim_size <- as.numeric(commandArgs(TRUE)[2])
# RE_type <- as.character(commandArgs(TRUE)[3])
# print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",mix_num))
print(paste("Sim Seed:",sim_num,"HMM Num:",mix_num))

if(is.na(mix_num)){mix_num <- 5}
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

SimSurvival <- function(mixture_vec,beta_vec,beta_age,age_vec,lam = 1/20){
  
  num_of_people <- length(mixture_vec)
  
  failure_times <- numeric(num_of_people)
  censor_times <- numeric(num_of_people)
  
  for(i in 1:num_of_people){
    xb <- beta_vec[mixture_vec[i]] + beta_age*age_vec[i]
    
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
      vcovar_tran_list[[vcovar_ind]] <- Params2TranVectorT(mixture_ind,len,params_tran_array[,,vcovar_ind])
      
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
      
      vcovar_tran_list[[vcovar_ind]] <- TranByTimeVec(re_ind = mixture_ind,
                                                      params_tran = params_tran_array[,,vcovar_ind],
                                                      time_vec = time_vec)
      
    }
    mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
  }
  
  return(mixture_vcovar_tran_list)
}


SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,
                        emit_act,emit_light,corr_mat,
                        lod_act,lod_light, nu,nu2,beta_vec_true,beta_age_true,missing_perc){
  
  mix_num <- dim(emit_act)[3]
  
  age_vec <- floor(runif(num_of_people,0,81))
  pi_l_true <- CalcPi(nu_true,nu2_true,age_vec)
  mixture_vec <- rMultinom(pi_l_true,1)
  
  vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48))
  vcovar_vec <- rep(vcovar_vec,day_length/192)
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
  
  surv_list <- SimSurvival(mixture_vec,beta_vec_true,beta_age_true,age_vec)
  cox_mod <- coxph(Surv(surv_list[[1]], surv_list[[2]]) ~age_vec + as.factor(mixture_vec))
  beta_vec_long_emp <- as.vector(as.vector(cox_mod$coefficients))
  
  
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
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,mixture_vec,age_vec,vcovar_mat,surv_list))
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
                           corr_mat, lintegral_mat, pi_l,vcovar_mat){
  
  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  
  tran_vals_re_array <- array(NA,c(2,2,len - 1,num_people,num_re))
  
  for(init_state in 1:2){
    for(new_state in 1:2){
      for (clust_i in 1:num_re){
        tran_vals_re_array[init_state,new_state,,,clust_i] <- CalcTranHelperC(init_state = init_state-1, new_state = new_state-1,act = act,
                                                          light = light,tran_list_mat = tran_list_mat,
                                                          emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
                                                          emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
                                                          ind_like_vec = ind_like_vec,
                                                          alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                                          corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                                          clust_i = clust_i-1, vcovar_mat = vcovar_mat)
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

  
CalcTranC <- function(alpha,beta,act,light,params_tran_array,emit_act,emit_light,corr_mat,pi_l,lod_act,lod_light,lintegral_mat, vcovar_mat){
  
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
                                       vcovar_mat = vcovar_mat[-1,])
  
  
  cos_vec <- cos(2*pi*c(2:(len))/96)
  sin_vec <- sin(2*pi*c(2:(len))/96)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals <- tran_vals_re_array[init_state,new_state,,,]
      
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
  
  grad_desc <- TRUE
  if(iter_count > NR_count){grad_desc <- FALSE}
  
  grad_hess_list <- CalcGradHess(gradient,hessian_vec,cos_part_vec,sin_part_vec,cos_sin_part)
  grad_array <- grad_hess_list[[1]]
  hess_array <- grad_hess_list[[2]]
  
  return(LM(grad_array,hess_array,params_tran_array))

}

LM <- function(grad_array,hess_array,params_tran_array){
  params_tran_array_new <- params_tran_array
  new_likelihood <- -Inf
  # while (new_likelihood - likelihood < -1e-5 ){
  
  for (re_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      inf_fact <- solve(hess_array[,,re_ind,vcovar_ind],grad_array[,re_ind,vcovar_ind])
      
      step_size <- 1
      while(max(abs(inf_fact)) > 3){
        step_fact <- matrix(0,6,6)
        diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) * step_size
        
        inf_fact <- solve(hess_array[,,re_ind,vcovar_ind]+step_fact,grad_array[,re_ind,vcovar_ind])
        step_size <- step_size * 10
      }
      params_tran_array_new[re_ind,,vcovar_ind] <- params_tran_array_new[re_ind,,vcovar_ind] - inf_fact
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
  
  log_like <- logClassificationC(act_vec,light_vec,
                                 mu_act,
                                 sig_act,
                                 mu_light,
                                 sig_light,
                                 lod_act,lod_light,bivar_corr,lintegral)
  
  log_like[log_like == -Inf] <- -9999

  return(-sum(log_like * weights_mat[vcovar_vec_indicator]))
}

CalcActMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat){
  
  mu_act <- optimize(EmitLogLike, c(log(.001),10), act = act, light = light, 
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
  
  mu_act <- optimize(EmitLogLike, c(-100,10), act = act, light = light, 
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
  
  mu_act <- optimize(EmitLogLike, c(0.01,30), act = act, light = light, 
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
  
  mu_act <- optimize(EmitLogLike, c(-.98,.98), act = act, light = light, 
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

CalcBLHaz <- function(beta_age,vec, re_prob,surv_event,surv_time,age_vec){
  n <- length(surv_event)
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)
  
  for (time_ind in 1:n){
    risk_set <- surv_time >= surv_time[time_ind]
    denom <- sum((re_prob[risk_set,] %*% exp(beta_vec)) * exp(beta_age * age_vec[risk_set]))
    bline_vec[time_ind] <- surv_event[time_ind]/denom
  }
  
  for(time_ind in 1:n){
    anti_risk_set <- surv_time <= surv_time[time_ind]
    cbline_vec[time_ind] <- sum(bline_vec[anti_risk_set])
  }
  
  return(list(bline_vec,cbline_vec))
}

SurvLike <- function(beta_vec_long){
  
  beta_age <- beta_vec_long[1]
  beta_vec <- beta_vec_long[-1]
  
  beta_vec <- c(0,beta_vec)
  
  
  loglike <- sum(log(bline_vec^surv_event) + 
                   ((re_prob %*% beta_vec)+(beta_age * age_vec)) * surv_event - 
                   cbline_vec*exp((re_prob %*% beta_vec)+(beta_age * age_vec)))
  
  return(-loglike)
}

CalcBetaHelper <- function(beta_vec_long){
  
  gradient <- grad(SurvLike, beta_vec_long)
  hessian <- hessian(SurvLike, beta_vec_long)
  
  nr_fact <- solve(hessian,gradient)
  step_size <- 1
  step_mat <- matrix(0,dim(hessian)[1],dim(hessian)[2])
  while(max(abs(nr_fact)) > .5 | all(nr_fact == -1)){
    diag(step_mat) <- (diag(hessian) * (step_size-1)) + 1
    nr_fact <- SolveCatch(hessian + step_mat,gradient)
    step_size <- step_size * 10
  }
  
  beta_vec_long <- beta_vec_long - nr_fact
  return(beta_vec_long)
}

CalcBeta <- function(beta_vec_long){
  l2norm <- 2
  beta_it <- 1
  while (l2norm > 1e-1 & beta_it < 50){
    beta_vec_long_new <- CalcBetaHelper(beta_vec_long)
    l2norm <- sum(sqrt((beta_vec_long_new - beta_vec_long)^2))
    beta_vec_long <- beta_vec_long_new
    beta_it <- beta_it+1
  }
  
  
  
  
  return(beta_vec_long)
}

Vec2Mat <- function(vect){
  mat <- matrix(0,nrow = length(vect),ncol = max(vect))
  
  for(i in 1:length(vect)){
    mat[i,vect[i]] <- 1
  }
  return (mat)
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

ParamsArray2DF <- function(params_tran_array){
  
  tran_df <- data.frame(prob = c(),
                        type = c(),
                        time = c(),
                        age = c(),
                        weekend = c(),
                        mixture = c())
  for (re_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      tran_mat <- Params2TranVectorTresid(re_ind,96,params_tran_array[,,vcovar_ind])
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

CalcPiHelper <- function(nu, nu2,age){
  pi_ind <- exp(nu*age + nu2*age^2)/sum(exp(nu*age + nu2*age^2))
  if (any(is.na(pi_ind))){
    pi_ind[is.na(pi_ind)] <- 1
    pi_ind <- pi_ind/sum(pi_ind)
  }
  return(pi_ind)
}

CalcPi <- function(nu,nu2,age_vec){
  return(t(sapply(age_vec/10,CalcPiHelper,nu=nu, nu2=nu2)))
}

CalcNu <- function(nu,nu2,re_prob,age_vec){
  age_vec <- age_vec/10
  gradient_nu <- numeric(mix_num)
  hess_nu <- matrix(0,mix_num*2,mix_num*2)
  for (ind in 1:num_of_people){
    age_ind <- age_vec[ind]
    age_ind_vec <- c(age_ind,age_ind^2)
    age_ind_mat <- age_ind_vec %*% t(age_ind_vec)
    
    num <- exp(nu*age_ind + nu2*age_ind^2)
    denom <- sum(num)
    p_vec <- num/denom
    
    gradient_nu <- gradient_nu + (re_prob[ind,] - p_vec) %x% age_ind_vec
    
    p_vec[1]*(1-p_vec[1]) * age_ind_mat
    
    p_mat <- p_vec %x% t(p_vec)
    diag(p_mat) <- -p_vec * (1-p_vec)
    hess_nu <- hess_nu + p_mat %x% age_ind_mat
  }
  
  hess_nu <- hess_nu[c(-1,-2),c(-1,-2)]
  gradient_nu <- gradient_nu[c(-1,-2)]
  
  inf_mat <- SolveCatch(hess_nu,gradient_nu)
  step_size <- .01
  while(max(abs(inf_mat)) > 2 | all(inf_mat == -1)){
    
    step_fact <- matrix(0,dim(hess_nu)[1],dim(hess_nu)[2])
    
    diag(step_fact) <- step_size
    
    # diag(step_fact) <- diag(hess_nu) * step_size
    # if (any(diag(step_fact) == 0)){diag(step_fact) <- step_size}
    
    
    inf_mat <- SolveCatch(hess_nu+step_fact,gradient_nu)
    step_size <- step_size * 10
  }
  
  
  nu_long <- as.vector(rbind(nu[-1],nu2[-1])) - inf_mat
  nu_mat <-matrix(nu_long,nrow = 2)
  nu <- c(0,nu_mat[1,])
  nu2 <- c(0,nu_mat[2,])
  
  return(list(nu,nu2))
}

NuLogLike <- function(nu_long){
  nu_mat <- matrix(nu_long,nrow = 2, byrow = T)
  nu_mat <- cbind(c(0,0),nu_mat)
  nu <- nu_mat[1,]
  nu2 <- nu_mat[2,]
  age_vec <- age_vec/10
  loglike <- 0
  for (ind in 1:num_of_people){
    age_ind_vec <- c(age_vec[ind],age_vec[ind]^2)
    for (re_ind in 1:mix_num){
      cont <- (age_ind_vec %*% nu_mat[,re_ind]) - log(sum(exp(age_ind_vec %*% nu_mat)))
      loglike_ind <- re_prob[ind,re_ind]*cont
      loglike <- loglike + loglike_ind
    }
  }
  return(loglike)
}

CalcNuNum <- function(nu,nu2){
  nu_long <- c(nu[-1],nu2[-1])
  grad_nu_num <- grad(NuLogLike,nu_long)  
  hess_nu_num <- hessian(NuLogLike,nu_long) 
  
  inf_mat <- solve(hess_nu_num,grad_nu_num)
  
  step_size <- 1.1
  if (max(abs(inf_mat)) > .25){
    diag(hess_nu_num) <- diag(hess_nu_num) * step_size
    inf_mat <- solve(hess_nu_num,grad_nu_num)
    step_size <- step_size * 10
  }
  
  
  nu_long_new <- nu_long - inf_mat
  nu_mat <- matrix(nu_long_new,nrow = 2, byrow = T)
  # nu_mat[,1] <- 0
  nu <- c(0,nu_mat[1,])
  nu2 <- c(0,nu_mat[2,])
  return(list(nu,nu2))
}

ViterbiInd <- function(ind){
  
  
  clust_i <- which.max(re_prob[ind,])
  tran_list_clust <- tran_list[[clust_i]]
  
  emit_act_week <- emit_act[,,,1]
  emit_act_weekend <- emit_act[,,,2]
  
  emit_light_week <- emit_light[,,,1]
  emit_light_weekend <- emit_light[,,,2]
  
  vcovar_vec <- vcovar_mat[,ind]
  
  
  log_class_0_week <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_week[1,1,clust_i], 
                                          emit_act_week[1,2,clust_i], 
                                          emit_light_week[1,1,clust_i],
                                          emit_light_week[1,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,1,1], lintegral_mat[clust_i,1,1])
  
  log_class_1_week <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_week[2,1,clust_i], 
                                          emit_act_week[2,2,clust_i], 
                                          emit_light_week[2,1,clust_i],
                                          emit_light_week[2,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,2,1], lintegral_mat[clust_i,2,1])  
  
  log_class_0_weekend <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_weekend[1,1,clust_i], 
                                          emit_act_weekend[1,2,clust_i], 
                                          emit_light_weekend[1,1,clust_i],
                                          emit_light_weekend[1,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,1,2], lintegral_mat[clust_i,1,2])
  
  log_class_1_weekend <- logClassificationC( act[,ind], light[,ind], 
                                          emit_act_weekend[2,1,clust_i], 
                                          emit_act_weekend[2,2,clust_i], 
                                          emit_light_weekend[2,1,clust_i],
                                          emit_light_weekend[2,2,clust_i], 
                                          lod_act, lod_light, corr_mat[clust_i,2,2], lintegral_mat[clust_i,2,2])
  
  
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


readCpp( "cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )

################## EM Setup ################## 

###### True Settings ###### 
day_length <- 96 * 2
num_of_people <- 1000
missing_perc <- 0

init_true <- matrix(NA,ncol = 2,nrow = mix_num)
init_true[,1] <- seq(.1,.9,length.out = mix_num)
init_true[,2] <- 1 - init_true[,1]

params_tran_week_true <- matrix(rep(c(-3,-.5,-.7,.2,1.6,.4),mix_num),ncol = 6,byrow = T) 
params_tran_weekend_true <- matrix(rep(c(-5,3.5,2.2,-1,.4,-.4),mix_num),ncol = 6,byrow = T) 
params_tran_array_dim <- c(mix_num,6,vcovar_num)
params_tran_array_true <- array(NA,dim = params_tran_array_dim)
params_tran_array_true[,,1] <- params_tran_week_true 
params_tran_array_true[,,2] <- params_tran_weekend_true 

emit_act_week_true <- array(NA, c(2,2,mix_num))
emit_act_week_true[1,1,] <- seq(4,5,length.out = mix_num)
emit_act_week_true[1,2,] <- seq(1,2,length.out = mix_num)
emit_act_week_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_act_week_true[2,2,] <- seq(2,1,length.out = mix_num)

emit_act_weekend_true <- array(NA, c(2,2,mix_num))
emit_act_weekend_true[1,1,] <- seq(5,7,length.out = mix_num)
emit_act_weekend_true[1,2,] <- seq(2,3,length.out = mix_num)
emit_act_weekend_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_act_weekend_true[2,2,] <- seq(2,1,length.out = mix_num)

emit_act_true <- array(NA, c(2,2,mix_num,2))
emit_act_true[,,,1] <- emit_act_week_true
emit_act_true[,,,2] <- emit_act_weekend_true


emit_light_week_true <- array(NA, c(2,2,mix_num))
emit_light_week_true[1,1,] <- seq(6,7,length.out = mix_num)
emit_light_week_true[1,2,] <- seq(1,2,length.out = mix_num)
emit_light_week_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_light_week_true[2,2,] <- seq(1,2,length.out = mix_num)

emit_light_weekend_true <- array(NA, c(2,2,mix_num))
emit_light_weekend_true[1,1,] <- seq(7,9,length.out = mix_num)
emit_light_weekend_true[1,2,] <- seq(2,3,length.out = mix_num)
emit_light_weekend_true[2,1,] <- seq(1,0,length.out = mix_num)
emit_light_weekend_true[2,2,] <- seq(1,2,length.out = mix_num)

emit_light_true <- array(NA, c(2,2,mix_num,2))
emit_light_true[,,,1] <- emit_light_week_true
emit_light_true[,,,2] <- emit_light_weekend_true

corr_mat_true <- array(NA, c(mix_num,2,2))
corr_mat_true[,1,1] <- seq(.8,-.6,length.out = mix_num)
corr_mat_true[,2,1] <- seq(.6,-.8,length.out = mix_num)
corr_mat_true[,1,2] <- seq(.7,-.8,length.out = mix_num)
corr_mat_true[,2,2] <- seq(.5,-.9,length.out = mix_num)

beta_vec_true <- c(0,seq(-.75,1,length.out = mix_num-1))
beta_age_true <- 0.025

nu_true <- c(0,seq(-.4,.5,length.out = (mix_num-1)))
nu2_true <- c(0,seq(.05,-.15,length.out = (mix_num-1)))
# nu2_true <- 0


init <- init_true
params_tran_array <- params_tran_array_true
emit_act <- emit_act_true
emit_light <- emit_light_true
corr_mat <- corr_mat_true
beta_vec <- beta_vec_true
beta_age <- beta_age_true
nu <- nu_true
nu2 <- nu2_true

###### Simulate Data ###### 
if (!real_data){
  lod_act_true <- 0
  lod_light_true <- 0
  
  lod_act <- lod_act_true
  lod_light <- lod_light_true
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                               init=init_true,params_tran_array = params_tran_array_true,
                               emit_act = emit_act_true,emit_light = emit_light_true,
                               corr_mat = corr_mat_true,
                               lod_act = lod_act_true,lod_light = lod_light_true,
                               nu = nu_true,nu2 = nu2_true,beta_age_true = beta_age_true,
                               missing_perc = missing_perc, beta_vec_true = beta_vec_true)
  mc <- simulated_hmm[[1]]
  act <- simulated_hmm[[2]]
  light <- simulated_hmm[[3]]
  mixture_mat <- simulated_hmm[[4]]
  age_vec <- simulated_hmm[[5]]
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
  lod_act <- min(act[act>0],na.rm = T)
  act0 <- act == 0
  act <- log(act)
  act[act0] <- log(lod_act * .99)
  
  light_G <- wave_data_G[[2]]
  light_H <- wave_data_H[[2]]
  light <- rbind(light_G,light_H)
  light <- t(light[,-1])
  lod_light <- min(light[light>0],na.rm = T)
  light0 <- light == 0
  light <- log(light)
  light[light0] <- log(lod_light * .99)
  
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
  
  fcovar_num <- max(id$age_disc,na.rm = T)
  
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
  fcovar_mat <- Vec2Mat(fcovar_vec)
  surv_event <- surv_event[to_keep_inds]
  surv_time <- surv_time[to_keep_inds]
  
  
  day_length <- dim(act)[1]
  num_of_people <- dim(act)[2]
  
  age_vec <-id$age
  
  
}

###### Initial Settings ###### 

runif_tol <- .3

init <- matrix(rep(.5,mix_num*2),ncol = 2)
params_tran_array <- params_tran_array_true + runif(unlist(length(params_tran_array_true)),-runif_tol,runif_tol)



emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-runif_tol,runif_tol)
emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-runif_tol,runif_tol)


corr_mat <- corr_mat_true #+ runif(length(unlist(corr_mat_true)),-runif_tol,runif_tol)
beta_vec <- beta_vec + runif(mix_num,-runif_tol,runif_tol)
beta_vec[1] <- 0
beta_age <- 0

nu <- numeric(mix_num)
nu2 <- numeric(mix_num)

log_sweights_vec <- numeric(dim(act)[2])
time_vec <- c()
pi_l <- CalcPi(nu,nu2,age_vec)
re_prob <- pi_l

##########

mod_name <- paste0("JMHMM",mix_num,"SeedNA.rda")
load(mod_name)
init <- to_save[[2]][[1]]
params_tran_array <- to_save[[2]][[2]]
emit_act <- to_save[[2]][[3]]
emit_light <- to_save[[2]][[4]]
corr_mat <- to_save[[2]][[5]]
nu <- to_save[[2]][[6]]
nu2 <- to_save[[2]][[7]]
beta_vec <- to_save[[2]][[8]]
beta_age <- to_save[[2]][[9]]
re_prob <- to_save[[2]][[11]]



#########


################## EM ##################
bhaz_vec <- CalcBLHaz(beta_age,beta_vec,re_prob,surv_event,surv_time,age_vec)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)

tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
         emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
         lod_act = lod_act, lod_light = lod_light, 
         corr_mat = corr_mat, beta_vec = beta_vec, beta_age = beta_age,
         event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
         lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
         age_vec = age_vec, vcovar_mat = vcovar_mat)

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                  emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
                  emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
                  lod_act = lod_act, lod_light =  lod_light, 
                  corr_mat = corr_mat,lintegral_mat = lintegral_mat,
                  vcovar_mat = vcovar_mat)
         

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
iter_count <- 1

while(abs(like_diff) > 1e-3){

  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####

  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)


  params_tran_array <- CalcTranC(alpha,beta,act,light,params_tran_array,
                                      emit_act,emit_light,corr_mat,
                                      pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat)

  ##### Weights  #####
  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])
  
  
  
  ##### Mixing Proportion  #####
  if (beta_bool){
    nu_list  <- CalcNu(nu,nu2,re_prob,age_vec)
    # nu_list  <- CalcNuNum(nu,nu2)
    nu <- nu_list[[1]]
    nu2 <- nu_list[[2]]

    pi_l <- CalcPi(nu,nu2,age_vec)
  }
  
  re_prob <- CalcProbRE(alpha,pi_l)
    
  
  ##### Bivariate Normal Est  #####
  
  ##### Mixture Normal Param
  
  corr_mat[,1,] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  corr_mat[,2,] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  emit_act[1,1,,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_slee,vcovar_mat+1)
  emit_act[2,1,,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  emit_light[1,1,,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  emit_light[2,1,,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  emit_act[1,2,,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  #HAD ISSUES WITH THIS VAR
  emit_act[2,2,,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  emit_light[1,2,,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
  emit_light[2,2,,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

  
  ##### Survival ####
  if(beta_bool){
  
    beta_vec_long <- c(beta_age,beta_vec[-1])

    beta_vec_long <- CalcBeta(beta_vec_long)
    beta_age <- beta_vec_long[1]
    beta_vec <- c(0,beta_vec_long[-1])

    bhaz_vec <- CalcBLHaz(beta_age,beta_vec,re_prob,surv_event,surv_time,age_vec)
    bline_vec <- bhaz_vec[[1]]
    cbline_vec <- bhaz_vec[[2]]
  }
  
  
  ##### Reorder #####
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  if (mix_num > 1){
    # reord_inds <- c(1,order(beta_mat[2:mix_num,1])+1)
    # emit_act <- emit_act[,,reord_inds]
    # emit_light <- emit_light[,,reord_inds]
    # pi_l <- pi_l[reord_inds]
    # re_prob <- re_prob[,reord_inds]
    # params_tran_array <- params_tran_array[reord_inds,,]
    # dim(params_tran_array) <- params_tran_array_dim
    # corr_mat <- corr_mat[reord_inds,]
    # init <- init[reord_inds,]
    # beta_mat <- beta_mat[reord_inds,]
    # if(!real_data){mixture_mat <- mixture_mat[,reord_inds]}
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
  
  ##### #####

  # Calculate tran_list after reordering 
  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
                    emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
                    lod_act = lod_act, lod_light = lod_light, 
                    corr_mat = corr_mat, beta_vec = beta_vec, beta_age = beta_age,
                    event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                    lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]),
                    age_vec = age_vec, vcovar_mat = vcovar_mat)
  
  beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                    emit_act_week = emit_act[,,,1],emit_light_week = emit_light[,,,1],
                    emit_act_weekend = emit_act[,,,2],emit_light_weekend = emit_light[,,,2],
                    lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat,
                    vcovar_mat = vcovar_mat)
  
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",mix_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  # if(iter_count > 5 & min(surv_event %*% re_prob) > 1e-10){beta_bool <- T}
  if(iter_count > 3){beta_bool <- T}
  
  
  
  iter_count <- iter_count + 1
  
  if (iter_count %% 20 == 0){
    
    tran_df <- ParamsArray2DF(params_tran_array)
    
    if (!real_data){
      tran_df_true <- ParamsArray2DF(params_tran_array_true)
      
      tran_df <- tran_df %>% mutate(truth = tran_df_true[,1])
      tran_df <- tran_df %>% mutate(resid = prob - truth)
    }
    
    
    true_params <- list(init_true,params_tran_array_true,
                        emit_act_true,emit_light_true,
                        corr_mat_true,nu_true,nu2_true,
                        beta_vec_true,beta_age_true)
    
    est_params <- list(init,params_tran_array,
                       emit_act,emit_light,
                       corr_mat,nu,nu2,
                       beta_vec,beta_age,
                       tran_df,
                       re_prob,
                       weights_array_list)
    
    bic <- CalcBIC(new_likelihood,mix_num,act,light)
    to_save <- list(true_params,est_params,bic)
    
    save(to_save,file = paste0("InterJMHMM",mix_num,"Seed",sim_num,".rda"))
  }
}


tran_df <- ParamsArray2DF(params_tran_array)
decoded_mat <- sapply(c(1:num_of_people), ViterbiInd)


if (!real_data){
  tran_df_true <- ParamsArray2DF(params_tran_array_true)
  
  tran_df <- tran_df %>% mutate(truth = tran_df_true[,1])
  tran_df <- tran_df %>% mutate(resid = prob - truth)
}


true_params <- list(init_true,params_tran_array_true,
                    emit_act_true,emit_light_true,
                    corr_mat_true,nu_true,nu2_true,
                    beta_vec_true,beta_age_true)

est_params <- list(init,params_tran_array,
                   emit_act,emit_light,
                   corr_mat,nu,nu2,
                   beta_vec,beta_age,
                   tran_df,
                   re_prob,
                   weights_array_list,
                   decoded_mat)

# init <- to_save[[2]][[1]]
# params_tran_array <- to_save[[2]][[2]]
# emit_act <- to_save[[2]][[3]]
# emit_light <- to_save[[2]][[4]]
# corr_mat <- to_save[[2]][[5]]
# nu <- to_save[[2]][[6]]
# nu2 <- to_save[[2]][[7]]
# beta_vec <- to_save[[2]][[8]]
# beta_age <- to_save[[2]][[9]]
# re_prob <- to_save[[2]][[11]]

bic <- CalcBIC(new_likelihood,mix_num,act,light)
to_save <- list(true_params,est_params,bic)


if (like_diff > -1){
  save(to_save,file = paste0("JMHMM",mix_num,"Seed",sim_num,".rda"))
}
