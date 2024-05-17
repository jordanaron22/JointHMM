################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)

RE_type <- "norm"

real_data <- F
set_seed <- T

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 99}
if (set_seed){set.seed(sim_num)}

RE_num <- as.numeric(commandArgs(TRUE)[1])
sim_size <- as.numeric(commandArgs(TRUE)[2])
RE_type <- as.character(commandArgs(TRUE)[3])
print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",RE_num))

if(is.na(RE_num)){RE_num <- 2}
if(is.na(sim_size)){sim_size <- 0}
if(is.na(RE_type)){RE_type <- "norm"}

################## Functions ################## 

sourceCpp( "cFunctions.cpp" )

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

Params2TranVectorT <- function(index,len,params_tran){
  return(t(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,index=index)))
}

TranByTimeVec <- function(index, params_tran, time_vec){
  return(lapply(time_vec, Params2Tran, params_tran = params_tran,index=index))
}

Param2TranHelper <- function(p12,p21){
  tran <- matrix(0,2,2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1- tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  return(tran)
}

Params2Tran <- function(params_tran,time,index){
  
  harmonic_ind <- ((index*2)-1)+6
  param_matrix <- matrix(params_tran,ncol=18,nrow=2, byrow = T)
  if (index == 1){
    tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,7]*cos(2*pi*time/96)+param_matrix[1,8]*sin(2*pi*time/96),
                             param_matrix[2,1]+param_matrix[2,7]*cos(2*pi*time/96)+param_matrix[2,8]*sin(2*pi*time/96))
  } else {
    tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,index]+
                               (param_matrix[1,7]+param_matrix[1,harmonic_ind])*cos(2*pi*time/96)+
                               (param_matrix[1,8]+param_matrix[1,harmonic_ind+1])*sin(2*pi*time/96),
                             param_matrix[2,1]+param_matrix[2,index]+
                               (param_matrix[2,7]+param_matrix[2,harmonic_ind])*cos(2*pi*time/96)+
                               (param_matrix[2,8]+param_matrix[2,harmonic_ind+1])*sin(2*pi*time/96))
  }
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

ChooseTran <- function(covar_tran_bool){
  covar_ind <- which(covar_tran_bool== 1)
  if (length(covar_ind) == 2){
    return(covar_ind[2])
  } else {
    return(1)
  }
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

SimulateHMM <- function(day_length,num_of_people,init,params_tran,emit_act,
                        emit_light,corr_mat,lod_act,lod_light, pi_l){

  
  mixture_num <- dim(emit_act)[3]

  
  
  # covar_mat <- t(rmultinom(num_of_people,1,rep(1/mixture_num,mixture_num)))
  covar_mat <- t(rmultinom(num_of_people,1,pi_l_true))
  covar_mat <- cbind((numeric(num_of_people) + 1),covar_mat[,2:mixture_num])
  
  tran_list <- lapply(c(1:dim(covar_mat)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:day_length))
  covar_vec <- apply(covar_mat,1,ChooseTran)
  
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

  
  
  
  light_matrix[light_matrix<lod_light] <- lod_light
  activity_matrix[activity_matrix<lod_act] <- lod_act
  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,covar_mat))
}

CondMarginalize <- function(alpha,beta,pi_l,wakesleep){
  
  if(wakesleep == "wake"){wakesleep_ind <- 1}
  if(wakesleep == "sleep"){wakesleep_ind <- 2}
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      # alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(re_weights[ind,re_ind])
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[re_ind])
    }
  }
  
  ind_like_mat <- apply(alpha_beta,c(1,4),logSumExp)
  
  weight_array <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  for (ind in 1:dim(alpha_beta)[4]){
    for (t in 1:dim(alpha_beta)[1]){
      weight_array[t,ind,] <- alpha_beta[t,wakesleep_ind,,ind] - ind_like_mat[t,ind]
    }
  }
  
  return(weight_array)
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
  
  gradient <- matrix(0,2,18)
  hessian1 <- matrix(0,18,18)
  hessian2 <- matrix(0,18,18)
  hessian <- matrix(0,36,36)
  hessian_vec <- matrix(0,2,18)
  cos_part_vec <- matrix(0,2,6)
  sin_part_vec <- matrix(0,2,6)
  cos_sin_part <- matrix(0,2,6)
  
  
  # tran_mat <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  tran_list_mat <- lapply(c(1:dim(covar_mat_tran)[2]),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  tran_vals_re_00 <- CalcTranHelperC(init_state = 0,new_state = 0,act = act,
                                     light = light,tran_list_mat = tran_list_mat,tran_ind_vec = tran_ind_vec,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_01 <- CalcTranHelperC(init_state = 0,new_state = 1,act = act,
                                     light = light,tran_list_mat = tran_list_mat,tran_ind_vec = tran_ind_vec,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_10 <- CalcTranHelperC(init_state = 1,new_state = 0,act = act,
                                     light = light,tran_list_mat = tran_list_mat,tran_ind_vec = tran_ind_vec,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  tran_vals_re_11 <- CalcTranHelperC(init_state = 1,new_state = 1,act = act,
                                     light = light,tran_list_mat = tran_list_mat,tran_ind_vec = tran_ind_vec,
                                     emit_act = emit_act,emit_light = emit_light,ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                     corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals_re <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals_re <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals_re <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals_re <- tran_vals_re_11}
      
      
      tran_vals <- apply(tran_vals_re, c(1,2), sum)
      
      for (ind in 1:length(alpha)){
        
        tran_ind <- tran_ind_vec[ind]
        harmonic_ind <- ((tran_ind*2)-1)+6
        tran_vec <- tran_list_mat[[tran_ind]]
        
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
        
        #Maybe change this back to 1:len-1
        cos_vec <- cos(2*pi*c(2:(len))/96)
        sin_vec <- sin(2*pi*c(2:(len))/96)
        
        
        gradient[init_state,tran_ind] <- gradient[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime)
        gradient[init_state,harmonic_ind] <- gradient[init_state,harmonic_ind] + sum(tran_vals[,ind]*tran_prime*cos_vec)
        gradient[init_state,harmonic_ind+1] <- gradient[init_state,harmonic_ind+1] + sum(tran_vals[,ind]*tran_prime*sin_vec)
        
        if (tran_ind != 1){
          gradient[init_state,1] <- gradient[init_state,1] + sum(tran_vals[,ind]*tran_prime)
          gradient[init_state,7] <- gradient[init_state,7] + sum(tran_vals[,ind]*tran_prime*cos_vec)
          gradient[init_state,8] <- gradient[init_state,8] + sum(tran_vals[,ind]*tran_prime*sin_vec)
        }
        
        
        
        
        hessian_vec[init_state,tran_ind] <- hessian_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime)
        
        hessian_vec[init_state,harmonic_ind] <- hessian_vec[init_state,harmonic_ind] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec^2)
        hessian_vec[init_state,harmonic_ind+1] <- hessian_vec[init_state,harmonic_ind+1] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec^2)
        
        cos_part_vec[init_state,tran_ind] <- cos_part_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec)
        sin_part_vec[init_state,tran_ind] <- sin_part_vec[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec)
        
        cos_sin_part[init_state,tran_ind] <- cos_sin_part[init_state,tran_ind] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec*sin_vec)
        
        if (tran_ind != 1){
          hessian_vec[init_state,1] <- hessian_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime)
          
          hessian_vec[init_state,7] <- hessian_vec[init_state,7] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec^2)
          hessian_vec[init_state,8] <- hessian_vec[init_state,8] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec^2)
          
          cos_part_vec[init_state,1] <- cos_part_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec)
          sin_part_vec[init_state,1] <- sin_part_vec[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime*sin_vec)
          
          cos_sin_part[init_state,1] <- cos_sin_part[init_state,1] + sum(tran_vals[,ind]*tran_prime_prime*cos_vec*sin_vec)
        }
        
        
        
        
        
      }
      
      
    }
  }
  
  gradient <- as.vector(t(gradient))
  
  #### HESS 1
  
  diag(hessian1) <- c(hessian_vec[1,])
  hessian1[1,1:6] <- c(hessian_vec[1,1:6])
  hessian1[1:6,7] <- cos_part_vec[1,]
  hessian1[1:6,8] <- sin_part_vec[1,]
  hessian1[1,c(7,9,11,13,15,17)] <- cos_part_vec[1,]
  hessian1[1,c(8,10,12,14,16,18)] <- sin_part_vec[1,]
  hessian1[2,c(9,10)] <- c(cos_part_vec[1,2],sin_part_vec[1,2])
  hessian1[3,c(11,12)] <- c(cos_part_vec[1,3],sin_part_vec[1,3])
  hessian1[4,c(13,14)] <- c(cos_part_vec[1,4],sin_part_vec[1,4])
  hessian1[5,c(15,16)] <- c(cos_part_vec[1,5],sin_part_vec[1,5])
  hessian1[6,c(17,18)] <- c(cos_part_vec[1,6],sin_part_vec[1,6])
  
  hessian1[7,c(7,9,11,13,15,17)] <- hessian_vec[1,c(7,9,11,13,15,17)]
  hessian1[7,c(8,10,12,14,16,18)] <- cos_sin_part[1,]
  
  hessian1[8,c(9,11,13,15,17)] <- cos_sin_part[1,2:6]
  hessian1[8,c(8,10,12,14,16,18)] <- hessian_vec[1,c(8,10,12,14,16,18)]
  
  hessian1[7,8] <- cos_sin_part[1,1]
  hessian1[9,10] <- cos_sin_part[1,2]
  hessian1[11,12] <- cos_sin_part[1,3]
  hessian1[13,14] <- cos_sin_part[1,4]
  hessian1[15,16] <- cos_sin_part[1,5]
  hessian1[17,18] <- cos_sin_part[1,6]
  
  hessian1 <- hessian1+t(hessian1)
  diag(hessian1) <- diag(hessian1)/2
  
  #### HESS 2
  
  diag(hessian2) <- c(hessian_vec[2,])
  hessian2[1,1:6] <- c(hessian_vec[2,1:6])
  hessian2[1:6,7] <- cos_part_vec[2,]
  hessian2[1:6,8] <- sin_part_vec[2,]
  hessian2[1,c(7,9,11,13,15,17)] <- cos_part_vec[2,]
  hessian2[1,c(8,10,12,14,16,18)] <- sin_part_vec[2,]
  hessian2[2,c(9,10)] <- c(cos_part_vec[2,2],sin_part_vec[2,2])
  hessian2[3,c(11,12)] <- c(cos_part_vec[2,3],sin_part_vec[2,3])
  hessian2[4,c(13,14)] <- c(cos_part_vec[2,4],sin_part_vec[2,4])
  hessian2[5,c(15,16)] <- c(cos_part_vec[2,5],sin_part_vec[2,5])
  hessian2[6,c(17,18)] <- c(cos_part_vec[2,6],sin_part_vec[2,6])
  
  hessian2[7,c(7,9,11,13,15,17)] <- hessian_vec[2,c(7,9,11,13,15,17)]
  hessian2[7,c(8,10,12,14,16,18)] <- cos_sin_part[2,]
  
  hessian2[8,c(9,11,13,15,17)] <- cos_sin_part[2,2:6]
  hessian2[8,c(8,10,12,14,16,18)] <- hessian_vec[2,c(8,10,12,14,16,18)]
  
  
  hessian2[7,8] <- cos_sin_part[2,1]
  hessian2[9,10] <- cos_sin_part[2,2]
  hessian2[11,12] <- cos_sin_part[2,3]
  hessian2[13,14] <- cos_sin_part[2,4]
  hessian2[15,16] <- cos_sin_part[2,5]
  hessian2[17,18] <- cos_sin_part[2,6]
  
  hessian2 <- hessian2+t(hessian2)
  diag(hessian2) <- diag(hessian2)/2
  
  #######
  
  
  hessian[1:18,1:18] <- hessian1
  hessian[19:36,19:36] <- hessian2
  
  grad <- -gradient
  hess <- -hessian
  
  if (!real_data){
    null_inds <- grad != 0
    params_tran[null_inds] <- params_tran[null_inds] - solve(hess[null_inds,null_inds],grad[null_inds])
  } else {
    params_tran <- params_tran - solve(hess,grad)
  }
  
  return(params_tran)
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

################## EM Setup ################## 

obs_per_day <- 96

init_true <- matrix(c(.15,.85,.3,.7),2,2,byrow = T)
params_tran_true <- c(-3,.5,-.25,0,0,0,.9,-.8,-1,.1,.2,1.5,0,0,0,0,0,0,
                      -2.2,.3,-.3,0,0,0,-.75,.8,.75,.2,.3,.75,0,0,0,0,0,0)

emit_act_true <- array(NA, c(2,2,2))
emit_act1 <- matrix(c(1.5,.8,0,.5),2,2,byrow = T)
emit_act2 <- matrix(c(2.5,1.2,.25,.75),2,2,byrow = T)
emit_act_true[,,1] <- emit_act1
emit_act_true[,,2] <- emit_act2



emit_light_true <- array(NA, c(2,2,2))
emit_light1 <- matrix(c(7,2.5,2,1.5),2,2,byrow = T)
emit_light2 <- matrix(c(8,3,3,2),2,2,byrow = T)
emit_light_true[,,1] <- emit_light1
emit_light_true[,,2] <- emit_light2

corr_mat_true <- matrix(c(.6,.8,-.5,-.7),2,2,byrow = T)

lod_act_true <- 0
lod_light_true <- 0
pi_l_true <- c(1/3,2/3)
simulated_hmm <- SimulateHMM(day_length = obs_per_day,num_of_people = 1000,
                              init=init_true,params_tran = params_tran_true,
                              emit_act = emit_act_true,emit_light = emit_light_true,corr_mat = corr_mat_true,
                              lod_act = lod_act_true,lod_light = lod_light_true,pi_l = pi_l_true)

mc <- simulated_hmm[[1]]
act <- simulated_hmm[[2]]
light <- simulated_hmm[[3]]
covar_mat_tran <- simulated_hmm[[4]]
tran_ind_vec <- apply(covar_mat_tran,1,ChooseTran)



init <- init_true
params_tran <- params_tran_true

tran_list <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:96))


emit_act <- emit_act_true

emit_light <- emit_light_true 

corr_mat <- corr_mat_true
lod_act <- lod_act_true
lod_light <- lod_light_true

log_sweights_vec <- numeric(dim(act)[2])
time_vec <- c()
pi_l <- pi_l_true


# emit_act <- emit_act_true
# dim(emit_act)[3] <- 1
# emit_light <- emit_light_true
# dim(emit_light)[3] <- 1
# pi_l <- c(1)


################## EM ################## 

lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)

alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         tran_ind_vec = tran_ind_vec,lod_act = lod_act, lod_light =  lod_light, 
         corr_mat = corr_mat,lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]))

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
              emit_act = emit_act,emit_light = emit_light,
              tran_ind_vec = tran_ind_vec,lod_act = lod_act, lod_light =  lod_light, 
              corr_mat = corr_mat,lintegral_mat = lintegral_mat)
         

new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)

while(abs(like_diff) > 1e-3*1e-10){
  
  start_time <- Sys.time()
  likelihood <- new_likelihood

  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)

  params_tran <- CalcTranC(alpha,beta,act,light,params_tran,emit_act,emit_light,corr_mat,covar_mat_tran,pi_l,lod_act,lod_light,lintegral_mat)
  tran_list <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:obs_per_day))
  
  ################## Weights
  #Weights are prob currently in the wake state
  weights_array_wake <- exp(CondMarginalize(alpha,beta,pi_l,"wake"))
  weights_mat_wake <-rowSums(weights_array_wake, dims = 2)
  weights_vec_wake <- as.vector(weights_mat_wake)

  weights_array_sleep <- exp(CondMarginalize(alpha,beta,pi_l,"sleep"))
  weights_mat_sleep <-rowSums(weights_array_sleep, dims = 2)
  weights_vec_sleep <- as.vector(weights_mat_sleep)
  
  ################## RE
  
  re_prob <- CalcProbRE(alpha,pi_l)
  pi_l <- colSums(re_prob)/dim(act)[2]
  pi_l[pi_l<1e-100] <- 1e-100
  if (RE_num <= 1){pi_l <- c(1)}
  
  ################## Reorder
  
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  if (length(pi_l)>1){
    reord_inds <- order(emit_act[1,1,])
    emit_act <- emit_act[,,reord_inds]
    emit_light <- emit_light[,,reord_inds]
    pi_l <- pi_l[reord_inds]
  }
  
  # # Uncomment to use real mc states as weights
  # weights_array_wake <- 1-mc
  # weights_array_sleep <- mc
  # dim(weights_array_wake)[3] <- 1
  # dim(weights_array_sleep)[3] <- 1
  
  
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
  
  
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  
  alpha <- ForwardC(act = act,light = light,
                    init = init,tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    tran_ind_vec = tran_ind_vec,lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat,log_sweight = numeric(dim(act)[2]))
  
  beta <- BackwardC(act = act,light = light, tran_list = tran_list,
                    emit_act = emit_act,emit_light = emit_light,
                    tran_ind_vec = tran_ind_vec,lod_act = lod_act, lod_light =  lod_light, 
                    corr_mat = corr_mat,lintegral_mat = lintegral_mat)
  
  
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",RE_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  
}






