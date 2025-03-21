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

######## Set controls

#int (0-2), sets simulation size
sim_size <- 2
#bool, set T to set seed
set_seed <- 1
#bool, set T to use tobit approach
tobit <- 1
#bool, set T to include activity
incl_act <- 1
#bool, set T to run F-B within transition calculation to ensure likelihood increase
check_tran <- 0
#bool, set T to remove week data
weekend_only <- 0
#int, set to number of true clusters
misspecification <- 0
#bool, set T save less data and save space
save_space <- 0
#int, 0 -> all days, 1-7 -> subset to that day of the week only
single_day <- 0
#bool, set T only run survival portion
run_only_surv <- 0
#bool, set T include MIMs summary measures from Leroux in Cox model
tlogmims_bool <- 0
#bool, set T to subset data, currently age>=70
subset_data <- 1
#int, number of periods in a day
#24 hourly, 96 fifteen min, 1440 single min
period_len <- 24

vcovar_num <- 2

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(paste("Seed",sim_num))

if (is.na(sim_num)){sim_num <- 1}
if (set_seed){set.seed(sim_num)}

print("Command line arguments:")
print(commandArgs(T))


####### Set bash defined controls
#int, number of mixtures for MHMM
mix_num <- as.numeric(commandArgs(TRUE)[1])
print(paste("mixnum",mix_num))

#0 for two stage, 2 for JM
incl_surv <- as.numeric(commandArgs(TRUE)[2])
print(paste("Include Surv",incl_surv))

#0 for simulation, 1 for NHANES
real_data <- as.numeric(commandArgs(TRUE)[3])
print(paste("Real Data",real_data))

#0 for standard, 1 for bootstrap
bootstrap <- as.numeric(commandArgs(TRUE)[4])
print(paste("Bootstrap",bootstrap))

#0 for no randomization, set to max runif val to add (at most 1)
randomize_init <- as.numeric(commandArgs(TRUE)[5])
print(paste("Random start",randomize_init))

#0 for standard, 1 for leave out CV
leave_out <- as.numeric(commandArgs(TRUE)[6])
print(paste("Leave out",leave_out))

#0 for standard, 1 to load data for hot start
load_data <- as.numeric(commandArgs(TRUE)[7])
print(paste("load data",load_data))

#1 to include light data, 0 to exclude
incl_light <- as.numeric(commandArgs(TRUE)[8])
print(paste("include light",incl_light))

#0 for standard, 1 to only use predicted HMM state for LOCV
wake_sleep <- as.numeric(commandArgs(TRUE)[9])
print(paste("Wake/sleep only for leaveout",wake_sleep))



##### set bash controls if not actually running bash 
##### change these if not running on cluster
if(is.na(mix_num)){mix_num <- 2 }
if(is.na(incl_surv)){incl_surv <-0}
if(is.na(real_data)){real_data <-1}
if(is.na(bootstrap)){bootstrap <- 0}
if(is.na(randomize_init)){randomize_init <-.2}
if(is.na(load_data)){load_data <- 1}
if(is.na(leave_out)){leave_out <- 0}
if(is.na(incl_light)){incl_light <- 1}
if(is.na(wake_sleep)){wake_sleep <- 0}

#Start by not estimation survival coef for stability, unless doing hot start
if (load_data | leave_out){
  beta_bool <- 1
} else {
  beta_bool <- 0
}
#overides leave_out
if (bootstrap){beta_bool <- 0}

#switch to order clusters from best to worst survival
relabel_reset <- F



print(paste("Sim Seed:",sim_num,"HMM Num:",mix_num))

#compiles full model nage to individualize save files
model_name <- "JMHMM"
if (!real_data){model_name <- paste0(model_name,"SimSize",sim_size)}
if (bootstrap){model_name <- paste0(model_name,"Bootstrap")}
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
if (subset_data){model_name <- paste0(model_name,"Subset",period_len)}

model_name <- paste0(model_name,"Mix",mix_num,"Seed",sim_num,".rda")
print(model_name)
################## Functions ##################
#Reads in rcpp file
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


#similar to below but transposes and compiles output into single matrix
#needed for faster derivations
Params2TranVectorTresid <- function(re_ind,len,params_tran){
  return(t(sapply(c(1:(len)),FUN = Params2Tran,params_tran = params_tran,re_ind=re_ind)))
}

#Calculates list of transition probabilities over all times
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

#takes vector of transition parameters and outputs transition matrix given current time and mixture
Params2Tran <- function(params_tran,time,re_ind){
  
  param_matrix <- matrix(params_tran[re_ind,],ncol=3,nrow=2, byrow = T)
  tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,2]*cos(2*pi*time/period_len)+param_matrix[1,3]*sin(2*pi*time/period_len),
                           param_matrix[2,1]+param_matrix[2,2]*cos(2*pi*time/period_len)+param_matrix[2,3]*sin(2*pi*time/period_len))


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

#Calcuates case where both activity and light are below LoD
Case4 <- function(act_obs,mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD){
  
  mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs)
  sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
  
  lognorm_dens <- dnorm(act_obs,mu_act,sig_act) * 
    pnorm(light_LOD,mu_light_cond,sig_light_cond)
  return(lognorm_dens)
}

#Calculates case 4 for given normal parameters and store in matrix for access later
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
  
  #work on log scale so -9999 is effectively -Inf
  lintegral_mat[lintegral_mat == -Inf] <- -9999
  
  return(lintegral_mat)
}

#Simulates hidden states
SimulateMC <- function(day_length,init,tran_list_ind,mixture_ind,vcovar_vec){
  hidden_states <- numeric(day_length)
  
  for (i in 1:day_length){
    
    tran <- tran_list_ind[[vcovar_vec[i]]][[(i-1)%%period_len+1]]
    
    if (i == 1) {
      hidden_states[1] <- rbinom(1,1,init[mixture_ind,2])
    } else {
      hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
    }
  }
  
  
  return(hidden_states)
}

#used to simulate survival/censor times
finv <- function(lam,time, randu, xb_ind){
  return((1-pexp(time,rate = lam))^(exp(xb_ind)) - randu)
}

#simulates survival by transforming runif var into rexponential
#uses: S(t) = 1-F(t)
#finds root of survival time: 1 - F(t|X) - u
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

#Generates transition matrices across time as one large matrix for faster computation
#Outer index mixture and then by week/weekend
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

#Generates big list of all transition matrices
#Outer index by mixture then by week/weekend then by time
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

#Simulates data 
SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,
                        emit_act,emit_light,corr_mat,
                        lod_act,lod_light, nu_mat,beta_vec_true,beta_age_true,beta_covar_sim,
                        missing_perc,lambda_act_mat,lambda_light_mat,
                        misspecification){
  
  mix_num <- dim(emit_act)[3]
  
  #simulates age and single hypothetical categorical variable with 3 categories
  age_vec <- floor(runif(num_of_people,10,81))
  surv_covar_sim_wide <- t(rmultinom(num_of_people,1,c(1/3,1/3,1/3)))
  surv_covar_sim <- apply(surv_covar_sim_wide, 1, which.max)
  
  #simulates stationary and moderact activity and then calulates nu matrix
  statact_vec <- floor(runif(num_of_people,0,16))
  modact_vec <- rbinom(num_of_people,1,.5)
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)
  
  pi_l_true <- CalcPi(nu_mat,nu_covar_mat)
  
  #if misspecified model, generate data according to true model
  if (misspecification == 1){
    pi_l_true <- cbind(pi_l_true[,1],pi_l_true[,1])
    pi_l_true[,1] <- 1
    pi_l_true[,2] <- 0
  } else if (misspecification > 1){
    pi_l_true <- pi_l_true[,1:misspecification]
    pi_l_true <- pi_l_true/ rowSums(pi_l_true)
  }
  
  #edge case for 1 mixture model
  if (dim(pi_l_true)[2] > 1){
    mixture_vec <- rMultinom(pi_l_true,1)
  } else {
    mixture_vec <- matrix(1,nrow = dim(pi_l_true)[1],ncol = 1)
  }
  
  #simulates week/weekend divide
  #only works for specific lengths
  if (day_length == 192){
    vcovar_vec <- c(rep(0,48),rep(1,period_len),rep(0,48))
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
  
  
  #actually simulates data here
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)
    
    mixture_ind <- mixture_vec[ind]
    
    tran_vcovar_list <- tran_list[[mixture_ind]]
    vcovar_vec_ind <- vcovar_mat[,ind]+1
    
    hidden_states <- SimulateMC(day_length,init,tran_vcovar_list,mixture_ind,vcovar_vec_ind)
    
    for (i in 1:day_length){
      
      #sets normal parameters for current mixture and week status
      mu_act <- emit_act[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]] 
      sig_act <- emit_act[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]] 
      mu_light <- emit_light[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]] 
      sig_light <- emit_light[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]] 
      bivar_corr <- corr_mat[mixture_ind,hidden_states[i] + 1,vcovar_vec_ind[i]] 
      
      #non tobit is semi-defunct
      #if tobit generates correlated data
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
  
  #simulates survival data
  surv_list <- SimSurvival(mixture_vec,beta_vec_true,beta_age_true,beta_covar_sim,age_vec,surv_covar_sim)


  #removes values below LoD
  if (tobit){
    light_matrix[light_matrix<lod_light] <- lod_light
    activity_matrix[activity_matrix<lod_act] <- lod_act
  }
  
  #removes some data
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)
  
  activity_matrix[act_missing==1] <- NA
  light_missing[light_missing==1] <- NA
  

  
  return(list(hidden_states_matrix,activity_matrix,light_matrix,mixture_vec,
              age_vec,nu_covar_mat,vcovar_mat,surv_list,surv_covar_sim))
}

#turns F-B into simple prob of wake/sleep arrays
CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  #need to add prob of being in each mixture
  #not included directly in F-B
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[ind,re_ind])
    }
  }
  
  #individual likelihood of being in specific mixture
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

#calculates initial probability
CalcInit <- function(alpha, beta,pi_l,log_sweights_vec){
  
  #setup
  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_1_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_mat <- matrix(0,dim(pi_l)[2],2)
  
  #individual likelihood vector
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)
  
  
  for(ind in 1:length(alpha)){ 
    ind_like <- ind_like_vec[ind]
    
    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l[ind,]) - ind_like + log_sweights_vec[ind]
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l[ind,]) - ind_like + log_sweights_vec[ind]

  }
  
  #normalizes
  for (re_ind in 1:(dim(pi_l)[2])){
    init_0 <- logSumExp(init_0_vec[,re_ind])
    init_1 <- logSumExp(init_1_vec[,re_ind])
    init_vec <- exp(c(init_0,init_1) - logSumExp(c(init_0,init_1)))
    init_mat[re_ind,] <- init_vec
  }
  
  return(init_mat)
  
}

#calculates probability of an individual being in each mixture
CalcProbRE <- function(alpha,pi_l){
  
  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)
  
  
  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      #sums over latent states for last time and specific mixture
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[ind,re_ind])
    }
    #normalizes
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
    
  }
  
  return(re_weights)
  
}

#Helper function for transition calculation
#Relies on helper in C but coded in R to organize
CalcTranHelper <- function(act, light, tran_list_mat, emit_act, emit_light, 
                           ind_like_vec, alpha, beta, lod_act, lod_light, 
                           corr_mat, lintegral_mat, pi_l,vcovar_mat,
                           lambda_act_mat, lambda_light_mat, tobit){
  
  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  
  tran_vals_re_array <- array(NA,c(2,2,len - 1,num_people,num_re))
  
  #cpp doesnt have 4d arrays so need to organize our 4d array 3d arrays for week/weekend
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

#Makes matrix symmetric
Symmetricize <- function(mat){
  mat <- mat + t(mat)
  diag(mat) <- diag(mat)/2
  return(mat)
}

#Calculates gradient and hessian for transition probabilities within transition helper
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

#heavy lifting of LM for transition
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
  
  
  cos_vec <- cos(2*pi*c(2:(len))/period_len)
  sin_vec <- sin(2*pi*c(2:(len))/period_len)
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      tran_vals <- tran_vals_re_array[init_state,new_state,,,]
      if (mix_num  == 1){dim(tran_vals)[3] <- 1}
      
      for (ind in 1:length(alpha)){
        
        vcovar_vec <- vcovar_mat[-1,ind]
        #need to add 1 bc C is base 0
        vcovar_vecR <- vcovar_vec + 1
        
        for(re_ind in 1:mix_num){
          
          #calculates transition over week/weekend
          if (vcovar_num == 1){
            tran_mat <- tran_list_mat[[re_ind]][[1]]
          } else if (vcovar_num == 2){
            tran_mat_week <- tran_list_mat[[re_ind]][[1]]
            tran_mat_weekend <- tran_list_mat[[re_ind]][[2]]
            tran_mat <- tran_mat_week * (1-vcovar_vec) + tran_mat_weekend * vcovar_vec
          }
            
          
          if(init_state == 1 & new_state == 1){
            #Left these in for debugging
            #moved to putting all transition values into large matrix to vectorize
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

            #grad and hessian calculations
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

}

#levenberg marquardt
#interpolates btwn newton and gradient descent
LM <- function(grad_array,hess_array,params_tran_array,check_tran,likelihood,pi_l, step_size = .01){
  params_tran_array_new <- params_tran_array
  new_likelihood <- -Inf
  
  for (re_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      #returns all -1 if non-invertible
      inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind],grad_array[,re_ind,vcovar_ind])
      
      step_size <- 1
      #runs until step isnt too big
      #increases step size effectively increases hyperparam for gradient descent -> smaller grad step
      while(max(abs(inf_fact)) > 2 | all(inf_fact == -1)){
        step_fact <- matrix(0,6,6)
        diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) * step_size
        if(all(inf_fact == -1)){diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) + step_size}
        
        inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind]+step_fact,grad_array[,re_ind,vcovar_ind])
        step_size <- step_size * 10
      }
      params_tran_array_new[re_ind,,vcovar_ind] <- params_tran_array_new[re_ind,,vcovar_ind] - inf_fact
    }
  }
  
  #Like decrease may happen here
  if (check_tran){
    tran_list <- GenTranList(params_tran_array_new,c(1:day_length),mix_num,vcovar_num)
    alpha <- Forward(act = act,light = light,
                     init = init,tran_list = tran_list,
                     emit_act = emit_act,emit_light = emit_light,
                     lod_act = lod_act, lod_light = lod_light, 
                     corr_mat = corr_mat, beta_vec = beta_vec, surv_coef = surv_coef,surv_covar_risk_vec = surv_covar_risk_vec,
                     event_vec = surv_event, bline_vec = bline_vec, cbline_vec = cbline_vec,
                     lintegral_mat = lintegral_mat,log_sweight = log_sweights_vec,
                     surv_covar = surv_covar, vcovar_mat = vcovar_mat,
                     lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit,incl_surv = incl_surv)
    new_like <- CalcLikelihood(alpha,pi_l)
    
    if (new_like < likelihood){
      # return(list(params_tran_array,alpha))
      print("Tran Like Dec")
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

#above parameters are for debugging
#calculates likelihood of emission dist
#used in direct optimization
EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,vcovar_mat,vcovar_ind,weights_mat){
  
  #lower should theoretically be -Inf, but had some divergence issues
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

#optimizes activity mean 
#all emission dist param calculated this way
#comment out which parameter currently being optimized
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
  
  mu_act <- optimize(EmitLogLike, c(-20,10), act = act, light = light, 
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

#Highest level for optimizing emission dist
#takes function as input, easier to process this way
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

#inverts matrix, throws error if singular
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

#calculates overall likelihood by summing individual likelihood
CalcLikelihood <- function(alpha,pi_l){
  return(sum(CalcLikelihoodIndVec(alpha,pi_l)))
}

#calculates likelihood vector, each individual as entry
CalcLikelihoodIndVec <- function(alpha,pi_l){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,pi_l,i,num_obs)))
    like_vec <- c(like_vec,ind_like)
  }
  return(like_vec)
}

#adds pi to F-B output
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

#Calculates likelihood of individual
IndLike <- function(alpha,pi_l,ind,len){
  likelihood <- logSumExp(SumOverREIndTime(alpha,pi_l,ind,len))
  return(likelihood)
}

#calculates BIC
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

#Defunct, old way of calculating baseline haz pre additional surv covar
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

#calculates linear value of risk due to sociodemo covar, does not incl mixture
SurvCovarRiskVec <- function(surv_covar,surv_coef){
  surv_covar[[1]] <- matrix(surv_covar[[1]],ncol = 1)
  surv_coef[[1]] <- matrix(surv_coef[[1]],ncol = 1)
  surv_covar_risk_vec <- rowSums(mapply(function(x, y) x %*% y,  surv_covar,surv_coef, SIMPLIFY = T))
  return(surv_covar_risk_vec)
}

#calculates non-parametric baseline haz using breslow estimator
CalcBLHaz <- function(surv_coef,beta_vec, re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar){
  n <- length(surv_event)
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)
  
  for (time_ind in 1:n){
    risk_set <- surv_time >= surv_time[time_ind]
    #old defunct way
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

#Calculates likelihood only due to survival component
SurvLike <- function(beta_vec,surv_covar_risk_vec,surv_coef){
  #need to global assign for direct optimization
  bhaz_vec <<- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
  bline_vec <<- bhaz_vec[[1]]
  cbline_vec <<- bhaz_vec[[2]]
  
  
  loglike <- sum(log(bline_vec^surv_event) + 
                   ((re_prob %*% beta_vec)+surv_covar_risk_vec) * surv_event - 
                   cbline_vec*exp((re_prob %*% beta_vec)+surv_covar_risk_vec))
  
  return(-loglike)
}

#converts long vector of survival coef into beta and sociodemo list
OutofBetaSurvCoef <- function(beta_surv_coef,surv_coef_len){
  if (mix_num > 1){
    beta_vec <- c(0,beta_surv_coef[2:mix_num])
  } else {
    beta_vec <- 0
  }
    
  surv_coef_new <- list()
  surv_coef_new <- append(surv_coef_new,beta_surv_coef[[1]])
  
  
  surv_coef_len_alt <- cumsum(surv_coef_len-1)
  
  
  for (i in 1:(length(surv_coef_len)-1)) {
    coef_vec <- c(0,beta_surv_coef[(mix_num+1+surv_coef_len_alt[i]):(mix_num+surv_coef_len_alt[i+1])])
    surv_coef_new <- append(surv_coef_new,list(coef_vec))
  }
  
  return(list(beta_vec,surv_coef_new))
}

#input is beta vector and sociodemo coef list, output is long combined vector
#age first
IntoBetaSurvCoef <- function(beta_vec,surv_coef){
  beta_surv_coef_len <- mix_num + length(unlist(surv_coef)) - length(surv_coef)
  beta_surv_coef <- numeric(beta_surv_coef_len)
  
  beta_surv_coef[1] <- surv_coef[[1]]
  if(mix_num > 1){
    beta_surv_coef[2:mix_num] <- beta_vec[-1]
  }
    
  
  altered_surv_coef <- surv_coef[-1]
  for (i in 1:length(altered_surv_coef)){
    altered_surv_coef[[i]] <- altered_surv_coef[[i]][-1]
  }
  
  beta_surv_coef[(mix_num+1):beta_surv_coef_len] <- unlist(altered_surv_coef)
  
  return(beta_surv_coef)
}

#LM approach for calculating survival coefficients
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
      
      

      if (mix_num > 1){
        for (beta_ind in 2:mix_num){
          num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk))
          grad[beta_ind] <- grad[beta_ind] + (re_prob[ind,beta_ind] - num/denom)
          
          hess[beta_ind,beta_ind] <- hess[beta_ind,beta_ind] + (num/denom) - (num/denom)^2
          
          num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[1]][risk_set])
          hess[1,beta_ind] <- hess[1,beta_ind] + (denom*num_cross - num*num0)/denom^2
          hess[beta_ind,1] <- hess[beta_ind,1] + (denom*num_cross - num*num0)/denom^2
        }
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
        if (mix_num > 1){
          for (beta_ind in 2:mix_num){
            
            num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk))
            
            for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
              #sample code, vectorized this now
              # num_cross2 <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[2]][risk_set,2])
              # hess[5,beta_ind] <- hess[5,beta_ind] + (denom*num_cross2 - num*num_list[[2]][[2]])/denom^2
              
              num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1])
              hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] <- hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] + 
                (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2
              
              hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] <- hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] + 
                (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2
            }
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

    #LM aspect
    step_size <- .01
    step_mat <- matrix(0,dim(hess)[1],dim(hess)[2])
    nr_fact <- numeric(length(grad))-1
    max_val <- 6
    while(slike_diff < 0|all(nr_fact == -1) | max_val > 5){
      diag(step_mat) <- (diag(step_mat) + step_size - .01)
      nr_fact <- SolveCatch(hess + step_mat,-grad)
      step_size <- step_size * 10
    
      
      beta_surv_coef_new <- beta_surv_coef-nr_fact
      
      beta_vec_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len)[[1]]
      surv_coef_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len)[[2]]
      max_val <- abs(max(c(beta_vec_new,unlist(surv_coef_new))))
      
      
      new_slike <- SurvLike(beta_vec_new,surv_covar_risk_vec,surv_coef_new)
      slike_diff <- old_slike - new_slike
    }

    l2norm <- sum(sqrt((beta_surv_coef_new - beta_surv_coef)^2))



    beta_surv_coef <- beta_surv_coef_new

  }
  return(list(beta_surv_coef,sqrt(diag(solve(hess)))))
}

RemFirCol <- function(x){return(x[,-1])}

#Calculates beta, manual LM for JM, standard Cox for 2 stage
CalcBeta <- function(beta_surv_coef, combined_covar_mat,surv_covar_risk_vec ,incl_surv, stop_crit = .1){
  
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

#turns vector into dummy matrix
#useful for sociodemo covar
Vec2Mat <- function(vect){
  mat <- matrix(0,nrow = length(vect),ncol = max(vect))
  
  for(i in 1:length(vect)){
    mat[i,vect[i]] <- 1
  }
  return (mat)
}

#only used for 96 period len
#used in singleday
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

#determines week/weekend
FirstDay2WeekInd <- function(first_day){
  
  if (period_len == 96){
    weekday <- numeric(96)
    friday <- c(rep(0,68),rep(1,28))
    saturday <- numeric(96)+1
    sunday <- c(rep(1,68),rep(0,28))
  } else{
    weekday <- numeric(period_len)
    friday <- c(rep(0,period_len * 2 / 3),rep(1,period_len/3))
    saturday <- numeric(period_len)+1
    sunday <- c(rep(1,period_len * 2 / 3),rep(0,period_len/3))
  }
    
  
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

#turns transition matrices into dataframe for analyzing later
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
      
      tran_mat <- Params2TranVectorTresid(re_ind,period_len,params_tran)
      tosleep <- tran_mat[,3]
      towake <- tran_mat[,2]
      
      tran_df_working <- data.frame(prob = c(tosleep,towake),
                                    type = rep(c("Falling Asleep", "Waking"),each= period_len),
                                    time = rep(c(1:period_len)/4,2),
                                    weekend = vcovar_ind,
                                    mixture = re_ind)
      tran_df <- rbind(tran_df,tran_df_working)
    }
  }
  
  return(tran_df)
}

#Calculates prob of being in each mixture, just based on stationary and moderate activity
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

#Calculates ordinal logistic regression coef for pi
#Uses LM
CalcNu <- function(nu_mat,re_prob,nu_covar_mat){
  if(dim(re_prob)[2] == 1){return(nu_mat)}
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
  
  return(nu_mat_new)
}

#viterbi algorithm for global decoding
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

#Viterbi algorithm on an individual given their parameters and conditioning on mixture
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

#calculates non-tobit approach, similar to paper 2
#semi defunct
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

#calulates non-tobit mean
#these are all semi-defunct as now need to directly optimize to account for LoD
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


#forward algorithm but only uses decoded data, not activity/light
ForwardAlt <- function(post_decode_collapsed,init,tran_list,vcovar_mat){
  alpha_list <- list()
  
  for (ind in 1:dim(post_decode_collapsed)[2]){
    alpha_array <- array(NA,dim = c(dim(post_decode_collapsed)[1],2,mix_num))
    
    for (clust_i in 1:mix_num){
      alpha_array[,,clust_i] <- ForwardIndAltC(post_decode_collapsed[,ind],init[clust_i,],tran_list,clust_i-1,vcovar_mat[,ind])
      alpha_list[[ind]] <- alpha_array
    }
  }
  return(alpha_list)
  
}

#organizes forward output, mostly done in C
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


#organizes backward output, mostly done in C
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
      beta_array[,,clust_i] <- BackwardIndC(act[,ind], light[,ind], tran_list, emit_act_week, emit_light_week, 
                                            emit_act_weekend, emit_light_weekend,clust_i-1, lod_act, lod_light, corr_mat, 
                                            lintegral_mat, 
                                            vcovar_mat[,ind],lambda_act_mat, lambda_light_mat, tobit)
      beta_list[[ind]] <- beta_array
    }
  }
  return(beta_list)
  
}

#calculates survival function, mainy used for model selection in LOCV 
CalcS <- function(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec){
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

#needed to calculate prob of censoring for brier score 
Vec2StepPlot <- function(cens_dist,cens_dist_time,t_i){
  for (i in 1:(length(cens_dist_time)-1)){
    if ((t_i >= cens_dist_time[i]) & (t_i < cens_dist_time[i+1])){
      return(cens_dist[i])
    }
  }
  return(cens_dist[length(cens_dist_time)])
  
}

BrierScore <- function(bs_t,surv_event,surv_time,cens_dist,cens_dist_time,surv_mat_ind,event_time){
  bs_sum <- 0
  cens_dist[length(cens_dist)] <- cens_dist[length(cens_dist)-1]
  
  for (ind in 1:num_of_people){
    
    sprob <- Vec2StepPlot(surv_mat_ind[,ind],event_time,bs_t)
    
    num1 <- sprob^2 * (surv_time[ind] <= bs_t) * surv_event[ind]
    num2 <- (1-sprob)^2 * (surv_time[ind] > bs_t) 
    
    
    denom1 <- Vec2StepPlot(cens_dist,cens_dist_time,surv_time[ind])
    denom2 <-  Vec2StepPlot(cens_dist,cens_dist_time,bs_t)
    
    bs_sum <- bs_sum + (num1/denom1) + (num2/denom2)
    if (is.na((num1/denom1) + (num2/denom2))){
      #debugging NA
      # print(ind)
      # print(bs_t)
      # break
    }
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


CalcIBS <- function(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment,surv_covar_risk_vec){
  event_time <- unique(sort(surv_time))
  cbline_vec_new <- unique(sort(cbline_vec))
  # cbline_vec_new <- basehaz(cox1,centered = F)[,1]
  
  stime_vec <- c()
  for (stime in event_time){
    stime_vec <- c(stime_vec,which(surv_time == stime)[1])
  }
  
  cbline_vec_new <- cbline_vec[stime_vec]
  
  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec)
  
  
  
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

CalcIBSNew <- function(surv_time,surv_event,cbline_vec,beta_vec,re_prob,surv_covar_risk_vec) {
  km_fit <- survfit(Surv(surv_time, 1 - surv_event) ~ 1)
  cens_dist <- c(1,summary(km_fit)$surv)
  cens_dist_time <- c(0,summary(km_fit)$time)
  G <- stepfun(km_fit$time, c(1, km_fit$surv))
  
  event_time <- unique(sort(surv_time))
  cbline_vec_new <- unique(sort(cbline_vec))
  # cbline_vec_new <- basehaz(cox1,centered = F)[,1]
  
  stime_vec <- c()
  for (stime in event_time){
    stime_vec <- c(stime_vec,which(surv_time == stime)[1])
  }
  
  cbline_vec_new <- cbline_vec[stime_vec]
  
  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec)
  
  ibs_score <- integrate(vBrierScore,lower = 0,upper = max(event_time),
                         surv_event = surv_event,surv_time = surv_time,cens_dist = cens_dist,
                         cens_dist_time = cens_dist_time,surv_mat_ind = surv_mat_ind,event_time = event_time,
                         rel.tol=.05)
  return(ibs_score$value/max(event_time))
}

SurvCovar2Coef <- function(covar_mat,max_val = .1){
  return(seq(0,max_val,length.out = dim(covar_mat)[2]))
}

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

readCpp( "Scripting/cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )
################## EM Setup ################## 
###### True Settings ###### 

#Sets up simulation sizing
if (sim_size == 0){
  day_length <- 96 * 2
  num_of_people <- 1000
  missing_perc <- .2
} else if (sim_size== 1) {
  day_length <- 96 * 4
  num_of_people <- 4000
  missing_perc <- .2
} else if (sim_size== 2) {
  day_length <- 96 * 7
  num_of_people <- 7000
  missing_perc <- .2
}
  


if (misspecification){
  #load in "true" model to simulate under
  #but starts parameters at misspecified model
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
  
  lambda_act_mat_true <- to_save[[2]][[13]]
  lambda_light_mat_true <- to_save[[2]][[13]]
  
  #useful for simulations, dont actually need but makes it easier to run multiple at once
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
    #starting conditions for misspecified parameters
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
    
    ###semi defunct
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
  #if specified correctly
  #starting values
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
  params_tran_weekend_true[,1] <- seq(-3,-2.1,length.out = mix_num)
  params_tran_weekend_true[,2] <- seq(1,.7,length.out = mix_num)
  params_tran_weekend_true[,3] <- seq(.1,.9,length.out = mix_num)
  params_tran_weekend_true[,4] <- seq(-1.8,-2.4,length.out = mix_num)
  params_tran_weekend_true[,5] <- seq(-1.3,-.8,length.out = mix_num)
  params_tran_weekend_true[,6] <- seq(-.9,-1.1,length.out = mix_num)
  
  params_tran_array_dim <- c(mix_num,6,vcovar_num)
  params_tran_array_true <- array(NA,dim = params_tran_array_dim)
  params_tran_array_true[,,1] <- params_tran_week_true 
  params_tran_array_true[,,2] <- params_tran_weekend_true 
  
  emit_act_week_true <- array(NA, c(2,2,mix_num))
  emit_act_week_true[1,1,] <- seq(4,5,length.out = mix_num)
  emit_act_week_true[1,2,] <- seq(2,3,length.out = mix_num)
  emit_act_week_true[2,1,] <- seq(2,1,length.out = mix_num)
  emit_act_week_true[2,2,] <- seq(3,2,length.out = mix_num)
  
  emit_act_weekend_true <- array(NA, c(2,2,mix_num))
  emit_act_weekend_true[1,1,] <- seq(5,6,length.out = mix_num)
  emit_act_weekend_true[1,2,] <- seq(3,4,length.out = mix_num)
  emit_act_weekend_true[2,1,] <- seq(2,1,length.out = mix_num)
  emit_act_weekend_true[2,2,] <- seq(3,2,length.out = mix_num)
  
  emit_act_true <- array(NA, c(2,2,mix_num,2))
  emit_act_true[,,,1] <- emit_act_week_true
  emit_act_true[,,,2] <- emit_act_weekend_true  
  
  emit_light_week_true <- array(NA, c(2,2,mix_num))
  emit_light_week_true[1,1,] <- seq(-2,-1,length.out = mix_num)
  emit_light_week_true[1,2,] <- seq(8,6,length.out = mix_num)
  emit_light_week_true[2,1,] <- seq(-17,-19,length.out = mix_num)
  emit_light_week_true[2,2,] <- seq(13,15,length.out = mix_num)
  
  emit_light_weekend_true <- array(NA, c(2,2,mix_num))
  emit_light_weekend_true[1,1,] <- seq(-3,-2,length.out = mix_num)
  emit_light_weekend_true[1,2,] <- seq(9,7,length.out = mix_num)
  emit_light_weekend_true[2,1,] <- seq(-18,-21,length.out = mix_num)
  emit_light_weekend_true[2,2,] <- seq(15,18,length.out = mix_num)
  
  emit_light_true <- array(NA, c(2,2,mix_num,2))
  emit_light_true[,,,1] <- emit_light_week_true
  emit_light_true[,,,2] <- emit_light_weekend_true
  
  corr_mat_true <- array(NA, c(mix_num,2,2))
  corr_mat_true[,1,1] <- seq(0,.2,length.out = mix_num)
  corr_mat_true[,2,1] <- seq(.1,.3,length.out = mix_num)
  corr_mat_true[,1,2] <- seq(.2,.4,length.out = mix_num)
  corr_mat_true[,2,2] <- seq(.3,.5,length.out = mix_num)
  
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

#loads data in for a hot start
if (load_data){
  model_name_loadin <- "JMHMM"
  folder_name <- paste(mix_num)
  if (incl_surv==0){
    model_name_loadin <- paste0(model_name_loadin,"NoSurv")
    folder_name <- paste0("NS",folder_name)
  }
  
  model_name_loadin <- paste0(model_name_loadin,"Mix",mix_num,"Seed",".rda")
  
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
  
  if (length(to_save[[2]]) > 8){
    re_prob_true <- to_save[[2]][[10]]
    lambda_act_mat_true <- to_save[[2]][[13]]
    lambda_light_mat_true <- to_save[[2]][[13]]
    re_prob <- re_prob_true
  } 
  
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
  surv_list <- simulated_hmm[[8]]
  surv_covar_sim <- simulated_hmm[[9]]
  
  id_sim <- cbind(age_vec,surv_covar_sim-1)
  surv_covar <- list(age_vec,Vec2Mat(surv_covar_sim))
  surv_coef <- list(beta_age_true,beta_covar_sim)
  surv_coef_true <- surv_coef
  combined_covar_mat <- matrix(surv_covar_sim-1,nrow = num_of_people)
  combined_covar_mat <- as.factor(combined_covar_mat)
  
  surv_time <- surv_list[[1]]
  surv_event <- surv_list[[2]]
  
  #in simulated data sample weights are set to 0
  log_sweights_vec <- numeric(dim(act)[2])
  
  
} 
###### Read in Data ###### 

if (real_data) {
  #loads in NHANES
  setwd("Data/")
  load("NHANES_2011_2012_2013_2014.rda")
  nhanes1 <- NHANES_mort_list[[1]] %>% filter(eligstat == 1)
  nhanes2 <- NHANES_mort_list[[2]] %>% filter(eligstat == 1)
  lmf_data <- rbind(nhanes1,nhanes2)
  
  #depending on period len, loads in different data
  if (period_len == 24){
    load("Wavedata24_G.rda")
    load("Wavedata24_H.rda")
  } else if(period_len == 96){
    load("Wavedata_G.rda")
    load("Wavedata_H.rda")
  } else if(period_len == 1440){
    load("Wavedata1440_G.rda")
    load("Wavedata1440_H.rda")
  }
  
  setwd("..")
  
  #preps and combines 2 wakes of act data
  act_G <- wave_data_G[[1]]
  act_H <- wave_data_H[[1]]
  act <- rbind(act_G,act_H)
  act <- t(act[,-1])
  act0 <- act == 0
  act <- log(act)
  lod_act <- min(act[act!=-Inf],na.rm = T) - 1e-5
  act[act0] <- lod_act
  
  #preps and combines 2 wakes of light data
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
  
  #matches actigraphy data to public mortality data
  seqn_com_id <- id$SEQN %in% lmf_data$seqn
  seqn_com_lmf <- lmf_data$seqn %in% id$SEQN
  
  id <- id[seqn_com_id,]
  act <- act[,seqn_com_id]
  light <- light[,seqn_com_id]
  
  lmf_data <- lmf_data[seqn_com_lmf,]
  
  #sanity check
  if (sum(id$SEQN - lmf_data$seqn) != 0){print("LMF NOT LINKED CORRECTLY")}
  
  #sample weights for 2 waves
  log_sweights_vec <- log(id$sweights/2)
  
  id <- id %>% mutate(age_disc = case_when(age <=30 ~ 1,
                                           age <=50 & age > 30 ~ 2,
                                           age <=65 & age > 50 ~ 3,
                                           age > 65 ~ 4))
  
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  
  id$modact <- id$modact - 1
  
  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm
  
  #subset to smaller portion of data
  if(subset_data){
    subset_ind <- which(id$age >= 70)
    act <- act[,subset_ind]
    light <- light[,subset_ind]
    id <- id[subset_ind,]
    surv_event <- surv_event[subset_ind]
    surv_time <- surv_time[subset_ind]
  }
  
  #resample data for variance estimation
  if (bootstrap){
    boot_inds <- sample(dim(act)[2],dim(act)[2],T)
    act <- act[,boot_inds]
    light <- light[,boot_inds]
    id <- id[boot_inds,]
    surv_event <- surv_event[boot_inds]
    surv_time <- surv_time[boot_inds]
    
  }
  
  #saves original data before LOCV
  act_old <- act
  light_old <- light
  
  #cross validation
  if (leave_out){
    
    setwd("Data")
    #previously calculated who is left in/out for each seed
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
  
  #if single day can reduce memory used
  if (single_day != 0){
    single_day_mat <- sapply(first_day_vec,FirstDay2SingleDay,target_day = single_day)
    new_act <- matrix(NA,period_len,dim(act)[2])
    new_light <- matrix(NA,period_len,dim(light)[2])
    vcovar_mat <- matrix(0,period_len,dim(light)[2])
    
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
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)

  #sets of sociodemo covar list
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

  #initializes survival covariates
  surv_coef_true <- lapply(surv_covar[-1],SurvCovar2Coef)
  surv_coef_true <- append(list(.05),surv_coef_true)
  surv_coef_len <- unlist(lapply(surv_coef_true,length))
  surv_coef <- surv_coef_true
  
  #all sociodemo covar values
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
#if doing cv, load in full data values for hot start
if (leave_out){
  model_name_loadin <- "JMHMM"
  if (incl_surv==0){model_name_loadin <- paste0(model_name_loadin,"NoSurv")}
  if (incl_surv==1){model_name_loadin <- paste0(model_name_loadin,"HalfSurv")}
  model_name_loadin <- paste0(model_name_loadin,"Mix",mix_num,"Seed",".rda")
  
  print(paste("Loading",model_name_loadin))
  setwd("Data")
  load(model_name_loadin)
  setwd("..")
  
  #used if loading non-standard data
  # if (incl_surv == 2){
  #   model_name_loadin <- paste0("JMHMMMix",mix_num,"Seed",sim_num,".rda")
  #   foldername <- paste0(mix_num)
  # } else {
  #   model_name_loadin <- paste0("JMHMMNoSurvMix",mix_num,"Seed",sim_num,".rda")
  #   foldername <- paste0("NS",mix_num)
  # }
  # print(paste("Loading",model_name_loadin))
  # setwd(foldername)
  # load(model_name_loadin)
  # setwd("..")
  
  
  init_true <- to_save[[2]][[1]]
  params_tran_array_true <- to_save[[2]][[2]]
  emit_act_true <- to_save[[2]][[3]]
  emit_light_true <- to_save[[2]][[4]]
  corr_mat_true <- to_save[[2]][[5]]
  nu_mat_true <- to_save[[2]][[6]]
  pi_l_true <- CalcPi(nu_mat_true,nu_covar_mat)
  beta_vec_true <- to_save[[2]][[7]]
  surv_coef_true <- to_save[[2]][[8]]
  re_prob_true <- to_save[[2]][[10]]
  
  mix_assignment_true <- apply(to_save[[2]][[10]],1,which.max)
  re_prob <- to_save[[2]][[10]][-leave_out_inds,]
  
  lambda_act_mat <- to_save[[2]][[13]]
  lambda_light_mat <- to_save[[2]][[13]]
  
}

################## EM ##################

##### randomize starting parameters
# init <- matrix(rep(.5,mix_num*2),ncol = 2)
init <- init_true

params_tran_array <- params_tran_array_true + runif(unlist(length(params_tran_array_true)),-randomize_init*2,randomize_init*2)


emit_act <- emit_act_true + runif(length(unlist(emit_act_true)),-randomize_init,randomize_init)
emit_act[,2,,] <- abs(emit_act[,2,,])

emit_light <- emit_light_true + runif(length(unlist(emit_light_true)),-randomize_init*2,randomize_init*2)
emit_light[,2,,] <- abs(emit_light[,2,,])

#makes sure correlation makes sense
corr_mat <- corr_mat_true + runif(length(unlist(corr_mat_true)),-randomize_init/5,randomize_init/5)
corr_mat[corr_mat>.99] <- .99
corr_mat[corr_mat<-.99] <- -.99

#makes sure first val is always reference
if(misspecification == 1){beta_vec <- numeric(length(beta_vec))}
beta_vec <- beta_vec_true + runif(mix_num,-randomize_init,randomize_init)
beta_vec[1] <- 0

for (i in 1:length(surv_coef)){
  surv_coef[[i]] <-surv_coef_true[[i]]  +  runif(length(surv_coef_true[[i]]),-randomize_init/10,randomize_init/10)
  if (length(surv_coef[[i]]) != 1){surv_coef[[i]][1] <- 0} 
}
surv_coef[[1]] <- surv_coef_true[[1]] + runif(1,-randomize_init/100,randomize_init/100)

#dont randomize as these are very sensitive
nu_mat <- nu_mat_true
lambda_act_mat <- lambda_act_mat_true
lambda_light_mat <- lambda_light_mat_true

time_vec <- c()
pi_l <- CalcPi(nu_mat,nu_covar_mat)

#sets some controls so matrix sizing lines up
if (!leave_out & !load_data & misspecification == 0){re_prob <- pi_l}
if (subset_data & load_data){re_prob <- pi_l}
if (is.null(dim(re_prob))){re_prob <- matrix(re_prob,ncol = 1)}

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

#calculates baseline hazards
bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

#caluclates case4 probabilities ahead of time and transition list
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
if (incl_surv == 2 & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef)}
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood
#check to make sure all values are the same, simple sanity check
# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
iter_count <- 1
stop_crit <- 1e-2
if(mix_num > 8){stop_crit <- stop_crit * 10}
if(mix_num > 12){stop_crit <- stop_crit * 10}
if(mix_num > 15){stop_crit <- stop_crit * 5}

while((abs(like_diff) > stop_crit | iter_count < 5) & !run_only_surv){
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####
  
  #### Mixing Proportion  #####
  re_prob <- CalcProbRE(alpha,pi_l)
  
  ##### Survival ####
  #need model to fit a bit first otherwise may run into some instability
  if(beta_bool){

    nu_mat  <- CalcNu(nu_mat,re_prob,nu_covar_mat)
    pi_l <- CalcPi(nu_mat,nu_covar_mat)
    re_prob <- CalcProbRE(alpha,pi_l)
    
    if (incl_surv == 2){
      #calculates survival coef for JM 
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
  
  #### Weights  #####
  #calculates wake/sleep probabilities, needed for emission dist estimation
  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array_wake <- exp(weights_array_list[[1]])
  weights_array_sleep <- exp(weights_array_list[[2]])

  
  ##### Bivariate Normal Est  #####
  ##### Mixture Normal Param
  
  if (tobit){
    #calculates bivariate normal parameters wrt to LoD
    if(incl_light){

      emit_light[1,1,,] <- UpdateNorm(CalcLightMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_light[2,1,,] <- UpdateNorm(CalcLightMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

      emit_light[1,2,,] <- UpdateNorm(CalcLightSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_light[2,2,,] <- UpdateNorm(CalcLightSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    }

    if (incl_act){
      emit_act[1,2,,] <- UpdateNorm(CalcActSig,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_act[2,2,,] <- UpdateNorm(CalcActSig,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

      emit_act[1,1,,] <- UpdateNorm(CalcActMean,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      emit_act[2,1,,] <- UpdateNorm(CalcActMean,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
    }

    if (incl_act & incl_light){
      corr_mat[,1,] <- UpdateNorm(CalcBivarCorr,1,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)
      corr_mat[,2,] <- UpdateNorm(CalcBivarCorr,2,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array_wake,weights_array_sleep,vcovar_mat+1)

    }

  } else {
    #semi defunct
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
  #this only relies on normal parameters so calculate it now
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  
  #saves old transition values in case likelihood decrease
  params_tran_array_old <- params_tran_array
  #gradient and hessian for tran parameters
  tran_gradhess_list <- CalcTranCHelper(alpha,beta,act,light,params_tran_array,
                                        emit_act,emit_light,corr_mat,
                                        pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat,
                                        lambda_act_mat, lambda_light_mat, tobit, check_tran,likelihood)
  
  params_tran_array <- LM(tran_gradhess_list[[1]],tran_gradhess_list[[2]],params_tran_array,check_tran,likelihood,pi_l)
  
  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

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
  
  #if JM but during cold start, dont wan to actually include survial in likelihood yet
  if (incl_surv == 2 & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef)}
  
  like_diff <- new_likelihood - likelihood
  
  if (like_diff < 0){
    #all other parameters are either
    #1) closed form
    #2) we can quickly calculate likelihood difference
    #tran likelihood requires forward algorithm and thus much slower
    #thus any like decrease is from transition parameters
    #effectively just doesnt optimize tran in this step of EM
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
    if (incl_surv == 2 & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef)}
    like_diff <- new_likelihood - likelihood
  }
  
  print(paste("RE num:",mix_num,"Like:",round(like_diff,6)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  # after a few iterations don't have any issues with survival estimation stability
  if (iter_count == 3){
    beta_bool <- T
    print("Starting to Est Survival and Age Mixing Effect")
  }
  
  iter_count <- iter_count + 1
  
  #saves info every few iterations just in case
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
                       beta_vec,surv_coef,
                       0,
                       re_prob)
                       # tran_df,
                       # re_prob,
                       # new_likelihood,
                       #lambda_act_mat,lambda_light_mat,
                       # beta_se)
    
    # bic <- CalcBIC(new_likelihood,mix_num,act,light)
    to_save <- list(true_params,est_params)#,bic)
    
    if(!leave_out){
      save(to_save,file = paste0("Inter",model_name))
    }
    
      
  }
  
  
  #reorders clusters from best to worst survival
  if ((abs(like_diff) < stop_crit*10) & !relabel_reset & !bootstrap){
    relabel_reset <- TRUE
    relabel_bool <- 0
    print("Relabelling")
    print("Potential Soft Reset")
    #### Reorder #####
    #Reorder to avoid label switching
    #Cluster means go from small to large by activity
    
    if (incl_surv == 2){
      reord_inds <- order(beta_vec)
    } else if (incl_surv == 0){
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,surv_covar_risk_vec,incl_surv)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],surv_coef_len)
      beta_vec_temp <- beta_surv_coef_temp_list[[1]]
      reord_inds <- order(beta_vec_temp)
    }
    
    # reord_inds <- c(0,rev(order(beta_vec[-1])))+1
    if (!all(reord_inds == c(1:mix_num)) & !leave_out){
      print("Swapping Labels")
      relabel_bool <- 1
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
      re_prob <- re_prob[,reord_inds]

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]

    }
    
    for (re_ind in 1:mix_num){
      for (week_ind in 1:2){
        if (emit_act[2,1,re_ind,week_ind] > emit_act[1,1,re_ind,week_ind]){
          relabel_bool <- 1
          print(paste("Swapping wake/sleep for week_ind",week_ind,"Mixture",re_ind))

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

          #ISSUE HERE
          temp <- params_tran_array[re_ind,1:3,week_ind]
          params_tran_array[re_ind,1:3,week_ind] <- params_tran_array[re_ind,4:6,week_ind]
          params_tran_array[re_ind,4:6,week_ind] <- temp
          tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num)

          temp <- corr_mat[re_ind,1,week_ind]
          corr_mat[re_ind,1,week_ind] <- corr_mat[re_ind,2,week_ind]
          corr_mat[re_ind,2,week_ind] <- temp
        }
      }
    }
    
    if (relabel_bool){
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
      like_diff <- stop_crit * 1.1
    }
      
  }
  
  
  #Finally calls beta at very end of while loop
  beta <- Backward(act = act,light = light, tran_list = tran_list,
                   emit_act = emit_act,emit_light = emit_light,
                   lod_act = lod_act, lod_light =  lod_light, 
                   corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                   lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
  
}

#if 2-stage model, calculate survival here
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

  
    



#if LOCV we don't care about viterbi decoding
if(!leave_out){
  decoded_mat <- Viterbi(act,light,vcovar_mat)
} else {
  decoded_mat <- matrix(NA,2,2)
}

#transition parameters for later analysis
tran_df <- ParamsArray2DF(params_tran_array)
if (!real_data & !misspecification){
  tran_df_true <- ParamsArray2DF(params_tran_array_true,misspecification)
  tran_df_truth <- tran_df_true[,1]
  
  tran_df <- tran_df %>% mutate(truth = tran_df_truth)
  tran_df <- tran_df %>% mutate(resid = prob - truth)
}

#concatenates true and est parameters into lists to save
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

#if doing leave 100 out cross validation
#predict cluster assignment using varying levels of information
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
    empty_list[[i]] <- list()
  }
  
  empty_mat_list <- vector(mode = "list", length = 9)
  for (i in 1:9){
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
    
    new_act_working[inclexcl_ind_mat == 0] <- NA
    new_light_working[inclexcl_ind_mat == 0] <- NA
    
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
    
    
    if (wake_sleep){
      
      beta <- Backward(act = new_act_working,light = new_light_working, tran_list = tran_list,
                       emit_act = emit_act,emit_light = emit_light,
                       lod_act = lod_act, lod_light =  lod_light, 
                       corr_mat = corr_mat,lintegral_mat = lintegral_mat,vcovar_mat = vcovar_mat,
                       lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
      
      weights_array_list <- CondMarginalize(alpha,beta,new_pi_l)
      weights_array_wake <- exp(weights_array_list[[1]])
      weights_array_wake_collapsed <- apply(weights_array_wake,c(1,2),sum)
      post_decode_collapsed <- weights_array_wake_collapsed < .5
      
      alpha <- ForwardAlt(post_decode_collapsed,init,tran_list,new_vcovar_mat)
    }
   
    
    re_prob_new <- CalcProbRE(alpha,new_pi_l)
    mix_assignment_pred <- apply(re_prob_new,1,which.max)
    
    mix_assignment_pred <- c(mix_assignment_pred,c(1:mix_num))
    mix_assignment_true_ind <- c(mix_assignment_true[leave_out_inds],c(1:mix_num))
    conf_mat_ind <- table(mix_assignment_pred,mix_assignment_true_ind)
    diag(conf_mat_ind) <- diag(conf_mat_ind) - 1
    
    cindex <- CalcCindex(surv_time_new,surv_event_new,beta_vec,surv_coef,re_prob_new,new_surv_covar,surv_covar_risk_vec_new)
    # ibs <- CalcIBS(surv_time_new,surv_event_new,cbline_vec,beta_vec,surv_coef,new_surv_covar,re_prob_new,incl_surv,mix_assignment_pred,surv_covar_risk_vec_new)
    ibs <- CalcIBSNew(surv_time_new,surv_event_new,cbline_vec,beta_vec,re_prob_new,surv_covar_risk_vec_new)

    
    re_prob_new_list[[days_incl]] <- append(re_prob_new_list[[days_incl]],list(re_prob_new))
    mix_assignment_pred_list[[days_incl]] <- append(mix_assignment_pred_list[[days_incl]],list(mix_assignment_pred))
    conf_mat_list[[days_incl]] <- conf_mat_list[[days_incl]] + conf_mat_ind
    cindex_new_list[[days_incl]] <- c(cindex_new_list[[days_incl]],cindex)
    ibs_new_list[[days_incl]] <- c(ibs_new_list[[days_incl]],ibs)

    
  }
  
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
  # ibs <- CalcIBS(surv_time_new,surv_event_new,cbline_vec,beta_vec,surv_coef,new_surv_covar,re_prob_new,incl_surv,mix_assignment_pred,surv_covar_risk_vec_new)
  ibs <- CalcIBSNew(surv_time_new,surv_event_new,cbline_vec,beta_vec,re_prob_new,surv_covar_risk_vec_new)
  
  # leave_out_weekend_list <- list(leave_out_inds,re_prob_new,mix_assignment_pred,conf_mat,cindex,ibs)
  leave_out_weekend_list <- list(leave_out_inds,conf_mat,cindex,ibs)
  leave_out_to_save <- list(leave_out_list,leave_out_weekend_list)
} else {
  leave_out_to_save <- list()
}

#if not using simulated data or want to save space
if (real_data | save_space){
  simulated_hmm <- list()
}


#mixture predections
mix_assignment <- apply(re_prob,1,which.max)
if (!real_data){
  tab <- table(c(mixture_mat,c(1:mix_num)),c(mix_assignment,c(1:mix_num)))
  diag(tab) <- diag(tab) - 1
} else {
  tab <- table(c(mix_assignment,c(1:mix_num)),c(mix_assignment,c(1:mix_num)))
  diag(tab) <- diag(tab) - 1
}
  
#removes some data from saving
if (save_space){
  # est_params[[10]] <- 0
  est_params[[12]] <- 0
  est_params[[15]] <- 0
  est_params[[16]] <- 0
}

#diagnostics
ibs <- CalcIBSNew(surv_time,surv_event,cbline_vec,beta_vec,re_prob,surv_covar_risk_vec)
cindex <- CalcCindex(surv_time,surv_event,beta_vec,surv_coef,re_prob,surv_covar,surv_covar_risk_vec)
diagnostics <- list(cindex,ibs,tab)
#save everything
bic <- CalcBIC(new_likelihood,mix_num,act,light)
to_save <- list(true_params,est_params,bic,leave_out_to_save,simulated_hmm,diagnostics)
setwd("/gpfs/gsfs12/users/aronjr/JM/Routputs")
#"~/JM/Routputs"
save(to_save,file = model_name)