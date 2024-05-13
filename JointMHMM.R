library(Rcpp)
library(RcppArmadillo)
library(matrixStats)

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

CalcLintegralMat <- function(){
  re_num <- dim(emit_act)[3]
  
  lintegral_mat <- matrix(NA,nrow = re_num, ncol = 2)
  for (i in 1:re_num){
    lintegral_mat[i,1] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                        emit_act[1,1,i],emit_act[1,2,i],emit_light[1,1,i],emit_light[1,2,i],corr_vec[1],light_LOD)[[1]])
    
    lintegral_mat[i,2] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                        emit_act[2,1,i],emit_act[2,2,i],emit_light[2,1,i],emit_light[2,2,i],corr_vec[1],light_LOD)[[1]])
  }
  
  return(lintegral_mat)
}

logClassification <- function(time,act,light,mu_act,sig_act,mu_light,sig_light,lod_act,light_LOD,bivar_corr,lintegral){

  #CASE 1
  if (act[time] > lod_act & light[time] > light_LOD){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      dnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 2 
  if (act[time] > lod_act & light[time] == light_LOD){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      pnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 3
  if (act[time] == lod_act & light[time] > light_LOD){
    
    mu_act_cond <- CalcCondMean(mu_act,sig_act,mu_light,sig_light,bivar_corr,light[time])
    sig_act_cond <- CalcCondSig(sig_act,bivar_corr)
    
    lognorm_dens <- pnorm(act[time],mu_act_cond,sig_act_cond,log = T) + 
      dnorm(light[time],mu_light,sig_light,log = T)
  }
  
  
  #CASE 4 
  if (act[time] == lod_act & light[time] == light_LOD){
    
    lognorm_dens <- log(integrate(Case4,lower = -Inf,upper = lod_act,
                                  mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD)[[1]])
  }

  return(lognorm_dens)
}


readCpp("cFunctions.cpp")


act <- cbind(rbinom(1000,2,.5),rbinom(1000,2,.5))
light <- cbind(rbinom(1000,2,.5),rbinom(1000,2,.5))


init <- c(.3,.7)

tran_covar_num <- 3
covar_mat_tran <- t(rmultinom(1000,1,rep(1/tran_covar_num,tran_covar_num)))
covar_mat_tran <- cbind((numeric(1000) + 1),covar_mat_tran[,2:tran_covar_num])

params_tran <- c(-3,.5,-.25,0,0,0,.9,-.8,-1,.1,.2,1.5,0,0,0,0,0,0,
                      -2.2,.3,-.3,0,0,0,-.75,.8,.75,.2,.3,.75,0,0,0,0,0,0)
tran_list <- lapply(c(1:dim(covar_mat_tran)[2]),TranByTimeVec, params_tran = params_tran, time_vec = c(1:96))


emit_act <- array(NA, c(2,2,2))
emit_act1 <- matrix(c(2,0,1,.5),2,2)
emit_act2 <- matrix(c(3,1,1,.5),2,2)
emit_act[,,1] <- emit_act1
emit_act[,,2] <- emit_act2

emit_light <- array(NA, c(2,2,2))
emit_light1 <- matrix(c(7,-1,3,2),2,2)
emit_light2 <- matrix(c(8,0,3,2),2,2)
emit_light[,,1] <- emit_light1
emit_light[,,2] <- emit_light2

corr_vec <- c(.6,.8)
clust_i <- 1
lod_act <- 0
light_LOD <- 0



lintegral_mat <- CalcLintegralMat()


x <- ForwardIndC(act_ind = act[,1],light_ind = light[,1],
            init = init,tran_list = tran_list,
            emit_act = emit_act,emit_light = emit_light,
            tran_ind = 3,clust_i = 0,lod_act = 0, lod_light =  0, 
            corr_vec = corr_vec,lintegral_mat = lintegral_mat,log_sweight = 0)

y <- BackwardIndC(act_ind = act[,1],light_ind = light[,1],tran_list = tran_list,
            emit_act = emit_act,emit_light = emit_light,
            tran_ind = 3,clust_i = 0,lod_act = 0, lod_light =  0, 
            corr_vec = corr_vec,lintegral_mat = lintegral_mat)

alpha <- ForwardC(act = act,light = light,
         init = init,tran_list = tran_list,
         emit_act = emit_act,emit_light = emit_light,
         tran_ind_vec = c(1,1),lod_act = 0, lod_light =  0, 
         corr_vec = corr_vec,lintegral_mat = lintegral_mat,log_sweight = c(0,0))

beta <- BackwardC(act = act,light = light, tran_list = tran_list,
              emit_act = emit_act,emit_light = emit_light,
              tran_ind_vec = c(1,1),lod_act = 0, lod_light =  0, 
              corr_vec = corr_vec,lintegral_mat = lintegral_mat)
         

apply(alpha[[2]][,,1]+beta[[2]][,,1],1,logSumExp)






