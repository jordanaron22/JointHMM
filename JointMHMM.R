library(Rcpp)
library(RcppArmadillo)

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

logClassification <- function(time,current_state,act,light,emit_act,emit_light,corr_vec,clust_i,act_LOD,light_LOD){
  
  mu_act <- emit_act[current_state+1,1,clust_i]
  sig_act <- emit_act[current_state+1,2,clust_i]
  
  mu_light <- emit_light[current_state+1,1,clust_i]
  sig_light <- emit_light[current_state+1,2,clust_i]
  
  bivar_corr <- corr_vec[current_state+1]

  #CASE 1
  if (act[time] > act_LOD & light[time] > light_LOD){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      dnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 2 
  if (act[time] > act_LOD & light[time] == light_LOD){
    
    mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act[time])
    sig_light_cond <- CalcCondSig(sig_light,bivar_corr)
    
    lognorm_dens <- dnorm(act[time],mu_act,sig_act,log = T) + 
      pnorm(light[time],mu_light_cond,sig_light_cond,log = T)
  }
  
  #CASE 3
  if (act[time] == act_LOD & light[time] > light_LOD){
    
    mu_act_cond <- CalcCondMean(mu_act,sig_act,mu_light,sig_light,bivar_corr,light[time])
    sig_act_cond <- CalcCondSig(sig_act,bivar_corr)
    
    lognorm_dens <- pnorm(act[time],mu_act_cond,sig_act_cond,log = T) + 
      dnorm(light[time],mu_light,sig_light,log = T)
  }
  
  
  #CASE 4 
  if (act[time] == act_LOD & light[time] == light_LOD){
    
    lognorm_dens <- log(integrate(Case4,lower = -Inf,upper = act_LOD,
                                  mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD)[[1]])
  }

  return(lognorm_dens)
}


readCpp("cFunctions.cpp")


act <- c(1,2,0,0)
light <- c(5,0,10,0)


act <- rbinom(1000,2,.5)
light <- rbinom(1000,2,.5)


emit_act <- matrix(c(2,0,1,.5),2,2)
emit_light <- matrix(c(8,0,3,2),2,2)

dim(emit_act) <- c(2,2,1)
dim(emit_light) <- c(2,2,1)

corr_vec <- c(.6,.8)
clust_i <- 1
act_LOD <- 0
light_LOD <- 0




lintegral <- log(integrate(Case4,lower = -Inf,upper = act_LOD,
                           emit_act[1,1,1],emit_act[1,2,1],emit_light[1,1,1],emit_light[1,2,1],corr_vec[1],light_LOD)[[1]])

x <- logClassificationCtest(act,light,
                            emit_act[1,1,1],emit_act[1,2,1],emit_light[1,1,1],emit_light[1,2,1],
                            0,0,corr_vec[1],lintegral)
x






