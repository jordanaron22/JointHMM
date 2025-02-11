//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
double CalcCondMeanC(double mu1, double sig1, double mu2, double sig2, double bivar_corr, double obs2){
	double cond_mu1 = mu1 + bivar_corr*(sig1/sig2)*(obs2-mu2);
	return (cond_mu1);
}

// [[Rcpp::export]]
double CalcCondSigC(double sig1, double bivar_corr){
	double cond_sig1 = sig1*sqrt(1-pow(bivar_corr,2));
	return cond_sig1;
}

// [[Rcpp::export]]
NumericVector vectorEqBool(NumericVector vec, double lod) {
  NumericVector vecBool;
  for(NumericVector::iterator i = vec.begin(); i != vec.end(); ++i) {
    if (*i == lod) {
      vecBool.push_back(1);
    } else {
      vecBool.push_back(0);
    }
  }
  return vecBool;
}

// [[Rcpp::export]]
vec logClassificationCnonTobit(NumericVector act_obs, NumericVector light_obs, double mu_act, double sig_act, double mu_light, double sig_light, 
                               double lod_act, double lod_light, double lambda_act, double lambda_light) {
  vec temp_act;
  vec temp_light;
  
  NumericVector vec_eq_act = vectorEqBool(act_obs, lod_act);
  NumericVector vec_eq_light = vectorEqBool(light_obs, lod_light);
  
  temp_act = ((log(1-lambda_act) + Rcpp::dnorm( act_obs, mu_act, sig_act, true )) * (1 - vec_eq_act)) + (log(lambda_act)*vec_eq_act);	
  temp_light = ((log(1-lambda_light) + Rcpp::dnorm( light_obs, mu_light, sig_light, true )) * (1 - vec_eq_light)) + (log(lambda_light)*vec_eq_light);	
  
  LogicalVector nan_vec_act = is_na(act_obs);
  LogicalVector nan_vec_light = is_na(light_obs);
  int n = act_obs.length();
  
  for (int i = 0; i < n; i++) {
    if(nan_vec_act[i] == true){
      temp_act[i] = 0;
    }
    
    if(nan_vec_light[i] == true){
      temp_light[i] = 0;
    }
  }
  
  return temp_act + temp_light;
}


// [[Rcpp::export]]
vec logClassificationCTobit(NumericVector act, NumericVector light, double mu_act, double sig_act, double mu_light, double sig_light, 
                            double lod_act, double lod_light, double bivar_corr, double lintegral) {

	int n = act.length();
	vec log_dens_vec(n);

	LogicalVector nan_vec_act = is_na(act);
	LogicalVector nan_vec_light = is_na(light);
	
	for (int i = 0; i < n; i++) {
		double act_obs = act(i);
		double light_obs = light(i);

		// Observe both act and light
		if (nan_vec_act(i) != true && nan_vec_light(i) != true){
			// CASE 1
			if(act_obs > lod_act & light_obs > lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::dnorm( light_obs, mu_light_cond, sig_light_cond, true ) ;
			}

			//CASE 2
			if(act_obs > lod_act & light_obs <= lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::pnorm( light_obs, mu_light_cond, sig_light_cond, true, true ) ;
			}

			//CASE 3
			if(act_obs <= lod_act & light_obs > lod_light){

				double mu_act_cond = CalcCondMeanC(mu_act,sig_act,mu_light,sig_light,bivar_corr,light_obs);
				double sig_act_cond = CalcCondSigC(sig_act,bivar_corr);

				log_dens_vec(i) = R::pnorm( act_obs, mu_act_cond, sig_act_cond, true, true ) + R::dnorm( light_obs, mu_light, sig_light, true ) ;
			}

			//CASE 4
			if(act_obs <= lod_act & light_obs <= lod_light){
				log_dens_vec(i) = lintegral;
			}

			// Light Missing
		} else if (nan_vec_act(i) != true && nan_vec_light(i) == true){
			if(act_obs > lod_act){	
				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true );
			} else {
				log_dens_vec(i) = R::pnorm( act_obs, mu_act, sig_act, true, true );
			}
			
			// Act missing
		} else if((nan_vec_act(i) == true && nan_vec_light(i) != true)){
			if(light_obs > lod_light){	
				log_dens_vec(i) = R::dnorm( light_obs, mu_light, sig_light, true );
			} else {
				log_dens_vec(i) = R::pnorm( light_obs, mu_light, sig_light, true, true );
			}

			// Both Missing
		} else {
			log_dens_vec(i) = 0;
		}

		
	}

return log_dens_vec;

}

// [[Rcpp::export]]
vec logClassificationC(NumericVector act, NumericVector light, double mu_act, double sig_act, double mu_light, double sig_light, 
                       double lod_act, double lod_light, double bivar_corr, double lintegral, double lambda_act, double lambda_light,
                       bool tobit) {
  
  vec log_class;
  
  if (tobit){
    log_class = logClassificationCTobit(act, light, mu_act, sig_act, mu_light, sig_light, 
                                        lod_act, lod_light, bivar_corr, lintegral);
  } else {
    log_class = logClassificationCnonTobit(act, light, mu_act, sig_act, mu_light, sig_light, 
                                           lod_act, lod_light, lambda_act, lambda_light);
  }
  return log_class;
}




/* This is from the seqHMM github*/
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#  define EXPL expl
#else
#  define LDOUBLE double
#  define EXPL exp
#endif

// [[Rcpp::export]]
double logSumExpC(const arma::vec& x) {
  unsigned int maxi = x.index_max();
  LDOUBLE maxv = x(maxi);
  if (!(maxv > -arma::datum::inf)) {
    return -arma::datum::inf;
  }
  LDOUBLE cumsum = 0.0;
  for (unsigned int i = 0; i < x.n_elem; i++) {
    if ((i != maxi) && (x(i) > -arma::datum::inf)) {
      cumsum += EXPL(x(i) - maxv);
    }
  }
  
  return maxv + log1p(cumsum);
}

// [[Rcpp::export]]
mat ForwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, NumericVector init, Rcpp::List tran_list, 
                cube emit_act_week, cube emit_light_week,  cube emit_act_weekend, cube emit_light_weekend, 
                int clust_i, double lod_act, double lod_light, cube corr_vec, vec beta_vec, double covar_risk, int event, double bline, double cbline, cube lintegral_mat, double log_sweight, vec vcovar_vec,
                cube lambda_act_mat, cube lambda_light_mat, bool tobit, bool incl_surv){

	mat alpha( act_ind.length(), 2 );
  
  vec log_class_0_week = logClassificationC( act_ind, light_ind, 
                                             emit_act_week(0,0,clust_i), 
                                             emit_act_week(0,1,clust_i), 
                                             emit_light_week(0,0,clust_i),
                                             emit_light_week(0,1,clust_i), 
                                             lod_act, lod_light, corr_vec(clust_i,0,0), lintegral_mat(clust_i,0,0),
                                             lambda_act_mat(clust_i,0,0),lambda_light_mat(clust_i,0,0),tobit);
  
  vec log_class_1_week = logClassificationC( act_ind, light_ind, 
                                             emit_act_week(1,0,clust_i),
                                             emit_act_week(1,1,clust_i),
                                             emit_light_week(1,0,clust_i),
                                             emit_light_week(1,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,1,0), lintegral_mat(clust_i,1,0),
                                             lambda_act_mat(clust_i,1,0),lambda_light_mat(clust_i,1,0),tobit);
  
  vec log_class_0_weekend = logClassificationC( act_ind, light_ind, 
                                                emit_act_weekend(0,0,clust_i), 
                                                emit_act_weekend(0,1,clust_i), 
                                                emit_light_weekend(0,0,clust_i),
                                                emit_light_weekend(0,1,clust_i), 
                                                lod_act, lod_light, corr_vec(clust_i,0,1), lintegral_mat(clust_i,0,1),
                                                lambda_act_mat(clust_i,0,1),lambda_light_mat(clust_i,0,1),tobit);
  
  vec log_class_1_weekend = logClassificationC( act_ind, light_ind, 
                                                emit_act_weekend(1,0,clust_i),
                                                emit_act_weekend(1,1,clust_i),
                                                emit_light_weekend(1,0,clust_i),
                                                emit_light_weekend(1,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,1,1), lintegral_mat(clust_i,1,1),
                                                lambda_act_mat(clust_i,1,1),lambda_light_mat(clust_i,1,1),tobit);
	
	
	vec log_class_0 = (log_class_0_week % (1-vcovar_vec)) + (log_class_0_weekend % vcovar_vec); 
	vec log_class_1 = (log_class_1_week % (1-vcovar_vec)) + (log_class_1_weekend % vcovar_vec);

	double surv_comp = 0.0;

	if (incl_surv){
	  if (event == 1){
	    surv_comp = log(bline) + beta_vec[clust_i]+covar_risk - (cbline * exp(beta_vec[clust_i]+(covar_risk)));
	  } else {
	    surv_comp =  -cbline * exp(beta_vec[clust_i]+covar_risk);
	  }
	}
  	
	
	

	alpha(0,0) = log(init(0)) + log_class_0[0] + surv_comp;
	alpha(0,1) = log(init(1)) + log_class_1[0] + surv_comp;
	
	List tran_list_clust = tran_list[clust_i]; 
	
	for (int i = 1; i < act_ind.length(); i++) {
	  List tran_list_vcovar = tran_list_clust[vcovar_vec(i)];
	  
	  NumericMatrix tran = tran_list_vcovar[i%96];
	  
	  double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0[i];
	  
	  double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0[i];
	  
	  double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1[i];
	  
	  double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1[i];
	  
	  NumericVector fp_0 = NumericVector::create(fp_00,fp_10);
	  NumericVector fp_1 = NumericVector::create(fp_01,fp_11);
	  
	  alpha(i,0) = logSumExpC(fp_0);
	  alpha(i,1) = logSumExpC(fp_1);
	}

	alpha = alpha + log_sweight;

	return(alpha);
}

// [[Rcpp::export]]
mat BackwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, Rcpp::List tran_list, 
                 cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend,
                 int clust_i, double lod_act, double lod_light, cube corr_vec, cube lintegral_mat, vec vcovar_vec,
                 cube lambda_act_mat, cube lambda_light_mat, bool tobit){
  
  int n = act_ind.length(); 
  mat beta( n, 2 );
  
  vec log_class_0_week = logClassificationC( act_ind, light_ind, 
                                             emit_act_week(0,0,clust_i), 
                                             emit_act_week(0,1,clust_i), 
                                             emit_light_week(0,0,clust_i),
                                             emit_light_week(0,1,clust_i), 
                                             lod_act, lod_light, corr_vec(clust_i,0,0), lintegral_mat(clust_i,0,0),
                                             lambda_act_mat(clust_i,0,0),lambda_light_mat(clust_i,0,0),tobit);
  
  vec log_class_1_week = logClassificationC( act_ind, light_ind, 
                                             emit_act_week(1,0,clust_i),
                                             emit_act_week(1,1,clust_i),
                                             emit_light_week(1,0,clust_i),
                                             emit_light_week(1,1,clust_i),
                                             lod_act, lod_light, corr_vec(clust_i,1,0), lintegral_mat(clust_i,1,0),
                                             lambda_act_mat(clust_i,1,0),lambda_light_mat(clust_i,1,0),tobit);
  
  vec log_class_0_weekend = logClassificationC( act_ind, light_ind, 
                                                emit_act_weekend(0,0,clust_i), 
                                                emit_act_weekend(0,1,clust_i), 
                                                emit_light_weekend(0,0,clust_i),
                                                emit_light_weekend(0,1,clust_i), 
                                                lod_act, lod_light, corr_vec(clust_i,0,1), lintegral_mat(clust_i,0,1),
                                                lambda_act_mat(clust_i,0,1),lambda_light_mat(clust_i,0,1),tobit);
  
  vec log_class_1_weekend = logClassificationC( act_ind, light_ind, 
                                                emit_act_weekend(1,0,clust_i),
                                                emit_act_weekend(1,1,clust_i),
                                                emit_light_weekend(1,0,clust_i),
                                                emit_light_weekend(1,1,clust_i),
                                                lod_act, lod_light, corr_vec(clust_i,1,1), lintegral_mat(clust_i,1,1),
                                                lambda_act_mat(clust_i,1,1),lambda_light_mat(clust_i,1,1),tobit);
	
	
	vec log_class_0 = (log_class_0_week % (1-vcovar_vec)) + (log_class_0_weekend % vcovar_vec); 
	vec log_class_1 = (log_class_1_week % (1-vcovar_vec)) + (log_class_1_weekend % vcovar_vec); 
	  
  beta(n-1,0) = log(1);
  beta(n-1,1) = log(1);
  
  List tran_time_list = tran_list[clust_i]; 
  
  for (int i = n-2; i >= 0; i--) {
    
    List tran_list_vcovar = tran_time_list[vcovar_vec(i+1)];
    NumericMatrix tran = tran_list_vcovar[(i+1)%96];
    
    double bp_00 = log(tran(0,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_01 = log(tran(0,1)) + log_class_1[i+1] + beta(i+1,1);
    
    double bp_10 = log(tran(1,0)) + log_class_0[i+1] + beta(i+1,0);
    
    double bp_11 = log(tran(1,1)) + log_class_1[i+1] + beta(i+1,1);
    
    NumericVector bp_0 = NumericVector::create(bp_00,bp_01);
    NumericVector bp_1 = NumericVector::create(bp_10,bp_11);
    
    beta(i,0) = logSumExpC(bp_0);
    beta(i,1) = logSumExpC(bp_1);

  }
  
  return beta;

}

/*
// [[Rcpp::export]]
List ForwardC(const NumericMatrix& act, const NumericMatrix& light, NumericMatrix init, List tran_list, 
              cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend, 
              double lod_act, double lod_light, cube corr_mat, vec beta_vec, double beta_age, 
              vec event_vec, vec bline_vec, vec cbline_vec, cube lintegral_mat, NumericVector log_sweights_vec, 
              vec age_vec, mat vcovar_mat,cube lambda_act_mat, cube lambda_light_mat, bool tobit, bool incl_surv){ 
  
	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act_week.n_slices;
	List alpha_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
	  int age_ind = age_vec(ind);
	  
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);
		double log_sweight = log_sweights_vec(ind);
		
		vec vcovar_vec = vcovar_mat.col(ind);
		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			NumericVector init_vec = init.row(clust_i);

			Cube1.slice(clust_i) = ForwardIndC(act_ind, light_ind, init_vec, tran_list, emit_act_week, emit_light_week, 
               emit_act_weekend, emit_light_weekend,clust_i, lod_act, lod_light, corr_mat, 
               beta_vec, beta_age,age_ind,event_vec[ind], bline_vec[ind], cbline_vec[ind], lintegral_mat, 
               log_sweight, vcovar_vec,lambda_act_mat, lambda_light_mat, tobit, incl_surv);

		}

		alpha_list(ind) = Cube1;
	}
	return(alpha_list);
}
*/

// [[Rcpp::export]]
List BackwardC(const NumericMatrix& act, const NumericMatrix& light, List tran_list, 
               cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend, 
               double lod_act, double lod_light, cube corr_mat, cube lintegral_mat, mat vcovar_mat,
               cube lambda_act_mat, cube lambda_light_mat, bool tobit){

	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act_week.n_slices;
	List beta_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
	  
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);
		
		vec vcovar_vec = vcovar_mat.col(ind);
		
		for (int clust_i = 0; clust_i < num_re; clust_i++){

			Cube1.slice(clust_i) = BackwardIndC(act_ind, light_ind, tran_list, emit_act_week, emit_light_week, emit_act_weekend, 
               emit_light_weekend, clust_i, lod_act, lod_light, corr_mat, lintegral_mat, vcovar_vec,
               lambda_act_mat, lambda_light_mat, tobit);

		}

		beta_list(ind) = Cube1;
	}
	return(beta_list);
}  


// [[Rcpp::export]]
mat CalcTranHelperC(int init_state, int new_state, NumericMatrix act, NumericMatrix light, List tran_list_mat, 
                    cube emit_act_week, cube emit_light_week, cube emit_act_weekend, cube emit_light_weekend, 
                    NumericVector ind_like_vec, List alpha, List beta, double lod_act, double lod_light, 
                    cube corr_mat, cube lintegral_mat, mat pi_l, int clust_i, mat vcovar_mat,
                    cube lambda_act_mat, cube lambda_light_mat, bool tobit){
  int num_people = act.ncol();
  int len = act.nrow();
  
  mat tran_vals_re_mat( len-1, num_people );
  
  for (int ind = 0; ind < num_people; ind++) {
    
    vec vcovar_vec = vcovar_mat.col(ind);
    
    List tran_list_clust = tran_list_mat[clust_i];
    mat tran_mat = tran_list_clust(0);
    
    if (tran_list_clust.size() > 1){
      mat tran_mat_week = tran_list_clust(0);
      mat tran_mat_weekend = tran_list_clust(1);
      
      
      tran_mat_week.each_col() %= (1-vcovar_vec);
      tran_mat_weekend.each_col() %= vcovar_vec;
      tran_mat = tran_mat_week + tran_mat_weekend;
    } 
      

    
    //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
    int tran_vec_ind = init_state + (new_state * 2);
    
    arma::cube alpha_ind = alpha(ind); 
    arma::cube beta_ind = beta(ind);
    double likelihood = ind_like_vec(ind);
    
    NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
    NumericVector act_ind_m1 = act_ind.column(0); 
    
    NumericMatrix light_ind = light( Range(1,len-1) , Range(ind,ind) );
    NumericVector light_ind_m1 = light_ind.column(0); 
    
    vec class_vec_week = logClassificationC( act_ind, light_ind, 
                                             emit_act_week(new_state,0,clust_i), 
                                             emit_act_week(new_state,1,clust_i), 
                                             emit_light_week(new_state,0,clust_i), 
                                             emit_light_week(new_state,1,clust_i),
                                             lod_act, lod_light, corr_mat(clust_i,new_state,0), lintegral_mat(clust_i,new_state,0),
                                             lambda_act_mat(clust_i,new_state,0),lambda_light_mat(clust_i,new_state,0),tobit);
    
    vec class_vec_weekend = logClassificationC( act_ind, light_ind, 
                                                emit_act_weekend(new_state,0,clust_i), 
                                                emit_act_weekend(new_state,1,clust_i), 
                                                emit_light_weekend(new_state,0,clust_i), 
                                                emit_light_weekend(new_state,1,clust_i),
                                                lod_act, lod_light, corr_mat(clust_i,new_state,1), lintegral_mat(clust_i,new_state,1),
                                                lambda_act_mat(clust_i,new_state,1),lambda_light_mat(clust_i,new_state,1),tobit);
    
    vec class_vec = (class_vec_week % (1-vcovar_vec)) + (class_vec_weekend % vcovar_vec); 
    
    vec alpha_ind_slice = alpha_ind(span(0,len-2),span(init_state,init_state),span(clust_i,clust_i));
    vec beta_ind_slice = beta_ind(span(1,len-1),span(new_state,new_state),span(clust_i,clust_i));
    vec tran_vec_slice = tran_mat.col(tran_vec_ind);
    
    vec temp = alpha_ind_slice + beta_ind_slice + log(tran_vec_slice) + log(pi_l(ind,clust_i)) + class_vec - likelihood;
    vec tran_vals_re_ind = arma::exp(temp);
    
    tran_vals_re_mat.col(ind) = tran_vals_re_ind;
    
  }
  

  return tran_vals_re_mat;
}
	 
// [[Rcpp::export]]
mat ForwardIndAltC(mat decoded_mat, rowvec init, Rcpp::List tran_list, int clust_i, vec vcovar_vec){
  
  vec decoded_ind = decoded_mat.col(clust_i);
  mat alpha( decoded_ind.n_elem, 2 );
  
  alpha(0,0) = log(init(0));
  alpha(0,1) = log(init(1));
  
  List tran_list_clust = tran_list[clust_i]; 
  
  double log_class_0;
  double log_class_1;
  
  for (int i = 1; i < decoded_ind.n_elem; i++) {
    List tran_list_vcovar = tran_list_clust[vcovar_vec(i)];
    
    NumericMatrix tran = tran_list_vcovar[i%96];
    
    
    log_class_0 = 0;
    log_class_1 = 0;
    
    if (decoded_ind(i) == 0){
      log_class_0 = 0;
      log_class_1 = -9999999;
    }
    
    if (decoded_ind(i) == 1){
      log_class_0 = -9999999;
      log_class_1 = 0;
    }
    
    double fp_00 = alpha(i-1,0) + log(tran(0,0)) + log_class_0;
    
    double fp_10 = alpha(i-1,1) + log(tran(1,0)) + log_class_0;
    
    double fp_01 = alpha(i-1,0) + log(tran(0,1)) + log_class_1;
    
    double fp_11 = alpha(i-1,1) + log(tran(1,1)) + log_class_1;
    
    NumericVector fp_0 = NumericVector::create(fp_00,fp_10);
    NumericVector fp_1 = NumericVector::create(fp_01,fp_11);
    
    alpha(i,0) = logSumExpC(fp_0);
    alpha(i,1) = logSumExpC(fp_1);
  }
  
  return(alpha);
}


// [[Rcpp::export]]
List ForwardAltC(cube decoded_cube, mat init, List tran_list, mat vcovar_mat){
  int num_people = decoded_cube.n_cols;
  int len = decoded_cube.n_rows;
  int num_re = decoded_cube.n_slices;
  List alpha_list(num_people);
  
  
  for (int ind = 0; ind < num_people; ind++) {
    
    mat decoded_mat = decoded_cube.col(ind);
    
    arma::cube Cube1(len, 2, num_re);
    
    vec vcovar_vec = vcovar_mat.col(ind);
    
    for (int clust_i = 0; clust_i < num_re; clust_i++){
      
      arma::rowvec init_vec = init.row(clust_i);
      
      Cube1.slice(clust_i) = ForwardIndAltC(decoded_mat, init_vec, tran_list, clust_i, vcovar_vec);
      
    }
    
    alpha_list(ind) = Cube1;
  }
  return(alpha_list);
}
 




