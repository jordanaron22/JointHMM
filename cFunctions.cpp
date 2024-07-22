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
vec logClassificationC(NumericVector act, NumericVector light, double mu_act, double sig_act, double mu_light, double sig_light, double lod_act, double lod_light, double bivar_corr, double lintegral) {

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
			if(act_obs != lod_act & light_obs != lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::dnorm( light_obs, mu_light_cond, sig_light_cond, true ) ;
			}

			//CASE 2
			if(act_obs != lod_act & light_obs == lod_light){

				double mu_light_cond = CalcCondMeanC(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs);
				double sig_light_cond = CalcCondSigC(sig_light,bivar_corr);

				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true ) + R::pnorm( light_obs, mu_light_cond, sig_light_cond, true, true ) ;
			}

			//CASE 3
			if(act_obs == lod_act & light_obs != lod_light){

				double mu_act_cond = CalcCondMeanC(mu_act,sig_act,mu_light,sig_light,bivar_corr,light_obs);
				double sig_act_cond = CalcCondSigC(sig_act,bivar_corr);

				log_dens_vec(i) = R::pnorm( act_obs, mu_act_cond, sig_act_cond, true, true ) + R::dnorm( light_obs, mu_light, sig_light, true ) ;
			}

			//CASE 4
			if(act_obs == lod_act & light_obs == lod_light){
				log_dens_vec(i) = lintegral;
			}

			// Light Missing
		} else if (nan_vec_act(i) != true && nan_vec_light(i) == true){
			if(act_obs != lod_act){	
				log_dens_vec(i) = R::dnorm( act_obs, mu_act, sig_act, true );
			} else {
				log_dens_vec(i) = R::pnorm( act_obs, mu_act, sig_act, true, true );
			}
			
			// Act missing
		} else if((nan_vec_act(i) == true && nan_vec_light(i) != true)){
			if(light_obs != lod_light){	
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
mat ForwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, NumericVector init, Rcpp::List tran_list, cube emit_act, cube emit_light,int clust_i, double lod_act, double lod_light, NumericVector corr_vec, vec beta_vec, int event, double bline, double cbline, mat lintegral_mat, double log_sweight){

	mat alpha( act_ind.length(), 2 );

	vec log_class_0 = logClassificationC( act_ind, light_ind, emit_act(0,0,clust_i), emit_act(0,1,clust_i), emit_light(0,0,clust_i), emit_light(0,1,clust_i), lod_act, lod_light, corr_vec(0), lintegral_mat(clust_i,0));

	vec log_class_1 = logClassificationC( act_ind, light_ind, emit_act(1,0,clust_i), emit_act(1,1,clust_i), emit_light(1,0,clust_i), emit_light(1,1,clust_i), lod_act, lod_light, corr_vec(1), lintegral_mat(clust_i,1));

	double surv_comp = event * log(bline) + beta_vec[clust_i] - cbline * exp(beta_vec[clust_i]);

	alpha(0,0) = log(init(0)) + log_class_0[0] + surv_comp;
	alpha(0,1) = log(init(1)) + log_class_1[0] + surv_comp;
	
	List tran_time_list = tran_list[clust_i]; 
	
	for (int i = 1; i < act_ind.length(); i++) {
	  
	  NumericMatrix tran = tran_time_list[i%96];
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
mat BackwardIndC(const NumericVector& act_ind, const NumericVector& light_ind, Rcpp::List tran_list, cube emit_act, cube emit_light,int clust_i, double lod_act, double lod_light, NumericVector corr_vec, mat lintegral_mat){
  
  int n = act_ind.length(); 
  mat beta( n, 2 );
  
  vec log_class_0 = logClassificationC( act_ind, light_ind, emit_act(0,0,clust_i), emit_act(0,1,clust_i), emit_light(0,0,clust_i), emit_light(0,1,clust_i), lod_act, lod_light, corr_vec(0), lintegral_mat(clust_i,0));
  
  vec log_class_1 = logClassificationC( act_ind, light_ind, emit_act(1,0,clust_i), emit_act(1,1,clust_i), emit_light(1,0,clust_i), emit_light(1,1,clust_i), lod_act, lod_light, corr_vec(1), lintegral_mat(clust_i,1));
  
  beta(n-1,0) = log(1);
  beta(n-1,1) = log(1);
  
  List tran_time_list = tran_list[clust_i]; 
  
  for (int i = n-2; i >= 0; i--) {
    
    NumericMatrix tran = tran_time_list[(i+1)%96];
    
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

// [[Rcpp::export]]
List ForwardC(const NumericMatrix& act, const NumericMatrix& light, NumericMatrix init, List tran_list, cube emit_act, cube emit_light,double lod_act, double lod_light, NumericMatrix corr_mat, vec beta_vec, vec event_vec, vec bline_vec, vec cbline_vec, mat lintegral_mat, NumericVector log_sweights_vec){
	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act.n_slices;
	List alpha_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);
		double log_sweight = log_sweights_vec(ind);

		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			NumericVector corr_vec = corr_mat.row(clust_i);
			NumericVector init_vec = init.row(clust_i);

			Cube1.slice(clust_i) = ForwardIndC(act_ind, light_ind, init_vec, tran_list, emit_act, emit_light,clust_i, lod_act, lod_light, corr_vec, beta_vec, event_vec[ind], bline_vec[ind], cbline_vec[ind], lintegral_mat, log_sweight);

		}

		alpha_list(ind) = Cube1;
	}
	return(alpha_list);
}

// [[Rcpp::export]]
List BackwardC(const NumericMatrix& act, const NumericMatrix& light, List tran_list, cube emit_act, cube emit_light, double lod_act, double lod_light, NumericMatrix corr_mat, mat lintegral_mat){

	int num_people = act.ncol();
	int len = act.nrow();
	int num_re = emit_act.n_slices;
	List beta_list(num_people);
	
	for (int ind = 0; ind < num_people; ind++) {
		arma::cube Cube1(len, 2, num_re);
		NumericVector act_ind = act.column(ind);
		NumericVector light_ind = light.column(ind);
		
		for (int clust_i = 0; clust_i < num_re; clust_i++){
			NumericVector corr_vec = corr_mat.row(clust_i);

			Cube1.slice(clust_i) = BackwardIndC(act_ind, light_ind, tran_list, emit_act, emit_light, clust_i, lod_act, lod_light, corr_vec, lintegral_mat);

		}

		beta_list(ind) = Cube1;
	}
	return(beta_list);
}  


// [[Rcpp::export]]
cube CalcTranHelperC(int init_state, int new_state, NumericMatrix act, NumericMatrix light, List tran_list_mat, cube emit_act, cube emit_light, NumericVector ind_like_vec, List alpha, List beta, double lod_act, double lod_light, NumericMatrix corr_mat, mat lintegral_mat, vec pi_l){
  int num_people = act.ncol();
  int len = act.nrow();
  int num_re = emit_act.n_slices;

  arma::cube tran_vals_re_cube( len-1, num_people, num_re );

  for (int clust_i = 0; clust_i < num_re; clust_i++){

	mat tran_vals_re_mat( len-1, num_people );
	NumericVector corr_vec = corr_mat.row(clust_i);

    for (int ind = 0; ind < num_people; ind++) {
	  mat tran_mat = tran_list_mat(clust_i);
	  

	  //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
	  int tran_vec_ind = init_state + (new_state * 2);

	  arma::cube alpha_ind = alpha(ind); 
	  arma::cube beta_ind = beta(ind);
	  double likelihood = ind_like_vec(ind);
	  
	  NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
	  NumericVector act_ind_m1 = act_ind.column(0); 

	  NumericMatrix light_ind = light( Range(1,len-1) , Range(ind,ind) );
	  NumericVector light_ind_m1 = light_ind.column(0); 
	  
	  vec class_vec = logClassificationC( act_ind_m1, light_ind_m1, emit_act(new_state,0,clust_i), emit_act(new_state,1,clust_i), emit_light(new_state,0,clust_i), emit_light(new_state,1,clust_i), lod_act, lod_light, corr_vec(new_state), lintegral_mat(clust_i,new_state) );
	  
	  vec alpha_ind_slice = alpha_ind(span(0,len-2),span(init_state,init_state),span(clust_i,clust_i));
	  vec beta_ind_slice = beta_ind(span(1,len-1),span(new_state,new_state),span(clust_i,clust_i));
	  vec tran_vec_slice = tran_mat.col(tran_vec_ind);

	  vec temp = alpha_ind_slice + beta_ind_slice + log(tran_vec_slice) + log(pi_l(clust_i)) + class_vec - likelihood;
	  vec tran_vals_re_ind = arma::exp(temp);

	  tran_vals_re_mat.col(ind) = tran_vals_re_ind;
	
    }

	tran_vals_re_cube.slice(clust_i) = tran_vals_re_mat;
    
  }

  return tran_vals_re_cube;
}


// // [[Rcpp::export]]
// cube CalcTranIndHelperC(int init_state, int new_state, NumericMatrix act, List tran_list_mat, NumericVector tran_ind_vec, cube emit_act, NumericVector ind_like_vec, List alpha, List beta, double lepsilon, NumericVector act_light_binom, vec pi_l){
//   int num_people = act.ncol();
//   int len = act.nrow();
//   int num_re = 1;
//   int clust_i = 0;

//   arma::cube tran_vals_re_cube( len-1, num_people, num_re );


// 	mat tran_vals_re_mat( len-1, num_people );

//     for (int ind = 0; ind < num_people; ind++) {
      
//       int tran_ind = tran_ind_vec(ind);
// 	  mat tran_mat = tran_list_mat(tran_ind-1);
	  

// 	  //0,0->0 & 1,0->1 & 0,1->2 & 1,1->3
// 	  int tran_vec_ind = init_state + (new_state * 2);

// 	  arma::cube alpha_ind = alpha(ind); 
// 	  arma::cube beta_ind = beta(ind);
// 	  double likelihood = ind_like_vec(ind);
	  



// 	  NumericMatrix act_ind = act( Range(1,len-1) , Range(ind,ind) );
// 	  NumericVector act_ind_m1 = act_ind.column(0); 
	  
	  

