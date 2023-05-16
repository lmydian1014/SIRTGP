# include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
# include <math.h>
# include <list>
# include <iostream>
# include <vector>
# include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
double sample_cpp(vec x, int size, bool replace, vec prob){
	// Obtaining namespace of truncnorm package
	Environment pkg = Environment::namespace_env("base");
	// Picking up rtruncnorm() function from truncnorm package
	Function f = pkg["sample"];
	NumericVector r = f(x, size, replace, prob);
	return as<double>(r);
}

// [[Rcpp::export]]
vec quantile_cpp(vec x, double size, double low, double high){
	// Obtaining namespace of truncnorm package
	Environment pkg = Environment::namespace_env("stats");
	// Picking up rtruncnorm() function from truncnorm package
	Function f = pkg["quantile"];
	double r1 = as<double>(f(x, low));
	double r2 = as<double>(f(x, high));
	vec r = linspace(r1, r2, size);
	return r;
}

// [[Rcpp::export]]
double rtrunc_Rpkg(double a, double b, double mu, double sd){
	// Obtaining namespace of truncnorm package
	Environment pkg = Environment::namespace_env("truncnorm");

	// Picking up rtruncnorm() function from truncnorm package
	Function f = pkg["rtruncnorm"];
	NumericVector r = f(1, a, b, mu, sd);
	return as<double>(r);
}


// [[Rcpp::export]]
mat fun_mul(mat A, colvec x) 
    { 
        A.each_col() %= x;
        return A;
    }
// [[Rcpp::export]]
mat fun_dev(mat A, colvec x) 
    { 	
    	Rcout << size(A) << endl;
    	Rcout << size(x) << endl;
        A.each_col() /= x;
        return A;
    }
// [[Rcpp::export]]
vec G_thres_cpp(vec x, double thres){
	vec out(x.n_elem);
	for(int i=0; i < x.n_elem; i++){
		if(x(i) > thres){
			out(i) = 1;
		}
		else{
			out(i) = 0;
		}
	}
	return(out);
}

// [[Rcpp::export]]
double sample_e_cpp(int k, int l, int K, int rt, vec Z, mat X, mat Xmat, vec SR_other, vec S0, vec E_hat, double thres1, vec lambda, double tausq){
	// SM, SR matrix K*n
	vec S = S0 + SR_other;
	uvec ind = regspace<uvec>(k*rt,(k+1)*rt -1);  
	//vec T = sum(fun_mul(X.rows(ind), Xmat.rows(ind).eval().col(l) % G_thres_cpp(abs(E_hat(ind)), thres1)) ,0).t();

	vec T = sum(fun_mul(X.rows(ind), Xmat.col(l) % G_thres_cpp(abs(E_hat(ind)), thres1)) ,0).t()/(K*rt);

	double B = (-2) * 100000*lambda(l) * sum((Z-S) % T);
	double A = 100000*lambda(l) * sum(square(T)) + tausq;
	return(randn(distr_param(-B/(2.0*A), sqrt(tausq*100000*lambda(l)/A))));

}


// [[Rcpp::export]]
double sample_ic_cpp(double n, vec Z,  vec SR_other, vec S0, double tausq){
	// SM, SR matrix K*n
	vec S = S0 + SR_other;
	double B = (-2) * 0.001 * sum((Z-S));
	double A = 0.001 * n + tausq;
	return(randn(distr_param(-B/(2.0*A), sqrt(tausq*0.001/A))));

}


// [[Rcpp::export]]
double sample_eta_cpp(int u, int v, int K, int rt, vec Z, mat X, mat X0, vec SR, vec S0_other, vec E_hat, vec eta_hat, double thres1, double thres2, double tausq, double sigsq_eta){
	vec S = S0_other + SR;
	int ind0 = (2*K-3-u)*u/2 + v-1;
	double V0 = K*(K-1)/2;
	vec T1 = (abs(eta_hat(ind0)) > thres2) * X0.row(ind0).t();
	vec T = T1/V0;
	double B = (-2) * sigsq_eta * sum((Z-S) % T);
	double A = sigsq_eta * sum(square(T)) + tausq;
	return(randn(distr_param(-B/(2.0*A), sqrt(tausq * sigsq_eta/A))));
}



// [[Rcpp::export]]
double sample_E_hat_cpp(int k, int j, int K, int rt, vec Z, mat X, mat Xmat, vec SR_other, vec S0, mat e, mat eta_m, double thres1, double tausq, double sigsq){
	
	vec S = S0 + SR_other;
	uvec ind = regspace<uvec>(k*rt,(k+1)*rt -1); 

	//double mu = as_scalar(Xmat.rows(ind).eval().row(j) * e.col(k)) + sum(eta_m.row(k));
	double mu = as_scalar(Xmat.row(j) * e.col(k));

	vec T = mu * X.rows(ind).eval().row(j).t()/(K*rt);
	double res = 0;


	if(thres1 > 0){
		vec L_trans(3, fill::zeros);
		vec L(3, fill::zeros);
		double Lmin = R::pnorm((-1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lnegthres = R::pnorm((-thres1 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lposthres = R::pnorm((thres1 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lmax = R::pnorm((1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		//Rcout << "hello" << endl;
		L[1] = Lposthres + log(1-exp(Lnegthres - Lposthres)) - sum(square(Z - S))/(2*tausq);
		L[0] = Lnegthres + log(1-exp(Lmin - Lnegthres)) - sum(square(T - (Z - S)))/(2*tausq);
		L[2] = Lmax + log(1-exp(Lposthres - Lmax)) - sum(square(T - (Z - S)))/(2*tausq);
		double L0 = max(L);
		L_trans = exp(L - L0);

		vec W = L_trans/sum(L_trans);
		double rand = conv_to<double>::from(randu(1));
		vec interval(4);
		interval(0) = R_NegInf; interval(1) = -thres1; interval(2) = thres1; interval(3) = R_PosInf;
		
		if(rand < W(0)){
			res = rtrunc_Rpkg(interval(0), interval(1), mu, sqrt(sigsq));
		}
		else{
			if(rand < W(0) + W(1)){
				res = rtrunc_Rpkg(interval(1), interval(2), mu, sqrt(sigsq));
			}
			else{
				res = rtrunc_Rpkg(interval(2), interval(3), mu, sqrt(sigsq));
			}
		}
	}
	else{
		vec L_trans(2, fill::zeros);
		vec L(2, fill::zeros);
		double Lmin = R::pnorm((-1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lposthres = R::pnorm((thres1 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lmax = R::pnorm((1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		L[0] = Lposthres + log(1-exp(Lmin - Lposthres)) - sum(square(T - (Z - S)))/(2*tausq);
		L[1] = Lmax + log(1-exp(Lposthres - Lmax)) - sum(square(T - (Z - S)))/(2*tausq);
		double L0 = max(L);
		L_trans = exp(L - L0);
		vec W = L_trans/sum(L_trans);
		double rand = conv_to<double>::from(randu(1));
		vec interval(3);
		interval(0) = R_NegInf; interval(1) = thres1; interval(2) = R_PosInf;
		
		if(rand < W(0)){
			res = rtrunc_Rpkg(interval(0), interval(1), mu, sqrt(sigsq));
		}
		else{
			res = rtrunc_Rpkg(interval(1), interval(2), mu, sqrt(sigsq));	
		}
	}
	//Rcout << res << endl;
	return res;
	//}
	
}	


// [[Rcpp::export]]
double sample_eta_hat_cpp(int u, int v, int K, int rt, vec Z, mat X, mat X0, mat Xmat, vec SR, vec S0_other, vec eta, double thres2, double tausq, double sigsq){
	Environment pkg = Environment::namespace_env("truncnorm");
	Function f_sample = pkg["rtruncnorm"];
	vec S = SR + S0_other;
	int ind0 = (2*K-3-u)*u/2 + v-1;
	double V0 = K*(K-1)/2;

	double mu = eta(ind0);
	vec T = mu * X0.row(ind0).t()/V0;
	

	double res = 0;
	if(thres2 > 0){
		vec L_trans(3, fill::zeros);
		vec L(3, fill::zeros);
		double Lmin = R::pnorm((-1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lnegthres = R::pnorm((-thres2 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lposthres = R::pnorm((thres2 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lmax = R::pnorm((1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		L[1] = Lposthres + log(1-exp(Lnegthres - Lposthres)) - sum(square(Z - S))/(2*tausq);
		L[0] = Lnegthres + log(1-exp(Lmin - Lnegthres)) - sum(square(T - (Z - S)))/(2*tausq);
		L[2] = Lmax + log(1-exp(Lposthres - Lmax)) - sum(square(T - (Z - S)))/(2*tausq);
		double L0 = max(L);
		L_trans = exp(L - L0);
		vec W = L_trans/sum(L_trans);
		double rand = conv_to<double>::from(randu(1));
		vec interval(4);
		//interval(0) = -thres2-1e9; interval(1) = -thres2; interval(2) = thres2; interval(3) = thres2+1e9;
		interval(0) = R_NegInf; interval(1) = -thres2; interval(2) = thres2; interval(3) = R_PosInf;
		
		if(rand < W(0)){
			//res = as<double>(rtruncnorm(1, mu, sqrt(sigsq), interval(0), interval(1)));
			//res = as<double>(f_sample(1, interval(0), interval(1), mu, sqrt(sigsq)));
			res = rtrunc_Rpkg(interval(0), interval(1), mu, sqrt(sigsq));
		}
		else{
			if(rand < W(0) + W(1)){
				//res = as<double>(rtruncnorm(1, mu, sqrt(sigsq), interval(1), interval(2)));
				res = rtrunc_Rpkg(interval(1), interval(2), mu, sqrt(sigsq));
			}
			else{
				//res = as<double>(rtruncnorm(1, mu, sqrt(sigsq), interval(2), interval(3)));
				res = rtrunc_Rpkg(interval(2), interval(3), mu, sqrt(sigsq));
			}
		}
	}
	else{
		vec L_trans(2, fill::zeros);
		vec L(2, fill::zeros);
		double Lmin = R::pnorm((-1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lposthres = R::pnorm((thres2 - mu)/sqrt(sigsq), 0, 1, true, true);
		double Lmax = R::pnorm((1e9 - mu)/sqrt(sigsq), 0, 1, true, true);
		L[0] = Lposthres + log(1-exp(Lmin - Lposthres)) - sum(square(T - (Z - S)))/(2*tausq);
		L[1] = Lmax + log(1-exp(Lposthres - Lmax)) - sum(square(T - (Z - S)))/(2*tausq);
		double L0 = max(L);
		L_trans = exp(L - L0);
		vec W = L_trans/sum(L_trans);
		//Rcout << W << endl;
		//Rcout << mu << endl;
		double rand = conv_to<double>::from(randu(1));
		vec interval(3);
		interval(0) = R_NegInf; interval(1) = thres2; interval(2) = R_PosInf; ;
		
		if(rand < W(0)){
			//res = as<double>(rtruncnorm(1, mu, sqrt(sigsq), interval(0), interval(1)));
			res = rtrunc_Rpkg(interval(0), interval(1), mu, sqrt(sigsq));
		}
		else{
			//res = as<double>(rtruncnorm(1, mu, sqrt(sigsq), interval(1), interval(2)));
			res = rtrunc_Rpkg(interval(1), interval(2), mu, sqrt(sigsq));
		}
	}
	
	return res;
}





// [[Rcpp::export]]
double loglike_cpp(int n, vec Y, vec SR, vec S0, double tausq){
	vec mu = SR + S0;
	vec loglike = Y % log(normcdf(mu)) + ((1-Y) % log(1 - normcdf(mu)));
	return (sum(loglike));
}



// [[Rcpp::export]]
vec sample_Z(int n, vec Y, vec SR, vec S0){
	Environment pkg = Environment::namespace_env("truncnorm");
	Function f_sample = pkg["rtruncnorm"];
	vec mu = SR + S0;
	vector<double> Z;
	double r;
	for(int i=0; i < n; i++){
		
		if(Y(i) > 0){
			r = rtrunc_Rpkg(0, R_PosInf, mu(i), 1.0);
			Z.push_back(r);
		}
		else{
			r = rtrunc_Rpkg(R_NegInf, 0, mu(i), 1.0);
			Z.push_back(r);
		}
	}

	return(Z);
}



// [[Rcpp::export]]
vec sample_thres1(int n, int K, vec Z, mat X, mat Xmat, vec val, vec prob1, mat eta_m, mat e, vec E_hat, vec S0, double tausq, double rt){
	vec prob_new(val.n_elem, fill::zeros);
	vec mu;
	vec SR(n, fill::zeros);
	vec SM(n, fill::zeros);
		
	for(int cand = 0; cand < val.n_elem; cand ++){
		
		double temp = val(cand);
		
		for(int iter = 0; iter < K; iter++){
			uvec ind = regspace<uvec>(iter*rt,(iter+1)*rt -1); 
			//SR = SR + sum(fun_mul(X.rows(ind), (Xmat.rows(ind) * e.col(iter)) % G_thres_cpp(abs(E_hat(ind)), thres)).t(), 1);
			SR = SR + sum(fun_mul(X.rows(ind), (Xmat * e.col(iter)) % G_thres_cpp(abs(E_hat(ind)), temp)).t(), 1);		
			//SM = SM + sum(eta_m.row(iter)) * sum(fun_mul(X.rows(ind), G_thres_cpp(abs(E_hat(ind)), temp)) ,0).t();
		}
		SR = SR/(K*rt);

		mu = S0 + SR;
		prob_new(cand) = (-sum(square(Z-mu))/(2*tausq)) + log(prob1(cand));
		//prob_new(cand) = (exp(-sum(square(Z-mu))/(2*tausq)) * prob1(cand));
	}
	double max_prob_new = max(prob_new);
	vec prob_update = exp(prob_new - max_prob_new);
	prob_update = prob_update/sum(prob_update);
	return(prob_update);
}


// [[Rcpp::export]]
vec sample_thres2(int n, int K, vec Z, mat X0, vec val, vec prob2, mat eta, vec eta_hat, vec SR, double tausq){
	vec prob_new(val.n_elem, fill::zeros);
	vec S0;
	vec mu;
	double V0 = K*(K-1)/2;
	for(int cand = 0; cand < val.n_elem; cand ++){
		double temp = val(cand);
		S0 = sum(fun_mul(X0, eta % G_thres_cpp(abs(eta_hat), temp)), 0).t();
		mu = S0/V0 + SR;
		prob_new(cand) = -sum(square(Z-mu))/(2*tausq) - log(1/prob2(cand));
	}
	
	double max_prob_new = max(prob_new);
	vec prob_update = exp(prob_new - max_prob_new);
	
	prob_update = prob_update/sum(prob_update);
	return(prob_update);
}

// [[Rcpp::export]]
List gibbs_sample_cpp(int T, int K, int L, int n, vec Y, mat X, mat X0, mat Xmat, vec eta, mat eta_m, mat e, vec E_hat, vec eta_hat, double beta0,
	double thres1, double thres2, vec lambda, double tausq, double sigsq, double sigsq_eta, int rt, double alpha1, double alpha2, vec prob1, vec prob2, vec val1, vec val2){

	//Clock clock;
	uvec ind;
	double V = X.n_rows;
	double V0 = K*(K-1)/2;
	vec SR(n, fill::zeros);

	mat SR_record(n, T);
	mat S0_record(n, T);

	for(int iter = 0; iter < K; iter++){
		uvec ind = regspace<uvec>(iter*rt,(iter+1)*rt -1); 
		SR = SR + sum(fun_mul(X.rows(ind), (Xmat * e.col(iter)) % G_thres_cpp(abs(E_hat(ind)), thres1)).t(), 1);		
	}
	vec S0 = sum(fun_mul(X0, eta % G_thres_cpp(abs(eta_hat), thres2)), 0).t();
	SR = SR/(K*rt) + beta0;
	//SM = SM/(K*rt);
	S0 = S0/V0;

	SR_record.col(0) = SR;
	S0_record.col(0) = S0;

	
	vec Z_temp = SR+S0;
	//vec Z_temp = Z;
	mat Z_record(n,T);
	Z_record.col(0) = Z_temp;
	/////////////////////  initialize parameters
	
	mat e_record(L,T*K);
	mat eta_record(V0,T);
	mat E_hat_record(V,T);
	mat eta_hat_record(V0,T);
	vector<double> loglike_record;
	vector<double> tausq_record;
	vector<double> thres1_record;
	vector<double> thres2_record;
	vector<double> beta0_record;

	vec eta_temp = eta;
	eta_record.col(0) = eta;

	mat e_temp = e;
	e_record.cols(span(0,K-1)) = e;

	vec E_hat_temp = E_hat;
	E_hat_record.col(0) = E_hat;

	mat eta_m_temp = eta_m;

	vec eta_hat_temp = eta_hat;
	eta_hat_record.col(0) = eta_hat;

	double temp_loglike = loglike_cpp(n, Y, SR, S0, tausq);
	loglike_record.push_back(temp_loglike);

	double tausq_temp = tausq; 
	tausq_record.push_back(tausq_temp);

	double thres1_temp = thres1; 
	thres1_record.push_back(thres1_temp);
	double thres2_temp = thres2; 
	thres2_record.push_back(thres2_temp);


	double beta0_temp = beta0;
	beta0_record.push_back(beta0_temp);

	vec prob1_temp = prob1;
	vec prob2_temp = prob2;
	vec val1_temp = val1;
	vec val2_temp = val2;

	Rcout << "hello1" << endl;
	Rcout << temp_loglike << endl;
	Rcout << "hello2" << endl;

	//clock.tick("time");
	//vec SM_other;
	vec S0_other;
	vec SR_other;
	
	for(int t=0; t<(T-1); t++){
		if(t % 100 == 0){
			Rcout << t << endl;
		}
		
		for(int k=0; k<K; k++){
			
			ind = regspace<uvec>(k*rt,(k+1)*rt -1); 
			//------------------------************** Sample e **************------------------------//
			for(int l=0; l<L; l++){
				//Rcout << l << endl;
        	    vec SR_temp = sum(fun_mul(X.rows(ind), Xmat.col(l) % G_thres_cpp(abs(E_hat_temp(ind)), thres1_temp)).t(), 1);
        	    SR_other = SR - SR_temp * e_temp(l,k)/(K*rt);
        		e_temp(l,k) = sample_e_cpp(k, l, K, rt, Z_temp, X, Xmat, SR_other, S0, E_hat_temp, thres1_temp, lambda, tausq_temp);
        		e_record(l, (t+1)*K + k) = e_temp(l,k);
        		SR = SR_other + SR_temp * e_temp(l,k)/(K*rt);
			}
			//------------------------************** Sample E_hat **************------------------------//
			for(int j=0; j<rt; j++){
			    SR_other = SR - as_scalar((Xmat.row(j) * e_temp.col(k)) * (abs(E_hat_temp(k*rt + j)) > thres1_temp))/(K*rt) *  X.row(k*rt + j).t(); //should be a scalar * vec
			    E_hat_temp(k*rt + j) = sample_E_hat_cpp(k, j, K, rt, Z_temp, X, Xmat, SR_other, S0, e_temp, eta_m_temp, thres1_temp, tausq_temp, sigsq);
			    E_hat_record(k*rt + j, t+1) = E_hat_temp(k*rt + j);
				SR = SR_other + as_scalar((Xmat.row(j) * e_temp.col(k)) * (abs(E_hat_temp(k*rt + j)) > thres1_temp))/(K*rt) *  X.row(k*rt + j).t(); //should be a scalar * vec 	
			}
		
		}
		for(int u = 0; u < (K-1); u++){
			for(int v = (u+1); v < K; v++){
				//------------------------************** Sample eta **************------------------------//
				uvec ind_u = regspace<uvec>(u*rt,(u+1)*rt -1); 
				uvec ind_v = regspace<uvec>(v*rt,(v+1)*rt -1); 
				int index = ((2*K-3-u)*u)/2 + v-1;
				vec SM_tmp2 = (abs(eta_hat_temp(index)) > thres2_temp) * X0.row(index).t();
				S0_other = S0 - eta_temp(index)/V0 * SM_tmp2;
				eta_temp(index) = sample_eta_cpp(u, v, K, rt, Z_temp, X, X0, SR, S0_other, E_hat_temp, eta_hat_temp, thres1_temp, thres2_temp, tausq_temp, sigsq_eta);
				S0 = S0_other + eta_temp(index)/V0 * SM_tmp2;

				eta_record(index, t+1) = eta_temp(index);
				eta_m_temp(u,v) = eta_temp(index);// Update eta_m_temp
				eta_m_temp(v,u) = eta_temp(index);// Update eta_m_temp
				
				//------------------------************** Sample eta_hat **************------------------------//
				S0_other = S0 - eta_temp(index)/V0 * (abs(eta_hat_temp(index)) > thres2_temp) * X0.row(index).t();
				eta_hat_temp(index) = sample_eta_hat_cpp(u, v, K, rt, Z_temp, X, X0, Xmat, SR, S0_other, eta_temp, thres2_temp, tausq_temp, sigsq);
				eta_hat_record(index,t+1) = eta_hat_temp(index);
				S0 = S0_other + eta_temp(index)/V0 * (abs(eta_hat_temp(index)) > thres2_temp) * X0.row(index).t();
			}
		}
		//if(t < 200){
		//SR_other = SR - beta0_temp;
		//beta0_temp = sample_ic_cpp(n, Z_temp, SR_other, S0, tausq_temp);
		beta0_record.push_back(beta0_temp);
		//SR = SR_other + beta0_temp;
		/*}
		else{
			beta0_record.push_back(beta0_temp);
		}*/
		
		//------------------------************** Sample Z **************------------------------//
		Z_temp = sample_Z(n, Y, SR, S0);
		Z_record.col(t+1) = Z_temp;
		//------------------------************** Sample thres1 **************------------------------//
		if(t > 100){
			val1_temp = quantile_cpp(abs(E_hat_temp), prob1_temp.n_elem, 0.25, 0.9);
			//Rcout << val1_temp <	< endl;
			prob1_temp = sample_thres1(n, K, Z_temp, X, Xmat, val1_temp, prob1_temp, eta_m_temp, e_temp, E_hat_temp, S0, tausq_temp, rt);
			double thres1_temp = sample_cpp(val1_temp, 1, TRUE, prob1_temp);
			thres1_record.push_back(thres1_temp);
			for(int i=0; i<n; i++){
				SR(i) = 0;
				//SM(i) = 0;
			}
			for(int iter = 0; iter < K; iter++){
				uvec ind = regspace<uvec>(iter*rt,(iter+1)*rt -1); 
				SR = SR + sum(fun_mul(X.rows(ind), (Xmat * e_temp.col(iter)) % G_thres_cpp(abs(E_hat_temp(ind)), thres1_temp)).t(), 1);		
			}
			SR = SR/(K*rt);
			//------------------------************** Sample thres2 **************------------------------//
			val2_temp = quantile_cpp(abs(eta_hat_temp), prob2_temp.n_elem, 0.25, 0.9);
			//Rcout << val2_temp << endl;
			prob2_temp = sample_thres2(n, K, Z_temp, X0, val2_temp, prob2_temp, eta_temp, eta_hat_temp, SR, tausq_temp);
			double thres2_temp = sample_cpp(val2_temp, 1, TRUE, prob2_temp);
			thres2_record.push_back(thres2_temp);
			S0 = sum(fun_mul(X0, eta_temp % G_thres_cpp(abs(eta_hat_temp), thres2_temp)), 0).t();
			S0 = S0/V0;
			/*thres1_temp = 40;
			thres2_temp = 40;
			thres1_record.push_back(thres1_temp);
			thres2_record.push_back(thres2_temp);*/
		}
		else{
			thres1_record.push_back(thres1_temp);
			thres2_record.push_back(thres2_temp);
		}
		
		//------------------------************** loglike **************------------------------//
		double temp_loglike = loglike_cpp(n, Y, SR, S0, tausq_temp);
		loglike_record.push_back(temp_loglike);

		if(t % 100 == 0){
			Rcout << temp_loglike << endl;
		}

		SR_record.col(t+1) = SR;
		S0_record.col(t+1) = S0;

	}
	/*clock.tock("time");
	clock.stop("naptimes");*/
	Rcout << "finish" << endl;
	List chain = List::create(Named("e") = e_record , _["eta"] = eta_record, _["E_hat"] = E_hat_record, _["eta_hat"] = eta_hat_record, _["loglike"] = loglike_record, 
		_["Z"] = Z_record, _["SR"] = SR_record, _["S0"] = S0_record, _["thres1"] = thres1_record, _["thres2"] = thres2_record, _["intecept"] = beta0_record);
	//e : L*(K*T); eta and eta_hat: V0(120)*T; E_hat: rt*K(416)*T; 
	return chain;

}

















