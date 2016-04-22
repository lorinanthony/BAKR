// load Rcpp
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
//random multivariate normal sample generator using RcppArmadillo
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
double rgammadouble(int a, double b, double c){
  Rcpp::NumericVector x = rgamma(a,b,1/c);
  return x(0);
}

// [[Rcpp::export]]
double rInvGamma(int n, double shape,double scale){
  double x = rgammadouble(n, shape, scale);
  return 1.0/x;
}

// [[Rcpp::export]]
double rScaledInvChiSq(int n, double nu, double tau2){
  double x = rInvGamma(n, nu/2,(nu * tau2)/2);
  return x;
}

// [[Rcpp::export]]
List BAKRGibbs(mat U, vec y, mat Lambda, int iter = 5000, int burn = 1000, double g = 1, double nu_alpha = 2, double phi_alpha = 1, double nu_e = 2, double phi_e = 1){
    
    int i;
    double q = U.n_cols;
    double n = U.n_rows;
    mat alphamat = zeros<mat>(iter-burn, q);
    vec sigma_e_vec = zeros(iter-burn), sigma_alpha_vec = zeros(iter-burn);
    
    // The rest of the code follows the R version
    rowvec alpha = zeros<rowvec>(q);
    double sigma_alpha = 0.5, sigma_e = 0.5;
    mat I_q(q,q); I_q.eye();
    arma::mat SS;
    
    for (i=0; i<iter; i++) {
        SS = inv(diagmat(sigma_e*Lambda+sigma_alpha*g*I_q));
        alpha = mvrnormArma(1,sigma_alpha*g*SS*trans(U)*y, (sigma_e*sigma_alpha*g)*SS);
        
        sigma_alpha = rScaledInvChiSq(1,nu_alpha + q,(1.0/(nu_alpha + q))*(nu_alpha*phi_alpha+as_scalar(alpha*Lambda*trans(alpha))/g));
        
        sigma_e = rScaledInvChiSq(1,nu_e + n,(1.0/(nu_e + n))*(nu_e*phi_e+as_scalar(trans(y - U*trans(alpha))*(y - U*trans(alpha)))));
        
        if(i>=burn){
            alphamat.row(i-burn) = alpha;
            sigma_alpha_vec(i-burn) = sigma_alpha;
            sigma_e_vec(i-burn) = sigma_e;
        }
    }
    return Rcpp::List::create(Rcpp::Named("alpha") = alphamat, Rcpp::Named("sigmaE") = sigma_e_vec, Rcpp::Named("sigmaA") = sigma_alpha_vec);
}

// [[Rcpp::export]]
arma::rowvec PostMean(arma::mat postdraws){
  arma::rowvec x = mean(postdraws);
  return x;
}

// [[Rcpp::export]]
arma::vec PostBeta(arma::mat GAMMA, arma::vec alpha){
  return GAMMA*alpha;
}

// [[Rcpp::export]]
arma::mat GetBeta(arma::mat GAMMA, arma::mat alphamat){
  int i;
  mat betamat(GAMMA.n_rows,alphamat.n_rows);
  for(i=0; (unsigned)i<betamat.n_cols; i++){
    betamat.col(i) = GAMMA*trans(alphamat.row(i));
  }
  return trans(betamat);
}

// [[Rcpp::export]]
arma::mat GetApproxKernel(arma::mat X, double iter){
  int i;
  double ncov = X.n_rows, samp_size = X.n_cols;
  mat zeta(iter,samp_size);
  for(i = 0; i<iter; i++){
    vec omega = arma::randn(ncov)*sqrt(2/ncov);
    vec b = runif(1,0,2*datum::pi);
    zeta.row(i) = sqrt(2/iter)*cos(trans(omega)*X+as_scalar(b));
  }
  mat K_hat = zeta.t()*zeta;
  return K_hat;
}

// [[Rcpp::export]]
arma::mat GetGaussKernel(arma::mat X, double h = 1){
    int i,j;
    double n = X.n_cols;
    double p = X.n_rows;
    mat K = zeros<mat>(n,n);
    for (i = 0; i<n; i++){
        for(j = 0; j<n; j++){
            if(i==j){
                break;
            }
            else{
                K(i,j) = exp(-h/(p)*sum(pow(X.col(i)-X.col(j),2)));
            }
        }
    }
    return K + trans(K);
}


// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat X){
    double p = X.n_rows;
    return trans(X)*X/p;
}

// [[Rcpp::export]]
arma::mat InverseMap(arma::mat X, arma::mat U){
    mat X_t = trans(X);
    mat GAMMA = pinv(X_t,1.490116e-08)*U;
    return GAMMA;
}

// [[Rcpp::export]]
double median_rcpp(NumericVector x) {
  NumericVector y = clone(x);
  int n, half;
  double y1, y2;
  n = y.size();
  half = n / 2;
  if(n % 2 == 1) {
    // median for odd length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    return y[half];
  } else {
    // median for even length vector
    std::nth_element(y.begin(), y.begin()+half, y.end());
    y1 = y[half];
    std::nth_element(y.begin(), y.begin()+half-1, y.begin()+half);
    y2 = y[half-1];
    return (y1 + y2) / 2.0;
  }
}

// [[Rcpp::export]]
double mad_rcpp(NumericVector x, double scale_factor = 1.4826) {
  // scale_factor = 1.4826; default for normal distribution consistent with R
  return median_rcpp(abs(x - median_rcpp(x))) * scale_factor;
}

// [[Rcpp::export]]
vec BAKRPredict(mat X_new, vec beta) {
  return trans(X_new)*beta;
}

// [[Rcpp::export]]
List BAKRSupervised(mat K, vec y, mat Sigma, int iter = 5000, int burn = 1000, double nu_alpha = 2, double phi_alpha = 1, double nu_e = 2, double phi_e = 1){
  
  int i;
  double q = K.n_cols;
  double n = K.n_rows;
  mat alphamat = zeros<mat>(iter-burn, q);
  mat yfitmat = zeros<mat>(iter-burn, n);
  
  // The rest of the code follows the R version
  rowvec alpha = zeros<rowvec>(q);
  double sigma_alpha = 0.5, sigma_e = 0.5;
  mat I_q(q,q); I_q.eye();
  
  arma::mat Lambda = Sigma;
  arma::mat U = K;
  
  for (i=0; i<iter; i++) {
    vec yfit = U*trans(alpha);
    vec ystar = yfit;
    ystar(find_finite(y)) = y(find_finite(y));
    
    alpha = mvrnormArma(1,sigma_alpha*inv(sigma_e*inv(Lambda)+sigma_alpha*I_q)*trans(U)*ystar, (sigma_e*sigma_alpha)*inv(sigma_e*inv(Lambda)+sigma_alpha*I_q));
    
    sigma_alpha = rScaledInvChiSq(1,nu_alpha + q,(1.0/(nu_alpha + q))*(nu_alpha*phi_alpha+as_scalar(alpha*inv(Lambda)*trans(alpha))));
    
    sigma_e = rScaledInvChiSq(1,nu_e + n,(1.0/(nu_e + n))*(nu_e*phi_e+as_scalar(trans(ystar - U*trans(alpha))*(ystar - U*trans(alpha)))));
    
    if(i>=burn){
      alphamat.row(i-burn) = alpha;
      yfitmat.row(i-burn) = trans(ystar);
    }
  }
  return List::create(Named("alpha") = alphamat,Named("YFit") = yfitmat);
}

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double norm_rs(double a, double b)
{
    double  x;
    x = Rf_rnorm(0.0, 1.0);
    while( (x < a) || (x > b) ) x = norm_rand();
    return x;
}

// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double half_norm_rs(double a, double b)
{
    double   x;
    x = fabs(norm_rand());
    while( (x<a) || (x>b) ) x = fabs(norm_rand());
    return x;
}

// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double unif_rs(double a, double b)
{
    double xstar, logphixstar, x, logu;
    
    // Find the argmax (b is always >= 0)
    // This works because we want to sample from N(0,1)
    if(a <= 0.0) xstar = 0.0;
    else xstar = a;
    logphixstar = R::dnorm(xstar, 0.0, 1.0, 1.0);
    
    x = R::runif(a, b);
    logu = log(R::runif(0.0, 1.0));
    while( logu > (R::dnorm(x, 0.0, 1.0,1.0) - logphixstar))
    {
        x = R::runif(a, b);
        logu = log(R::runif(0.0, 1.0));
    }
    return x;
}

// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double exp_rs(double a, double b)
{
    double  z, u, rate;
    
    //  Rprintf("in exp_rs");
    rate = 1/a;
    //1/a
    
    // Generate a proposal on (0, b-a)
    z = R::rexp(rate);
    while(z > (b-a)) z = R::rexp(rate);
    u = R::runif(0.0, 1.0);
    
    while( log(u) > (-0.5*z*z))
    {
        z = R::rexp(rate);
        while(z > (b-a)) z = R::rexp(rate);
        u = R::runif(0.0,1.0);
    }
    return(z+a);
}

// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

// [[Rcpp::export]]
double rnorm_trunc (double mu, double sigma, double lower, double upper)
{
    int change;
    double a, b;
    double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
    double z, tmp, lograt;
    
    change = 0;
    a = (lower - mu)/sigma;
    b = (upper - mu)/sigma;
    
    // First scenario
    if( (a == R_NegInf) || (b == R_PosInf))
    {
        if(a == R_NegInf)
        {
            change = 1;
            a = -b;
            b = R_PosInf;
        }
        
        // The two possibilities for this scenario
        if(a <= 0.45) z = norm_rs(a, b);
        else z = exp_rs(a, b);
        if(change) z = -z;
    }
    // Second scenario
    else if((a * b) <= 0.0)
    {
        // The two possibilities for this scenario
        if((R::dnorm(a, 0.0, 1.0,1.0) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1.0) <= logt1))
        {
            z = norm_rs(a, b);
        }
        else z = unif_rs(a,b);
    }
    // Third scenario
    else
    {
        if(b < 0)
        {
            tmp = b; b = -a; a = -tmp; change = 1;
        }
        
        lograt = R::dnorm(a, 0.0, 1.0, 1.0) - R::dnorm(b, 0.0, 1.0, 1.0);
        if(lograt <= logt2) z = unif_rs(a,b);
        else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
        else z = exp_rs(a,b);
        if(change) z = -z;
    }
    double output;
    output = sigma*z + mu;
    return (output);
}

// [[Rcpp::export]]
arma::mat BAKRGibbsProbit(mat U, vec y, mat Lambda, int iter = 5000, int burn = 1000, double nu_alpha = 2, double phi_alpha = 1){
    
    int i,j;
    double q = U.n_cols;
    double n = U.n_rows;
    mat alphamat = zeros<mat>(iter-burn, q);
    vec S = zeros<vec>(n);
    
    // The rest of the code follows the R version
    rowvec alpha = zeros<rowvec>(q);
    double sigma_alpha = 0.5;
    mat I_q(q,q); I_q.eye();
    
    for (i=0; i<iter; i++) {
        vec ms = U*trans(alpha);
        
        for(j=0; j < y.n_elem; j++){
            if(y(j)==0){
                S(j) = rnorm_trunc(ms(j),1,R_NegInf,0);
            }
            else{
                S(j) = rnorm_trunc(ms(j),1,0,R_PosInf);
            }
        }
        
        alpha = mvrnormArma(1,sigma_alpha*inv(inv(Lambda)+sigma_alpha*I_q)*trans(U)*S, sigma_alpha*inv(inv(Lambda)+sigma_alpha*I_q));
        
        sigma_alpha = rScaledInvChiSq(1,nu_alpha + q,pow(nu_alpha + q,-1)*(nu_alpha*phi_alpha+as_scalar(alpha*inv(Lambda)*trans(alpha))));
        
        if(i>=burn){
            alphamat.row(i-burn) = alpha;
        }
    }
    return alphamat;
}

// [[Rcpp::export]]
NumericVector BAKRProbitPredict(NumericMatrix X_new, NumericVector beta){
  NumericVector pred = Rcpp::pnorm(transpose(X_new)*beta,0.0,1.0,1,0);
    return pred;
}

// [[Rcpp::export]]
vec lfsr(mat beta_mat){
    int i;
    double M = beta_mat.n_rows;
    double P = beta_mat.n_cols;
    vec LFSR = zeros<vec>(P);
    
    for (i=0; i<P; i++) {
        double m = mean(beta_mat.col(i));
        if(m > 0){
            LFSR(i) = sum(beta_mat.col(i)<0)/M;
        }
        else{
            LFSR(i) = sum(beta_mat.col(i)>0)/M;
        }
    }
    return LFSR;
}

// [[Rcpp::export]]
NumericMatrix GetPPAAs(NumericMatrix betamat, double sigval){
    int i,j;
    int p = betamat.ncol();
    NumericMatrix PPAAmat(betamat.nrow(),p);
    vec qval = zeros<vec>(p);
    
    for(j=0; j<p; j++){
        double zz = j+1;
        qval(j) = sigval*(zz/p);
    }
        
    for(i=0; i<betamat.nrow(); i++){
        NumericVector PIP(p);
        vec pval = zeros(p);
        vec pvalsort = zeros(p);
        NumericVector beta = betamat.row(i);
        double sigmahat = mad_rcpp(beta)/0.674;
        pval = 2*(1-Rcpp::pnorm(abs(beta)/sigmahat,0.0,1.0,1,0));
        pvalsort = sort(pval);
        if(any(pvalsort<qval)){
            double pvalmax = max(pvalsort(find(pvalsort<qval)));
            double t = sigmahat*R::qnorm(1-pvalmax/2,0.0,1.0,1,0);
            PIP[abs(beta) > t] = 1;
            PPAAmat.row(i) = PIP;
        }
    }
    return PPAAmat;
}