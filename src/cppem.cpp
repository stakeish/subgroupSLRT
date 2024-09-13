#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List cppem_normal(arma::vec y, arma::mat x, arma::mat z, arma::vec d,
           double pi, arma::mat alpha_init, arma::vec beta1_init, arma::vec beta2_init,
           double sigma_init, int ninits, int maxit, double tol) {

  // Initialization
  // Rcout << "Initialization" << endl;
  int q1 = x.n_cols;
  int n = x.n_rows;
  arma::vec ll_set(ninits), beta1_set(ninits), beta2_set(ninits), sigma_set(ninits);
  arma::mat alpha_set(q1, ninits);

  // Rcout << "Initialization" << endl;

  arma::vec alpha_k(q1), coef_k(q1+2);
  double beta1_k, beta2_k, sigma_k, ll, ll_old, diff;

  // Rcout << "Initialization" << endl;

  arma::vec mean_1(n), mean_2(n), x_alpha(n), w1(n), w2(n), lik_vec(n);
  arma::mat wmat(n, 2), mat_equation_alpha(q1, q1+2), mat_equation_beta1(1, q1+2), mat_equation_beta2(1, q1+2);
  arma::mat mat_equation(q1+2, q1+2);
  arma::vec vec_equation_alpha(q1), vec_equation_beta1(1), vec_equation_beta2(1), vec_equation(q1+2);

  // Rcout << "Initialization" << endl;

  arma::mat xx = x.t()*x;

  for(int i1 = 0; i1 < ninits; i1++){

    // Rcout << "Initialization" << endl;

    alpha_k = alpha_init.col(i1);
    beta1_k = beta1_init(i1);
    beta2_k = beta2_init(i1);
    sigma_k = sigma_init;
    diff = 1;
    ll_old = -arma::datum::inf;
    x_alpha = x*alpha_k;

    // Rcout << "Initialization" << endl;

    for(int i2 = 0; i2 < maxit; i2++){

      // E-step

      // Rcout << "Initialization" << endl;

      mean_1 = y - x_alpha - d*beta1_k;
      mean_2 = y - x_alpha - d*beta2_k;
      w1 = pi * arma::normpdf(mean_1, 0, sigma_k);
      w2 = (1 - pi) * arma::normpdf(mean_2, 0, sigma_k);
      wmat = join_horiz(w1, w2);
      lik_vec = sum(wmat, 1);
      wmat = wmat.each_col() / lik_vec;

      // Rcout << "Initialization" << endl;

      ll = sum(log(lik_vec));

      diff = ll - ll_old;
      ll_old = ll;

      //if(diff < 0){
      //  Rcout << "Error: Decrease in log-likelihood" << endl;
      //  break;
      //}

      if(diff < tol){
        break;
      }

      // Rcout << "Initialization" << endl;

      // M-step
      mat_equation_alpha = join_horiz(xx, x.t() * (d % wmat.col(0)), x.t() * (d % wmat.col(1)));
      // Rcout << "Initialization" << endl;
      mat_equation_beta1 = join_horiz((x.t() * (d % wmat.col(0))).t(), d.t() * wmat.col(0), arma::zeros(1, 1));
      // Rcout << "Initialization" << endl;
      mat_equation_beta2 = join_horiz((x.t() * (d % wmat.col(1))).t(), arma::zeros(1, 1), d.t() * wmat.col(1));
      // Rcout << "Initialization" << endl;
      mat_equation = join_vert(mat_equation_alpha, mat_equation_beta1, mat_equation_beta2);

      // Rcout << "Initialization" << endl;

      vec_equation_alpha = x.t() * y;
      vec_equation_beta1 = y.t() * (d % wmat.col(0));
      vec_equation_beta2 = y.t() * (d % wmat.col(1));
      vec_equation = join_vert(vec_equation_alpha, vec_equation_beta1, vec_equation_beta2);

      // Rcout << "Initialization" << endl;

      coef_k = solve(mat_equation, vec_equation);
      alpha_k = coef_k.subvec(0, q1-1);
      beta1_k = coef_k(q1);
      beta2_k = coef_k(q1+1);

      // Rcout << "Initialization" << endl;

      x_alpha = x*alpha_k;
      sigma_k = sqrt(mean(wmat.col(0) % pow(y - x_alpha - d*beta1_k, 2)) + mean(wmat.col(1) % pow(y - x_alpha - d*beta2_k, 2)));
    }

    // Rcout << "Initialization" << endl;

    alpha_set.col(i1) = alpha_k;
    beta1_set(i1) = beta1_k;
    beta2_set(i1) = beta2_k;
    sigma_set(i1) = sigma_k;
    ll_set(i1) = ll;
  }

  // Rcout << "Initialization" << endl;

  return List::create(Named("alpha") = alpha_set,
                      Named("beta1") = beta1_set,
                      Named("beta2") = beta2_set,
                      Named("sigma") = sigma_set,
                      Named("ll") = ll_set);

}

// [[Rcpp::export]]
List EM_gamma_fixed_cpp(int n, arma::mat x, arma::vec d, arma::mat z, arma::vec y, arma::vec gamma, int q1, arma::vec theta_init, double epsilon, int maxit, arma::mat xx, arma::vec xy) {

  arma::vec alpha(q1), x_alpha(n), z_gamma(n), mean_1(n), mean_2(n), log_prop(n), l_vec(n), w1(n), w2(n), w1_d(n), w2_d(n), x_w1_d(q1), x_w2_d(q1), sum_w1_d(1), sum_w2_d(1), vec_equation(q1 + 2), coef(q1 + 2), theta(q1 + 3);
  arma::mat w_mat(n, 2), mat_equation_alpha(q1, q1 + 2), mat_equation_beta(1, q1 + 2), mat_equation_lambda(1, q1 + 2), mat_equation(q1 + 2, q1 + 2);
  double beta, lambda, sigma, loglik, diff;
  double oldloglik = -arma::datum::inf;

  alpha = theta_init.subvec(0, q1 - 1);
  beta = theta_init(q1);
  lambda = theta_init(q1 + 1);
  sigma = theta_init(q1 + 2);

  x_alpha = x * alpha;
  z_gamma = z * gamma;
  NumericVector z_gamma_nv = as<NumericVector>(wrap(z_gamma));
  NumericVector log_prop_nv = Rcpp::plogis(z_gamma_nv);
  log_prop = as<arma::vec>(wrap(log_prop_nv));

  mean_1 = y - x_alpha - d * beta;
  mean_2 = y - x_alpha - d * lambda;

  for(int i1 = 0; i1 < maxit; i1++){

    // E-step
    w_mat.col(0) = log_prop % normpdf(mean_1, 0, sigma);
    w_mat.col(1) = (1 - log_prop) % normpdf(mean_2, 0, sigma);
    l_vec = sum(w_mat, 1);
    w_mat = w_mat.each_col() / l_vec;

    // Exit loop or not
    loglik = sum(log(l_vec));
    diff = loglik - oldloglik;
    oldloglik = loglik;

    if(diff < epsilon || loglik == -arma::datum::inf){
      break;
    }

    // M-step
    w1 = w_mat.col(0);
    w2 = w_mat.col(1);
    w1_d = w1 % d;
    w2_d = w2 % d;
    x_w1_d = x.t() * w1_d;
    x_w2_d = x.t() * w2_d;
    sum_w1_d = sum(w1_d);
    sum_w2_d = sum(w2_d);
    mat_equation_alpha = join_horiz(xx, x_w1_d, x_w2_d);
    mat_equation_beta = join_horiz(x_w1_d.t(), sum_w1_d, arma::zeros(1));
    mat_equation_lambda = join_horiz(x_w2_d.t(), arma::zeros(1), sum_w2_d);
    mat_equation = join_vert(mat_equation_alpha, mat_equation_beta, mat_equation_lambda);

    vec_equation.subvec(0, q1 - 1) = xy;
    vec_equation(q1) = sum(y % w1_d);
    vec_equation(q1+1) = sum(y % w2_d);


    try {
      coef = solve(mat_equation, vec_equation, arma::solve_opts::no_approx);
      //double rank_mat = rank(mat_equation_outcome);
      //Rcout << rank_mat << std::endl;
    } catch(...) {
      break;
    }

    alpha = coef.subvec(0, q1 - 1);
    beta = coef(q1);
    lambda = coef(q1 + 1);
    x_alpha = x * alpha;

    mean_1 = y - x_alpha - d * beta;
    mean_2 = y - x_alpha - d * lambda;

    sigma = sqrt(mean(w1 % pow(mean_1, 2)) + mean(w2 % pow(mean_2, 2)));

  }

  theta.subvec(0, q1 - 1) = alpha;
  theta(q1) = beta;
  theta(q1 + 1) = lambda;
  theta(q1 + 2) = sigma;

  return List::create(Named("theta") = theta, Named("alpha") = alpha, Named("beta") = beta, Named("lambda") = lambda, Named("sigma") = sigma, Named("loglik") = loglik, Named("diff") = diff);

}
