#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


arma::mat get_logit_matrix(const arma::mat& x) {
  return 1 / (1 + exp(-x));
}


// Function for double input
double get_lambdastar_double(double x, double thetas, const arma::vec &lambdas) {
  double pstar0 = (1 - thetas) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = thetas * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  double lambdastar = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
  return lambdastar;
}


// Function for matrix input
arma::mat get_lambdastar_matrix(const arma::mat &x, const arma::vec &thetas, const arma::vec &lambdas) {

  arma::mat lambdastar_matrix(x.n_rows, x.n_cols);

  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword k = 0; k < x.n_cols; ++k) {
      double value = x(i, k);
      double theta = thetas[k];
      double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(value));
      double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(value));
      double pstar = pstar1 / (pstar0 + pstar1);
      lambdastar_matrix(i, k) = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
    }
  }

  return lambdastar_matrix;
}


double g(double x, double theta, double eta, const arma::vec &lambdas){
  double lambdastar=get_lambdastar_double(x,theta,lambdas);
  double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  return pow((lambdastar-lambdas[1]),2)+2/eta*log(pstar);
}


double get_delta(double theta, double eta, const arma::vec &lambdas){
  if (lambdas[0] == lambdas[1]){
    return eta * lambdas[1];
  } else {
    if (g(0, theta, eta, lambdas) > 0){
      double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * 0);
      double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * 0);
      double pstar = pstar1 / (pstar0 + pstar1);
      return sqrt(2 * eta * log(1/pstar)) + eta * lambdas[1];
    }  else {
      return eta * get_lambdastar_double(0,theta,lambdas);
    }
  }
}


double soft_thresholding(double x, double lambdastar, double eta, double delta) {
  double s = 0;
  if (x > 0) s = 1;
  else if (x < 0) s = -1;
  if (fabs(x) <= delta) {
    return(0);
  } else {
    double temp;
    temp = fabs(x) - eta * lambdastar;
    if (temp > 0) {
      return(temp * s);
    } else {
      return(0);
    }
  }
}


void update_A_B(arma::mat& A, arma::mat& B,
                const arma::mat& Y, double eta, double xi, const arma::vec& mu,
                const arma::mat& A_momentum, const arma::mat& B_momentum,
                const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                const arma::vec& tilde_thetas, const arma::vec& thetas) {

  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = A_momentum.n_cols;
  arma::vec one_vec = arma::ones(J);

  // Compute shared terms using A_momentum and B_momentum
  arma::mat W = (1 + xi * Y - Y) / (1 + arma::exp(- mu * one_vec.t() - A_momentum * B_momentum.t()));
  arma::mat dA = - xi * Y * B_momentum + W * B_momentum;
  arma::mat dB = - xi * Y.t() * A_momentum + W.t() * A_momentum;

  // Update A using A_momentum
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < I; ++i) {
      double lambdastar = get_lambdastar_double(A_momentum(i, k), tilde_thetas(k), tilde_lambdas);
      double delta = get_delta(tilde_thetas(k), eta, tilde_lambdas);
      A(i, k) = soft_thresholding(A_momentum(i, k) - eta * dA(i, k), lambdastar, eta, delta);
    }
  }

  // Update B using B_momentum
  for (int k = 0; k < K; ++k) {
    for (int j = 0; j < J; ++j) {
      double lambdastar = get_lambdastar_double(B_momentum(j, k), thetas(k), lambdas);
      double delta = get_delta(thetas(k), eta, lambdas);
      B(j, k) = soft_thresholding(B_momentum(j, k) - eta * dB(j, k), lambdastar, eta, delta);
    }
  }

}


void update_mu(arma::vec& mu,
               const arma::mat& A, const arma::mat& B,
               const arma::mat& Y, double eta, double xi){

  int I = Y.n_rows;
  int J = Y.n_cols;
  arma::vec one_vec = arma::ones(J);
  arma::mat p = 1 / (1 + arma::exp(-(mu * one_vec.t() + A * B.t())));  // Compute p_ij

  for (int i = 0; i < I; ++i) {
    double sum_xi_yij = 0.0;
    double sum_adjusted_pij = 0.0;
    double sum_xi_yij_plus_other = 0.0;

    for (int j = 0; j < J; ++j) {
      sum_xi_yij += xi * Y(i, j);
      sum_adjusted_pij += (xi * Y(i, j) + 1 - Y(i, j)) * p(i, j);
      sum_xi_yij_plus_other += xi * Y(i, j) + 1 - Y(i, j);
    }

    // Update mu_i
    mu[i] += 4 * (1 / sum_xi_yij_plus_other) * (sum_xi_yij - sum_adjusted_pij);
  }
}


void rescale_A_B(arma::mat& A, arma::mat& B) {

  int K = A.n_cols;
  int N = A.n_rows;
  int G = B.n_rows;
  double A_norm = 0;
  double B_norm = 0;
  arma::rowvec d = arma::ones<arma::rowvec>(K);

  for(int k = 0; k < K; k++) {
    for (int i = 0; i < N; i++) {
      A_norm += fabs(A(i, k));
    }
    for (int j = 0; j < G; j++) {
      B_norm += fabs(B(j, k));
    }
    if ((A_norm > 0) && (B_norm > 0)) {
      d(k) = pow(A_norm / B_norm, 0.5);
    }
    A_norm = 0;
    B_norm = 0;
  }
  A.each_row() /= d;
  B.each_row() %= d;

}


List get_thetas(const arma::mat& mat, double alpha, double beta, double tol = 0) {
  int K = mat.n_cols;
  int N = mat.n_rows;
  arma::vec thetas(K, arma::fill::zeros);
  arma::vec counts(K, arma::fill::zeros);

  for (int k = 0; k < K; ++k) {
    int q = sum(abs(mat.col(k)) > tol);
    counts[k] = q;
    thetas[k] = (alpha + q) / (alpha + beta + N);
  }

  return List::create(Named("thetas") = thetas, Named("counts") = counts);
}


List get_logLikelihood(const arma::mat& Y, double xi, const arma::vec mu,
                             const arma::mat& A, const arma::mat& B,
                             const arma::vec& tilde_thetas, const arma::vec& thetas,
                             const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                             double tilde_alpha, double tilde_beta, double alpha, double beta) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  arma::vec one_vec = arma::ones(J);

  // Compute the latent matrix M
  arma::mat M = mu * one_vec.t() + A * B.t();

  // Likelihood component
  arma::mat logLikelihood_matrix = xi * Y % M - (xi * Y + 1 - Y) % arma::log(1 + arma::exp(M));
  double logLikelihood_sum = arma::accu(logLikelihood_matrix);

  // Initialize accumulators
  double lambdastar_sum_A = 0.0;
  double lambdastar_sum_B = 0.0;

  // Compute lambdastar sums for A
  for (int i = 0; i < I; ++i) {
    for (int k = 0; k < A.n_cols; ++k) {
      lambdastar_sum_A += get_lambdastar_double(A(i, k), tilde_thetas(k), tilde_lambdas) * abs(A(i, k));
    }
  }

  // Compute lambdastar sums for B
  for (int j = 0; j < J; ++j) {
    for (int k = 0; k < B.n_cols; ++k) {
      lambdastar_sum_B += get_lambdastar_double(B(j, k), thetas(k), lambdas) * abs(B(j, k));
    }
  }

  // Combine all components
  double total_logLikelihood = logLikelihood_sum - lambdastar_sum_A - lambdastar_sum_B;

  // Additional adjustment term
  double BIC = -2 * total_logLikelihood + (std::log(I*J) * (1 + accu(A != 0) + accu(B != 0)));

  return List::create(Named("logLikelihood") = total_logLikelihood, Named("BIC") = BIC);
}

//' Main iterative updates for BiSSLB model using RcppArmadillo.
//'
//' This function performs the core iterative updates for the Binary Spike-and-Slab
//' Lasso Biclustering (BiSSLB) model. It alternates between updating the latent
//' row profile matrix \eqn{A}, column profile matrix \eqn{B}, and sparsity parameters
//' (\eqn{\theta}, \eqn{\tilde{\theta}}) under spike-and-slab priors. It supports
//' optional Indian Buffet Process (IBP) reordering and soft-thresholding.
//'
//' @param mu A numeric vector representing the baseline signal (location parameter).
//' @param A A numeric matrix representing initial row profiles.
//' @param B A numeric matrix representing initial column profiles.
//' @param tilde_thetas A numeric vector of initial sparsity parameters for \eqn{A}.
//' @param thetas A numeric vector of initial sparsity parameters for \eqn{B}.
//' @param Y A binary observation matrix to be approximated.
//' @param xi A positive scalar controlling the pseudo-likelihood scale.
//' @param eta A small positive scalar controlling soft-thresholding aggressiveness.
//' @param tilde_lambdas Spike-and-slab prior lambdas for \eqn{A}.
//' @param lambdas Spike-and-slab prior lambdas for \eqn{B}.
//' @param tilde_alpha Shape parameter of Beta prior on \eqn{\tilde{\theta}} (for A).
//' @param tilde_beta Rate parameter of Beta prior on \eqn{\tilde{\theta}} (for A).
//' @param alpha Shape parameter of Beta prior on \eqn{\theta} (for B).
//' @param beta Rate parameter of Beta prior on \eqn{\theta} (for B).
//' @param tol Convergence tolerance for relative change in log-likelihood.
//' @param max_iter Maximum number of EM-like iterations.
//' @param IBP Integer flag (1 = enable, 0 = disable) for IBP-style column reordering.
//' @param I Number of rows in \eqn{A} (features).
//' @param J Number of rows in \eqn{B} (samples).
//'
//' @return A list containing:
//' \describe{
//'   \item{A}{Final row profile matrix \eqn{A}.}
//'   \item{B}{Final column profile matrix \eqn{B}.}
//'   \item{mu}{Updated location parameter vector.}
//'   \item{tilde_thetas}{Updated sparsity vector for \eqn{A}.}
//'   \item{thetas}{Updated sparsity vector for \eqn{B}.}
//'   \item{delta_logLikelihood_save}{Vector of log-likelihood changes per iteration.}
//'   \item{logLikelihood_save}{Vector of log-likelihood values per iteration.}
//'   \item{BIC}{Final Bayesian Information Criterion (BIC) value.}
//' }
//'
//' @keywords internal
//'
// [[Rcpp::export]]
List main_iterations(arma::vec mu,
                     arma::mat A,
                     arma::mat B,
                     arma::vec tilde_thetas,
                     arma::vec thetas,
                     arma::mat Y,
                     double xi,
                     double eta,
                     arma::vec tilde_lambdas,
                     arma::vec lambdas,
                     double tilde_alpha,
                     double tilde_beta,
                     double alpha,
                     double beta,
                     double tol,
                     int max_iter,
                     int IBP,
                     int I,
                     int J) {

  // Prepare work
  arma::mat A_lag = A;
  arma::mat B_lag = B;
  arma::mat A_momentum = A;
  arma::mat B_momentum = B;

  double logLikelihood_old = get_logLikelihood(Y, xi, mu, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta)["logLikelihood"];
  double delta_logLikelihood_old = 1.0;
  arma::vec delta_logLikelihood_save;
  arma::vec logLikelihood_save;

  int K = thetas.size();

  for (int i = 1; i <= max_iter; ++i) {

    if (i > 1) {
      A_momentum = A + (i - 2.0) / (i + 1.0) * (A - A_lag);
      B_momentum = B + (i - 2.0) / (i + 1.0) * (B - B_lag);

      A_lag = A;  // Store the previous state of A
      B_lag = B;  // Store the previous state of B

      // Call the function to update A and B directly
      update_A_B(A, B, Y, eta, xi, mu, A_momentum, B_momentum, tilde_lambdas, lambdas, tilde_thetas, thetas);
    }

    // Update mu
    update_mu(mu, A, B, Y, eta, xi);

    // Update sparsity
    List A_sparse = get_thetas(A, tilde_alpha, tilde_beta, tol);
    List B_sparse = get_thetas(B, alpha, beta, tol);

    // Extract and sum counts
    arma::vec A_counts = as<arma::vec>(A_sparse["counts"]);
    arma::vec B_counts = as<arma::vec>(B_sparse["counts"]);
    arma::vec counts = A_counts;

    // Extract thetas
    tilde_thetas = as<arma::vec>(A_sparse["thetas"]);
    thetas = as<arma::vec>(B_sparse["thetas"]);

    // Re-order the columns
    arma::uvec re_order;
    if (IBP == 1) {
      re_order = sort_index(counts, "descend");

      A = A.cols(re_order);
      B = B.cols(re_order);
      A_lag = A_lag.cols(re_order);
      B_lag = B_lag.cols(re_order);
      tilde_thetas = tilde_thetas(re_order);
      thetas = thetas(re_order);
    }

    // remove zeros
    if (i % 100 == 0) {
      arma::vec A_zero = arma::zeros<arma::vec>(K);
      arma::vec B_zero = arma::zeros<arma::vec>(K);
      for (int k = 0; k < K; k++) {
        for (int i = 0; i < I; i++) {
          if (A(i, k) == 0) {
            A_zero(k)++;
          }
        }
        for (int j = 0; j < J; j++) {
          if (B(j, k) == 0) {
            B_zero(k)++;
          }
        }
      }

      arma::uvec keep = find(A_zero < I-1 && B_zero < J-1);

      if (keep.n_elem < K && keep.n_elem > 0) {
        A = A.cols(keep);
        B = B.cols(keep);
        A_lag = A_lag.cols(keep);
        B_lag = B_lag.cols(keep);
        tilde_thetas = tilde_thetas(keep);
        thetas = thetas(keep);
        K = A.n_cols;
      }

      if (keep.n_elem == 0) {
        K = 0;
        break;
      }

    }

    // Stop earlier by log likelihood
    if (i > 1) {
      double logLikelihood = get_logLikelihood(Y, xi, mu, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta)["logLikelihood"];

      if (std::isfinite(abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old))) {
        double delta_logLikelihood = abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old);

        logLikelihood_old = logLikelihood;
        delta_logLikelihood_old = delta_logLikelihood;

        delta_logLikelihood_save = arma::join_vert(delta_logLikelihood_save, arma::vec({delta_logLikelihood_old}));
        logLikelihood_save = arma::join_vert(logLikelihood_save, arma::vec({logLikelihood_old}));

        if (delta_logLikelihood < tol) {
          break;
        }
      }
    }

    // Re-scale A and B
    rescale_A_B(A, B);
    // rescale_A_B(A_lag, B_lag);

  }

  double BIC = get_logLikelihood(Y, xi, mu, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta)["BIC"];

  return List::create(
    Named("A") = A,
    Named("B") = B,
    Named("mu") = mu,
    Named("tilde_thetas") = tilde_thetas,
    Named("thetas") = thetas,
    Named("delta_logLikelihood_save") = delta_logLikelihood_save,
    Named("logLikelihood_save") = logLikelihood_save,
    Named("BIC") = BIC
  );
}


