
#################################
#' Get the logit values.
#'
#' @param x A number or a list.
#' @returns A number or a list.
#' @examples
#' get_logit(0)
get_logit <- function(x) {
  if (is.list(x)) {
    lapply(x, function(z) 1 / (1 + exp(-z)))
  } else {
    1 / (1 + exp(-x))
  }
}



#################################
#' Get the Laplace density
#'
#' @param x A number.
#' @param lambda A number for rate parameter (1/scale).
#' @param mu A number for location parameter (default = 0).
#' @returns A number.
#' @examples
#' get_Laplace(0, 1)
#' get_Laplace(0, 1, 2)
get_Laplace <- function(x, lambda, mu = 0) {
  lambda / 2 * exp(-lambda * abs(x - mu))
}


#################################
#' Bi-cluster a binary matrix using Binary Spike-and-Slab Lasso.
#'
#' This function performs bi-clustering on a binary matrix using a Bayesian
#' spike-and-slab Lasso prior (BiSSLB). The binary matrix is approximated by
#' the product of two latent matrices \eqn{A} and \eqn{B}, representing row-
#' and column-specific profiles, respectively. The method uses spike-and-slab
#' priors to enforce sparsity in the latent structures and iteratively estimates
#' the underlying patterns. This function acts as the R interface for the main
#' C++ engine and includes post-processing such as reordering and thresholding.
#'
#' @param Y A binary matrix to be bi-clustered (rows = features, columns = samples).
#' @param A Optional initialization for the row profile matrix \eqn{A}.
#' @param B Optional initialization for the column profile matrix \eqn{B}.
#' @param xi A positive value representing confidence dispersion (deprecated).
#' @param mu Optional initialization of the baseline signal (location parameter).
#' @param tilde_lambda_0 The "spike" (strong shrinkage) parameter for the prior on matrix \eqn{A}.
#' @param tilde_lambda_1 The "slab" (weak shrinkage) parameter for the prior on matrix \eqn{A}.
#' @param lambda_0 The "spike" (strong shrinkage) parameter for the prior on matrix \eqn{B}.
#' @param lambda_1 The "slab" (weak shrinkage) parameter for the prior on matrix \eqn{B}.
#' @param method A character string indicating how to initialize latent matrices; default is "svd".
#' @param tilde_alpha Shape parameter for the Beta prior on sparsity of \eqn{A}.
#' @param tilde_beta Rate parameter for the Beta prior on sparsity of \eqn{A}.
#' @param alpha Shape parameter for the Beta prior on sparsity of \eqn{B}.
#' @param beta Rate parameter for the Beta prior on sparsity of \eqn{B}.
#' @param K_init Initial number of latent factors (columns of \eqn{A}, \eqn{B}).
#' @param thisSeed An integer seed for reproducibility.
#' @param eta A small positive value controlling soft-thresholding smoothness.
#' @param IBP Logical flag (1 = TRUE, 0 = FALSE) to enable IBP-style feature selection and reordering.
#' @param tol Convergence tolerance for log-likelihood change.
#' @param max_iter Maximum number of EM iterations.
#' @param show_plot Logical; whether to plot convergence diagnostics.
#'
#' @return A list containing estimated matrices \eqn{A}, \eqn{B}, and other inference outputs.
#'
#' @examples
#' Y <- matrix(rbinom(100, 1, 0.5), 10, 10)
#' out <- BiSSLB(Y, K_init = 5)
#'
#' @export
BiSSLB <- function(Y, A = NULL, B = NULL, xi = 1, mu = NULL,
                   tilde_lambda_0 = 5, tilde_lambda_1 = 1,
                   lambda_0 = 5, lambda_1 = 1, method = "svd", # centering = TRUE,
                   tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1,
                   K_init = 50, thisSeed = 123, eta = 0.001, IBP = 1,
                   tol = 1e-10, max_iter = 500, show_plot = FALSE) {
  set.seed(thisSeed) # Ensure reproducibility

  I <- nrow(Y)
  J <- ncol(Y)
  K <- K_init

  Y <- as.matrix(Y)

  tilde_lambda_0 <- tilde_lambda_0 # Matrix A's spike parameter
  tilde_lambda_1 <- tilde_lambda_1 # Matrix A's slab parameter
  lambda_0 <- lambda_0 # Matrix B's spike parameter
  lambda_1 <- lambda_1 # Matrix B's slab parameter
  tilde_lambdas <- c(tilde_lambda_0, tilde_lambda_1)
  lambdas <- c(lambda_0, lambda_1)

  A_in <- A
  B_in <- B

  # Initialize A and B
  if (method == "svd") {
    Y_svd <- tryCatch(
      svd(Y),
      error = function(e) {
        message("SVD failed. Using random initialization.")
        NULL
      }
    )

    if (!is.null(Y_svd)) {
      if (length(Y_svd$d) < K_init) {
        # A <- Y_svd$u %*% diag(sqrt(Y_svd$d))
        # B <- Y_svd$v %*% diag(sqrt(Y_svd$d))
        A <- Y_svd$u %*% diag(Y_svd$d)
        B <- Y_svd$v
      } else {
        # A <- Y_svd$u[, 1:K_init] %*% diag(sqrt(Y_svd$d[1:K_init]))
        # B <- Y_svd$v[, 1:K_init] %*% diag(sqrt(Y_svd$d[1:K_init]))
        A <- Y_svd$u[, 1:K_init] %*% diag(Y_svd$d[1:K_init])
        B <- Y_svd$v[, 1:K_init]
      }
    } else {
      # Fallback to random initialization
      A <- matrix(rexp(I * K_init, rate = 1), nrow = I, ncol = K_init)
      B <- matrix(rexp(J * K_init, rate = 1), nrow = J, ncol = K_init)
    }
  } else {
    # Random initialization if method is not "svd"
    A <- matrix(rexp(I * K_init, rate = 1), nrow = I, ncol = K_init)
    B <- matrix(rexp(J * K_init, rate = 1), nrow = J, ncol = K_init)
  }
  if (!is.null(A_in)) A <- A_in
  if (!is.null(B_in)) B <- B_in

  # Initialize other variables
  if (is.null(mu)) mu <- rowMeans(Y)

  K <- K_init <- min(K_init, ncol(A))
  tilde_nus <- sort(rbeta(K_init, tilde_alpha, tilde_beta), decreasing = TRUE)
  tilde_thetas <- rep(0.5, K_init)
  nus <- sort(rbeta(K_init, alpha, beta), decreasing = TRUE)
  thetas <- rep(0.5, K_init)

  # The main loop: Use the Rcpp implementation
  result <- main_iterations(
    max_iter = max_iter,
    A = A,
    B = B,
    mu = mu,
    Y = Y,
    eta = eta,
    xi = xi,
    tilde_lambdas = tilde_lambdas,
    lambdas = lambdas,
    tilde_thetas = tilde_thetas,
    thetas = thetas,
    tilde_alpha = tilde_alpha,
    tilde_beta = tilde_beta,
    alpha = alpha,
    beta = beta,
    tol = tol,
    IBP = IBP,
    I = I,
    J = J
  )

  # Extract results
  A <- result$A
  B <- result$B
  mu <- result$mu
  tilde_thetas <- result$tilde_thetas
  thetas <- result$thetas
  logLikelihood_save <- result$logLikelihood_save
  delta_logLikelihood_save <- result$delta_logLikelihood_save
  BIC <- result$BIC

  # Prepare output
  out <- list(
    A = A,
    B = B,
    mu = mu,
    tilde_thetas = tilde_thetas,
    thetas = thetas,
    counts = sum(colSums(A != 0) + colSums(B != 0)),
    logLikelihood = logLikelihood_save,
    BIC = BIC
  )

  # Optionally show the plot
  if (show_plot) {
    plot(1:length(logLikelihood_save), logLikelihood_save,
      type = ifelse(length(logLikelihood_save) == 1, "o", "l"),
      xlab = "Iterations", ylab = "Log Likelihood"
    )
  }

  return(out)
}

#################################
#' Run BiSSLB across a ladder of spike parameters for model selection.
#'
#' This function calls \code{\link{BiSSLB}} multiple times using different values
#' of \code{tilde_lambda_0} and \code{lambda_0}, forming a "ladder" of spike
#' parameters. This allows users to evaluate the stability and sparsity of the
#' biclustering solution across a range of regularization strengths, and select
#' the best combination based on performance criteria such as BIC.
#'
#' @param Y A binary matrix to be bi-clustered (rows = features, columns = samples).
#' @param xi A positive value representing confidence dispersion (deprecated).
#' @param A Optional initialization for the row profile matrix \eqn{A}.
#' @param B Optional initialization for the column profile matrix \eqn{B}.
#' @param mu Optional initialization of the baseline signal (location parameter).
#' @param tilde_lambda_0s A numeric vector of spike values for matrix \eqn{A}.
#' @param tilde_lambda_1 The shared slab value for matrix \eqn{A}.
#' @param lambda_0s A numeric vector of spike values for matrix \eqn{B}.
#' @param lambda_1 The shared slab value for matrix \eqn{B}.
#' @param method A character string indicating how to initialize latent matrices; default is "svd".
#' @param tilde_alpha Shape parameter for the Beta prior on sparsity of \eqn{A}.
#' @param tilde_beta Rate parameter for the Beta prior on sparsity of \eqn{A}.
#' @param alpha Shape parameter for the Beta prior on sparsity of \eqn{B}.
#' @param beta Rate parameter for the Beta prior on sparsity of \eqn{B}.
#' @param K_init Initial number of latent factors (columns of \eqn{A}, \eqn{B}).
#' @param thisSeed An integer seed for reproducibility.
#' @param eta A small positive value controlling soft-thresholding smoothness.
#' @param IBP Logical flag (1 = TRUE, 0 = FALSE) to enable IBP-style feature selection and reordering.
#' @param tol Convergence tolerance for log-likelihood change.
#' @param max_iter Maximum number of EM iterations.
#' @param show_plot Logical; whether to plot convergence or BIC across ladders.
#'
#' @return A list containing the results of each BiSSLB run, typically including:
#' \describe{
#'   \item{results}{A list of BiSSLB outputs for each (\code{tilde_lambda_0}, \code{lambda_0}) pair.}
#'   \item{BICs}{A matrix or data frame of corresponding BIC values.}
#' }
#'
#' @examples
#' Y <- matrix(rbinom(100, 1, 0.5), 10, 10)
#' out <- BiSSLB_ladder(Y, tilde_lambda_0s = c(5, 10), lambda_0s = c(5, 10),
#'                      tilde_lambda_1 = 1, lambda_1 = 1)
#'
#' @export
BiSSLB_ladder <- function(Y, xi = 1, A = NULL, B = NULL, mu = NULL,
                          tilde_lambda_0s, tilde_lambda_1,
                          lambda_0s, lambda_1, method = "svd",
                          tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1,
                          K_init = 30, thisSeed = 123, eta = 0.001, IBP = 1,
                          tol = 1e-10, max_iter = 1000, show_plot = FALSE) {
  L <- length(lambda_0s)
  counts <- rep(0, L)
  BIC <- rep(0, L)
  update_tilde_lambda_0 <- 1
  update_lambda_0 <- 1
  tilde_lambda_0 <- tilde_lambda_0s[1]
  lambda_0 <- lambda_0s[1]

  for (l in 1:L) {
    print(glue::glue("Ladder {l}:
    tilde_lambdas = ({tilde_lambda_0},{tilde_lambda_1})
    lambdas = ({lambda_0},{lambda_1}) \n"))

    if (l == 1) {
      A <- A
      B <- B
      mu <- mu
      K_init <- K_init
      this_out <- BiSSLB(
        A = A, B = B, mu = mu,
        Y = Y, xi = xi,
        tilde_lambda_0 = tilde_lambda_0, tilde_lambda_1 = tilde_lambda_1,
        lambda_0 = lambda_0, lambda_1 = lambda_1,
        tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta,
        K_init = K_init, method = method,
        thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP,
        max_iter = floor(max_iter / L), show_plot = show_plot
      )
      counts[1] <- this_out$counts
      BIC[1] <- this_out$BIC
      print(glue::glue("K = {ncol(this_out$A)} and counts = {this_out$counts} \n"))
    } else {
      A <- this_out$A
      B <- this_out$B
      mu <- this_out$mu
      K_init <- ncol(A)

      if (update_tilde_lambda_0 == 1) {
        tilde_lambda_0 <- tilde_lambda_0s[l]
      }
      if (update_lambda_0 == 1) {
        lambda_0 <- lambda_0s[l]
      }

      this_out <- BiSSLB(
        A = A, B = B, mu = mu,
        Y = Y, xi = xi,
        tilde_lambda_0 = tilde_lambda_0, tilde_lambda_1 = tilde_lambda_1,
        lambda_0 = lambda_0, lambda_1 = lambda_1,
        tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta,
        K_init = K_init, method = method,
        thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP,
        max_iter = floor(max_iter / L), show_plot = show_plot
      )
      counts[l] <- this_out$counts
      BIC[l] <- this_out$BIC
      if (counts[l] > counts[l - 1]) {
        tilde_lambda_0 <- tilde_lambda_0s[l - 1]
        update_tilde_lambda_0 <- 0
        lambda_0 <- lambda_0s[l - 1]
        update_lambda_0 <- 0
      }

      print(glue::glue("K = {ncol(this_out$A)} and counts = {this_out$counts} \n"))
    }
  }

  out <- list(
    A = this_out$A,
    B = this_out$B,
    mu = this_out$mu,
    tilde_thetas = this_out$tilde_thetas,
    thetas = this_out$thetas,
    counts = this_out$counts,
    logLikelihood = this_out$logLikelihood,
    BIC = BIC
  )

  return(out)
}
