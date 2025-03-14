#' @title Quasi-Likelihood Estimation for a Spatial Negative Binomial Model
#' @description Iterative estimation of gamma, beta, rho1, and rho2 in a spatial negative binomial model.
#' @param y A numeric vector of observed counts.
#' @param z A matrix of covariates for `p_i`.
#' @param x A matrix of covariates for `lambda_i`.
#' @param W A spatial weight matrix (n x n).
#' @param gamma Initial coefficient vector for `p_i`.
#' @param beta Initial coefficient vector for `lambda_i`.
#' @param rho1 Initial spatial parameter for `p_i`.
#' @param rho2 Initial spatial parameter for `lambda_i`.
#' @param t Initial dispersion parameter of the negative binomial distribution.
#' @param max_iter Maximum number of iterations for convergence.
#' @param tol Convergence tolerance.
#' @return A list with estimated parameters `gamma`, `beta`, `rho1`, `rho2`, and `t`.
#' @examples
#' set.seed(123)
#' n <- 100
#' y <- rpois(n, lambda = 5)  
#' z <- matrix(runif(n * 3), n, 3)
#' x <- matrix(runif(n * 3), n, 3)
#' W <- matrix(runif(n * n), n, n)
#' diag(W) <- 0
#' W <- W/rowSums(W)
#' gamma <- rep(0, 3)
#' beta <- rep(0, 3)
#' rho1 <- 0.1
#' rho2 <- 0.1
#' t <- 1
#' result <- estimate_parameters(y, z, x, W, gamma, beta, rho1, rho2, t)
#' print(result)
#'
#' @export

estimate_parameters <- function(y, z, x, W, gamma, beta, rho1, rho2,
                                t, max_iter = 100, tol = 1e-6) {
  n <- length(y)
  
  # Initialize matrices
  A1 <- diag(n) - rho1 * W
  A2 <- diag(n) - rho2 * W
  
  for (iter in 1:max_iter) {
    # Update p and lambda
    p <- 1 / (1 + exp(-solve(A1) %*% z %*% gamma))
    lambda <- exp(solve(A2) %*% x %*% beta)
    
    # Compute u
    u <- (1 + ((1 - p) * (1 / (1 + t * lambda))^(1/t)) / p)^(-1) * (y == 0)
    
    # Compute gradients
    dp_dgamma <- t(z) %*% matrix((p * (1 - p)), ncol=1)
    dlambda_dbeta <- x * lambda
    
    # Compute V1 and V2
    V1 <- diag(sqrt(p * (1 - p))) %*% solve(t(A1) %*% A1) %*% diag(sqrt(p * (1 - p)))
    V2 <- diag(sqrt(lambda * (1 + t * lambda))) %*% solve(t(A2) %*% A2) %*% diag(sqrt(lambda * (1 + t * lambda)))
    
    # Update gamma
    gamma_update <- solve(t(dp_dgamma) %*% solve(V1) %*% dp_dgamma) %*% (t(dp_dgamma) %*% solve(V1) %*% (u - p))
    gamma <- gamma + gamma_update
    
    # Update beta
    beta_update <- solve(t(dlambda_dbeta) %*% solve(V2) %*% diag(1 - u) %*% dlambda_dbeta) %*% 
      (t(dlambda_dbeta) %*% solve(V2) %*% diag(1 - u) %*% (y - lambda))
    beta <- beta + beta_update
    
    # Update t
    t <- sum(lambda^2 * (1 - u)^2 * (y^2 - 2 * y * lambda - lambda^2 - lambda)) / sum((1 - u) * lambda^4)
    
    # Update rho1 and rho2 simultaneously using Newton-Raphson
    likelihood_grad <- function(rho) {
      rho1 <- rho[1]
      rho2 <- rho[2]
      A1 <- diag(n) - rho1 * W
      A2 <- diag(n) - rho2 * W
      p <- 1 / (1 + exp(-solve(A1) %*% z %*% gamma))
      lambda <- exp(solve(A2) %*% x %*% beta)
      log_likelihood <- sum(log(p + (1 - p) * dnbinom(y, size = 1/t, mu = lambda)))
      return(-log_likelihood)
    }
    
    hessian_func <- function(rho) {
      rho1 <- rho[1]
      rho2 <- rho[2]
      A1_inv <- solve(diag(n) - rho1 * W)
      A2_inv <- solve(diag(n) - rho2 * W)
      grad1 <- sum(A1_inv %*% z %*% gamma * (p * (1 - p)))
      grad2 <- sum(A2_inv %*% x %*% beta * lambda)
      H <- matrix(c(grad1, 0, 0, grad2), 2, 2)  # Hessian (approximate)
      return(H)
    }
    
    rho_old <- c(rho1, rho2)
    rho_new <- rho_old - solve(hessian_func(rho_old)) %*% c(likelihood_grad(rho_old))
    
    # Ensure rho remains within (-1,1)
    rho1 <- max(min(rho_new[1], 1), -1)
    rho2 <- max(min(rho_new[2], 1), -1)
    
    # Convergence check
    if (max(abs(gamma_update), abs(beta_update), abs(rho_new - rho_old)) < tol) {
      break
    }
  }
  
  return(list(gamma = gamma, beta = beta, rho1 = rho1, rho2 = rho2, t = t))
}
