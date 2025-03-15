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


set.seed(123)
n <- 100
y <- rpois(n, lambda = 2)  
z <- matrix(runif(n * 3), n, 3)
x <- matrix(runif(n * 3), n, 3)
W <- matrix(runif(n * n), n, n)
diag(W) <- 0
W <- W/rowSums(W)
result <- estimate_parameters(y, z, x, W, gamma, beta, rho1, rho2, t)
print(result)



estimate_parameters <- function(y, z, x, W, max_iter = 100, tol = 1e-6) {
  n <- length(y)
  rho <- c(0, 0)
  rho1 <- rho[1]
  rho2 <- rho[2]
  mod0 <- ZIM::zim.fit(y,x,z)
  tau <- 1
  beta <- mod0$para[1:(ncol(x))]
  gamma <- mod0$para[(ncol(x)+1):(ncol(x)+ncol(z))]
  # Initialize matrices
  A1 <- diag(n) - rho1 * W
  A2 <- diag(n) - rho2 * W
  
  for (iter in 1:max_iter) {
    # Update p and lambda
    p <- c(1 / (1 + exp(-solve(A1) %*% z %*% gamma)))
    lambda <- c(exp(solve(A2) %*% x %*% beta))
    
    # Compute u
    u <- (1 + ((1 - p) * (1 / (1 + tau * lambda))^(1/tau)) / p)^(-1) * (y == 0)
    
    # Compute gradients
    dp_dgamma <- z*matrix(p * (1 - p), nrow=n, ncol=ncol(z), byrow=FALSE)
    dlambda_dbeta <- x * matrix(lambda, nrow=n, ncol=ncol(x), byrow=FALSE)
    
    # Compute V1 and V2
    V1 <- diag(sqrt(p * (1 - p))) %*% solve(t(A1) %*% A1) %*% diag(sqrt(p * (1 - p)))
    V2 <- diag(sqrt(lambda * (1 + tau * lambda))) %*% solve(t(A2) %*% A2) %*%
      diag(sqrt(lambda * (1 + tau * lambda)))
    
    # Update gamma
    gamma_update <- solve(t(dp_dgamma) %*% solve(V1) %*% dp_dgamma) %*%
      (t(dp_dgamma) %*% solve(V1) %*% (u - p))
    gamma <- gamma + gamma_update
    
    # Update beta
    beta_update <- solve(t(dlambda_dbeta) %*% solve(V2) %*% diag(1 - u) %*% dlambda_dbeta) %*% 
      (t(dlambda_dbeta) %*% solve(V2) %*% diag(1 - u) %*% (y - lambda))
    beta <- beta + beta_update
    
    # Update t
    tau <- sum(lambda^2 * (1 - u)^2 * ((y- lambda)^2))/ sum((1 - u)^2 * lambda^4)
    
    # Update rho1 and rho2 simultaneously using Newton-Raphson
    likelihood_grad <- function(rho) {
      rho1 <- rho[1]
      rho2 <- rho[2]
      A1 <- diag(n) - rho1 * W
      A2 <- diag(n) - rho2 * W
      p <- 1 / (1 + exp(-solve(A1) %*% z %*% gamma))
      lambda <- exp(solve(A2) %*% x %*% beta)
      l0 <- p + (1 - p) * dnbinom(y, size = 1/tau, mu = lambda)
      l1 <-  (1 - p) * dnbinom(y, size = 1/tau, mu = lambda)
      log_likelihood <- sum(log(l0*(y==0)+l1*(y>0)))
      return(-log_likelihood)
    }
    
    parTemp <- optim(par=c(0,0), fn=likelihood_grad, method= "L-BFGS-B", lower=c(-0.99, -0.99),
          upper=c(0.99, 0.99), hessian=TRUE)
    
    rhoNew <- parTemp$par

    # Convergence check
    if (max(abs(gamma_update), abs(beta_update), abs(rhoNew - rho)) < tol) {
      break
    }
    rho =rhoNew
    print(gamma)
  }
  
  return(list(gamma = gamma, beta = beta, rho = rhoNew, 
              vcovrho=solve(parTemp$hessian), tau = tau))
}

estimate_parameters(y,x,z,W)
