#' SARARgamlss: Spatial Autoregressive Generalized Additive Model for Location Scale (GAMLSS)
#'
#' This function estimates a Spatial Autoregressive Generalized Additive Model for Location Scale 
#' (SARARgamlss) using GAMLSS. The model includes both spatial dependencies and the possibility of 
#' non-parametric terms in the formulas for the mean and variance. The function supports SAR, SARAR, 
#' and SEM model types and performs the estimation through an iterative process that updates spatial 
#' dependence parameters. The variance of the spatial parameters \(\hat{\rho}\) and \(\hat{\lambda}\) 
#' is estimated using the inverse of the Hessian matrix from the optimization.
#'
#' @param formula A formula specifying the mean structure of the model (response ~ explanatory variables).
#' @param sigma.formula A formula specifying the variance structure of the model (default: ~1).
#' @param W1 A spatial weights matrix for the SAR term (default: identity matrix).
#' @param W2 A spatial weights matrix for the SARAR term (default: identity matrix).
#' @param data A data.frame containing the variables used in the model.
#' @param tol Convergence tolerance (default: 1E-4).
#' @param maxiter Maximum number of iterations for optimization (default: 20).
#' @param type The type of spatial model to fit: one of "SAR", "SARAR", or "SEM".
#' @param weights Optional weights for the observations (default: NULL).
#' @param ... Additional arguments passed to `gamlss`.
#' 
#' @return A fitted GAMLSS model object with spatial autoregressive terms. The model object also includes
#' the variance of the spatial parameters \(\hat{\rho}\) and \(\hat{\lambda}\).
#' @export
#' @import gamlss splines
SARARgamlss <- function(formula = formula(data), sigma.formula = ~1,
                        W1 = diag(0, nrow(data)), W2 = diag(0, nrow(data)),
                        data, tol = 1E-4, maxiter = 20,
                        type = c("SAR", "SARAR", "SEM"),
                        weights = NULL, ...) {
  
  # Check for model type and handle accordingly
  if (type == "SAR" & sum(W1) == 0) {
    print("The SAR Model contains W1=0, it is a usual non-spatial GAMLSS")
  }
  if (type == "SAR") {
    W2 = 0 * W2
  }
  if (type == "SEM") {
    W1 = 0 * W1
  }
  if (type == "SEM" & sum(W2) == 0) {
    print("The SEM Model contains W2=0, it is a usual non-spatial GAMLSS")
  }
  if (type == "SARAR" & sum(W2) == 0 & sum(W1) == 0) {
    print("The SARAR Model contains W1=0 and W2=0, it is a usual non-spatial GAMLSS")
  }
  if (type == "SARAR" & sum(W2) == 0) {
    print("The SARAR Model assumes that W1 = W2")
  }
  
  # Initialize variables
  rho0 <- 0
  lambda <- 0
  n <- nrow(data)
  
  # Fit the initial GAMLSS model for mean (mu) and variance (sigma)
  m0 <- gamlss::gamlss(formula = formula, sigma.formula = sigma.formula, 
                       data = data, family = NO(), ...)
  
  Y <- matrix(m0$y, ncol = 1)
  namesX <- colnames(model.matrix(m0, what = "mu"))
  namesZ <- colnames(model.matrix(m0, what = "sigma"))
  
  # Initial variance estimation (sigma)
  var0 <- predict(m0, what = "sigma", type = "response")^2
  Xbeta <- predict(m0, what = "mu", type = "response")
  
  # Define the log-likelihood function
  loglik <- function(rholam, W1, W2, Xbeta, Y, var0) {
    AA <- diag(n) - rholam[1] * W1
    BB <- diag(n) - rholam[2] * W2
    VV <- BB %*% (AA %*% Y - Xbeta) / sqrt(var0)
    loglik <- -0.5 * sum(log(var0)) + log(det(AA)) + log(det(BB)) - 
      0.5 * sum(VV^2)
    return(-loglik)
  }
  
  # Optimization step to estimate spatial parameters (rho, lambda) with hessian
  p0 <- optim(par = c(0, 0), fn = loglik, method = "L-BFGS-B", W1 = W1, W2 = W2, 
              Xbeta = Xbeta, Y = Y, var0 = var0, 
              lower = c(-0.999, -0.999), upper = c(0.999, 0.99), hessian = TRUE)
  
  ## Hessian <- p0$hessian  # Extract Hessian matrix
  p0 <- p0$par
  
  tolTemp <- 1
  iter <- 1
  
  # Iteratively update spatial parameters and GAMLSS model
  while (tolTemp > tol & iter < maxiter) {
    p1 <- p0
    AA <- diag(n) - p1[1] * W1
    BB <- diag(n) - p1[2] * W2
    Ytemp <- as.matrix(BB %*% AA %*% Y)
    Xtemp <- as.matrix(BB %*% model.matrix(m0, what = "mu"))
    colnames(Xtemp) <- namesX
    Ztemp <- model.matrix(m0, what = "sigma")
    colnames(Ztemp) <- namesZ
    
    # Fit updated GAMLSS model with transformed data (dependent and independent)
    m1 <- gamlss::gamlss(Ytemp ~ Xtemp - 1, ~Ztemp - 1)
    var1 <- predict(m1, what = "sigma", type = "response")^2
    Xbeta <- predict(m1, what = "mu", type = "response")
    
    # Optimize again with updated variance
    p0 <- optim(par = c(0, 0), fn = loglik, method = "L-BFGS-B", W1 = W1, W2 = W2, 
                Xbeta = Xbeta, Y = Y, var0 = var1, 
                lower = c(-0.999, -0.999), upper = c(0.999, 0.99), hessian = TRUE)
    
    Hessian <- p0$hessian  # Extract Hessian matrix
    p0 <- p0$par
    
    var_cov_matrix <- -solve(Hessian)  # Inverse Hessian for variance-covariance
    
    
    
    # Update variance-covariance matrix
    var_rho <- var_cov_matrix[1, 1]
    var_lambda <- var_cov_matrix[2, 2]
    
    model_variances <- list(var_rho = var_rho, var_lambda = var_lambda)
    # Store variances in model object
    model_variances$var_rho <- var_rho
    model_variances$var_lambda <- var_lambda
    
    tolTemp <- sum(abs(p1 - p0))
    iter <- iter + 1
  }
  
  # Store updated model parameters
  namesTemp <- c("G.deviance", "residuals", "mu.fv", "mu.lp", "mu.wv", "mu.wt", 
                 "mu.qr", "mu.coefficients", "sigma.fv", "sigma.wv", "sigma.wt", 
                 "sigma.coefficients", "P.deviance", "aic", "sbc")
  
  for (names in namesTemp) {
    m0[[names]] <- m1[[names]]
  }
  
  # Return the final model object with spatial parameters and variance estimates
  if (type == "SARAR") {
    m0$rho <- p0[1]
    m0$lambda <- p0[2]
  } else if (type == "SAR") {
    m0$rho <- p0[1]
  } else {
    m0$lambda <- p0[2]
  }
  
  m0$variances <- model_variances  # Add variance estimates to the model
  class(m0) <- "SARARgamlss"
  return(m0)
}

