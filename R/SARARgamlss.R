#' SARARgamlss: Spatial Autoregressive Generalized Additive Model for Location Scale (GAMLSS)
#'
#' This function estimates a Spatial Autoregressive Generalized Additive Model for Location Scale 
#' (SARARgamlss) using GAMLSS. The model includes both spatial dependencies and the possibility of 
#' non-parametric terms in the formulas for the mean and variance. The function supports SAR, SARAR, 
#' and SEM model types and performs the estimation through an iterative process that updates spatial 
#' dependence parameters. The variance of the spatial parameters \eqn{\hat{\rho}} and \eqn{\hat{\lambda}}
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
#' 
#' @return A fitted GAMLSS model object with spatial autoregressive terms. The model object also includes
#' the variance of the spatial parameters \eqn{\hat{\rho}} and \eqn{\hat{\lambda}}
#' @references Toloza-Delgado, J. D., Melo, O. O., & Cruz, N. A.
#'  Joint spatial modeling of mean and non-homogeneous variance combining semiparametric SAR 
#'  and GAMLSS models for hedonic prices. Spatial Statistics, 65, 100864 (2025)
#'  @source https://doi.org/10.1016/j.spasta.2024.100864
#' @examples
#' library(spdep)
#' library(gamlss)
#' data(oldcol)
#' # Create spatial weight matrices W1 and W2
#' W1 <- spdep::nb2mat(COL.nb, style = "W")
#' W2 <- W1  # In this case, assume the same spatial weights for both
#' # Fit a SARARgamlss model
#' result <- SARARgamlss(formula = CRIME ~ INC + cs(HOVAL), 
#' sigma.formula = ~ INC + pb(HOVAL), W1 = W1, W2 = W2,data = COL.OLD, 
#' tol = 1E-4,  maxiter = 20, type = "SARAR")
#' summary_SAR(result)
#' gamlss::term.plot(result$gamlss, what="mu")
#' 
#' @export
#' @importFrom methods is
#' @importFrom stats .getXlevels binomial dpois fitted gaussian glm.fit 
#' @importFrom stats make.link model.matrix model.offset model.response 
#' @importFrom stats optim pchisq pnorm predict printCoefmat rpois update
#' @import gamlss
#' @import splines
#' @importFrom gamlss.dist NO

SARARgamlss <- function(formula, sigma.formula = ~1,
                        W1 = diag(0, nrow(data)), W2 = diag(0, nrow(data)),
                        data, tol = 1E-4, maxiter = 20,
                        type = c("SAR", "SARAR", "SEM"),
                        weights = NULL) {
  mf <- stats::model.frame(formula, data=data)
  if (type == "SAR") {
    W2 = 0 * W2
  }
  if (type == "SEM") {
    W1 = 0 * W1
  }
  if (type == "SARAR" & sum(W2) == 0) {
    W2 <- W1
  }
  
  # Initialize variables
  rho0 <- 0
  lambda <- 0
  n <- nrow(data)
  
  # Fit the initial GAMLSS model for mean (mu) and variance (sigma)
  m0 <- gamlss::gamlss(formula = formula, sigma.formula = sigma.formula, 
                       data = data, family = NO())
  
  Y <- matrix(m0$y, ncol = 1)
  
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
    Ztemp <- model.matrix(m0, what = "sigma")
    
    colnames(Xtemp) <- colnames(model.matrix(m0, what = "mu"))
    colnames(Ztemp) <- colnames(model.matrix(m0, what = "sigma"))
    
    
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
    tolTemp <- sum(abs(p1 - p0))
    iter <- iter + 1
  }
  
  # Store updated model parameters
  
  names_to_update <- c("G.deviance", "residuals", "mu.fv", "mu.lp", "mu.wv", 
                       "mu.wt", "mu.qr", "mu.coefficients", "sigma.fv",
                       "sigma.wv", "sigma.wt", "sigma.coefficients", "sigma.qr",
                       "P.deviance", "aic", "sbc")
  
  for(valu in names_to_update){
    names(m1[[valu]]) <- names(m0[[valu]])
  }
  
  for (names in names_to_update) {
    m0[[names]] <- m1[[names]]
  }
  
  #model.frame(m0, what="mu") <- 
  # Return the final model object with spatial parameters and variance estimates
  if (type == "SARAR") {
    spamu <- c(p0[1],p0[2])
    var_cov_matrix <- solve(Hessian)
    spacov = var_cov_matrix
  } else if (type == "SAR") {
    spamu <- c(p0[1], NA)
    var_cov_matrix <- matrix(c( solve(Hessian[1,1]),NA, NA,NA), ncol=2, nrow=2)
    spacov = var_cov_matrix
  } else {
    spamu <- c(NA, p0[2])
    var_cov_matrix <- matrix(c(NA, NA,NA, solve(Hessian[2,2])), ncol=2, nrow=2)
    spacov =var_cov_matrix 
  }
  # model.frame(m0)[,1] <- as.vector(model.frame(m1)[,1])
  m0$call$data <- data 
  out1 <- list(gamlss=m0, model = mf, data=data,
               spatial=list(spatial=spamu, sdspatial=spacov, type=type),
               gamlssAY=m1)
  class(out1) <- "SARARgamlss"
  return(out1)
}

