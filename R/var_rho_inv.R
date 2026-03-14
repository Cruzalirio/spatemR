#' Compute the Inverse Variance of the Spatial Autoregressive Parameter (rho)
#'
#' This function calculates the inverse of the variance of the spatial autoregressive parameter 
#' \eqn{\rho} in a generalized spatial autoregressive (GSAR) or GEE-SAR model. 
#' The calculation is based on the quasi-likelihood derivatives with respect to \eqn{\rho} 
#' for different exponential family distributions.
#'
#' @param A Matrix. The spatial transformation matrix \eqn{\mathbf{A} = \mathbf{I} - \rho \mathbf{W}}, 
#' typically of class `Matrix`.
#' @param W Matrix. Row-standardized spatial weights matrix \eqn{\mathbf{W}} of dimension \eqn{n \times n}.
#' @param X Matrix. Design matrix of covariates, dimension \eqn{n \times p}.
#' @param beta Numeric vector. Current estimates of regression coefficients \eqn{\pmb{\beta}}, length \eqn{p}.
#' @param family GLM family object. The response distribution family (e.g., `gaussian()`, `poisson()`, `binomial()`, `Gamma()`, `Negative Binomial()`).
#' @param weights Numeric vector. Observation weights \eqn{m_i} (e.g., number of trials for binomial data), length \eqn{n}.
#' @param phi Numeric. Dispersion parameter, used for `gaussian`, `Gamma`, or `Negative Binomial` families. Default is 1.
#' @param offs Numeric vector. Optional offset vector, length \eqn{n}. Default is 0.
#'
#' @return Numeric. The inverse of the variance of \eqn{\hat{\rho}} (\eqn{\operatorname{Var}(\hat{\rho})^{-1}}).
#'
#' @details
#' The function computes first and second derivatives of the mean \eqn{\mu_i} with respect to 
#' \eqn{\rho}, and then applies the appropriate formula for the inverse variance based on the 
#' selected family. This generalizes the quasi-likelihood derivations for spatially correlated 
#' generalized linear models.
#'
#' For binomial families with large \eqn{m_i}, it is recommended to truncate \eqn{\mu_i} within 
#' \eqn{[1e-10, 1-1e-10]} to avoid numerical instability.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' n <- 10
#' W <- Matrix(0,n,n)
#' diag(W[-1,]) <- 1
#' X <- matrix(rnorm(n*2), n, 2)
#' beta <- c(0.5, -0.2)
#' rho <- 0.3
#' A <- Diagonal(n) - rho*W
#' family <- binomial()
#' weights <- rep(1,n)
#' var_rho_inv(A, W, X, beta, family, weights)
#' }
#'
#' @export

var_rho_inv <- function(A, W, X, beta, family, offs=NULL, weights=NULL, phi=1){
  n <- nrow(X)
  Ainv <- solve(A)
  tildeX <- Ainv %*% X
  eta <- as.vector(tildeX %*% beta + offs)
  mu <- family$linkinv(eta)
  g1 <- family$mu.eta(eta)
  
  # primera derivada
  deta <- as.vector(Ainv %*% W %*% tildeX)
  dmu <- g1 * deta
  
  fam <- family$family
  inv_var <- 0
  
  if(fam=="gaussian"){
    inv_var <- sum(dmu^2)/phi
  } else if(fam=="poisson"|| fam=="ptfamily"){
    inv_var <- sum(dmu^2 / mu)
  } else if(fam=="binomial"){
    m <- weights
    if(is.null(m)) stop("Weights must be provided for binomial")
    inv_var <- sum(dmu^2 * (-m/mu^2 + (m-1)/(1-mu)^2))
  } else if(fam=="Gamma"){
    inv_var <- sum(dmu^2 * (1/mu^2 - 1/(mu*phi)))
  } else if(fam=="Negative Binomial"){
    inv_var <- sum(dmu^2 / (mu + mu^2/phi))
  } else stop("Family not implemented")
  
  return(inv_var)
}