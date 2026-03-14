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

var_rho_inv <- function(A, W, X, beta, family, weights = NULL, phi = 1, offs = 0) {
  # Número de observaciones
  n <- nrow(X)
  
  # X ajustada por A
  Xi <- Matrix::solve(A, X)
  
  # Linear predictor
  eta <- as.vector(Xi %*% beta + offs)
  
  # Media
  mu <- family$linkinv(eta)
  
  # Primera derivada de la link
  g1 <- family$mu.eta(eta)
  
  # Derivadas de eta respecto rho
  deta <- as.vector(Matrix::solve(A, W %*% Xi %*% beta))
  d2eta <- as.vector(2 * Matrix::solve(A, W %*% Matrix::solve(A, W %*% Xi %*% beta)))
  
  # Derivadas de mu respecto rho
  dmu <- g1 * deta
  # Segunda derivada aproximada de g
  eps <- 1e-6
  g2 <- (family$mu.eta(eta + eps) - family$mu.eta(eta - eps)) / (2 * eps)
  d2mu <- g2 * deta^2 + g1 * d2eta
  
  fam <- family$family
  
  # Inicializar inv_var
  inv_var <- 0
  
  # Variancia inversa según familia, incluyendo phi siempre
  if(fam == "gaussian") {
    inv_var <- sum(dmu^2) / phi
  } else if(fam == "poisson"|| fam =="ptfamily") {
    inv_var <- sum(dmu^2 / mu) / phi
  } else if(fam == "binomial") {
    m <- weights
    if(is.null(m)) stop("weights (m) must be provided for binomial")
    inv_var <- sum(dmu^2 * (-m / mu^2 + (m - 1) / (1 - mu)^2)) / phi
  } else if(fam == "Gamma") {
    inv_var <- sum(dmu^2 * (1 / mu^2 - 1 / (mu * phi))) / phi
  } else if(fam == "Negative Binomial") {
    inv_var <- sum(dmu^2 / (mu + mu^2 / phi)) / phi
  } else {
    stop("Family not implemented")
  }
  
  # Devolver -inv_var según formula de Var^-1
  return(-inv_var)
}