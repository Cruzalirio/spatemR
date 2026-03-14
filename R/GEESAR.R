#' Generalized Estimating Equations with Spatial Autoregressive Components
#'
#' @description
#' `GEESAR` estimates generalized estimating equations (GEE) incorporating spatial autoregressive (SAR) components.  
#' It extends GEE models to account for spatial dependence in the response variable.
#'
#' @param formula A formula specifying the model structure (response ~ predictors).
#' @param family A description of the error distribution and link function. Default is `gaussian()`.
#' @param weights Optional vector of prior weights. Must be positive.
#' @param data A data frame containing the variables in the model.
#' @param W A spatial weights matrix defining the spatial dependence structure.
#' @param start Optional starting values for parameter estimation.
#' @param toler Convergence tolerance for iterative optimization. Default is `1e-05`.
#' @param maxit Maximum number of iterations for model fitting. Default is `50`.
#' @param trace Logical; if `TRUE`, prints iteration details. Default is `FALSE`.
#' @param eps Minimun value for variance diference of zero
#'
#' @details
#' The function estimates a spatially autoregressive GEE model by iteratively updating the spatial dependence 
#' parameter (`rho`) and regression coefficients (`beta`). The estimation follows a quasi-likelihood approach 
#' using iterative weighted least squares (IWLS).  
#'
#' The function supports common GLM families (`gaussian`, `binomial`, `poisson`, `Gamma`, `inverse.gaussian`) and 
#' their quasi-likelihood equivalents.
#'
#' @return A list of class `"GEESAR"` containing:
#' \item{coefficients}{Estimated regression coefficients.}
#' \item{rho}{Estimated spatial autoregressive parameter.}
#' \item{fitted.values}{Predicted values from the model.}
#' \item{linear.predictors}{Linear predictor values (`X * beta`).}
#' \item{prior.weights}{Weights used in estimation.}
#' \item{y}{Observed response values.}
#' \item{formula}{Model formula.}
#' \item{call}{Function call used to fit the model.}
#' \item{data}{Data used in the model.}
#' \item{converged}{Logical indicating whether the algorithm converged.}
#' \item{logLik}{Quasi-log-likelihood of the fitted model.}
#' \item{deviance}{Residual deviance.}
#' \item{df.residual}{Residual degrees of freedom.}
#' \item{phi}{Dispersion parameter estimate.}
#' \item{R}{Robust Variance Estimation.}
#' \item{CIC}{Corrected Information Criterion.}
#' \item{RJC}{Robust Jackknife Correction.}
#'
#' @seealso
#' \code{\link{glm}}, \code{\link[gee]{gee}}, \code{\link[spdep]{spdep}}
#'
#' @references Cruz, N. A., Toloza-Delgado, J. D., & Melo, O. O. (2024). 
#' Generalized spatial autoregressive model. arXiv preprint arXiv:2412.00945.
#' @source https://doi.org/10.48550/arXiv.2412.00945
#'
#' @examples
#' \donttest{
#' library(spdep)
#' library(sp)
#' data(meuse)
#' sp::coordinates(meuse) <- ~x+y
#' W <- spdep::nb2mat(knn2nb(knearneigh(meuse, k=5)), style="W")
#' fit <- GEESAR(cadmium ~ dist + elev, family=poisson(), data=meuse, W=W)
#' summary_SAR(fit)
#'}
#' @export
#' @importFrom sphet spreg

GEESAR <- function (formula, family = gaussian(), weights=NULL, data, W,
                    start = NULL, 
                    toler = 1e-04, maxit = 200, trace = FALSE, eps=1e-6) {
  mf <- stats::model.frame(formula, data=data)
  y <- base::as.matrix(model.response(mf))
  if (is(family, "function")) 
    family <- stats::family()
  if (family$family %in% c("quasi", "quasibinomial", "quasipoisson")) {
    if (family$family == "quasi") {
      family$family <- switch(family$varfun, constant = "gaussian", 
                              `mu(1-mu)` = "binomial", mu = "poisson", `mu^2` = "Gamma", 
                              `mu^3` = "inverse.gaussian")
    }
    else {
      family$family <- switch(family$family, quasibinomial = "binomial", 
                              quasipoisson = "poisson")
    }
    family <- do.call(family$family, list(link = family$link))
    family2 <- family
  }else {
    if (family$family %in% c("gaussian", "binomial", "poisson", 
                             "Gamma", "inverse.gaussian")) {
      varfun <- switch(family$family, gaussian = "constant", 
                       binomial = "mu(1-mu)", poisson = "mu", Gamma = "mu^2", 
                       inverse.gaussian = "mu^3")
      family2 <- do.call("quasi", list(variance = varfun, 
                                       link = family$link))
      if (family$family == "binomial") 
        family2 <- do.call("quasibinomial", list(link = family$link))
      if (family$family == "poisson") 
        family2 <- do.call("quasipoisson", list(link = family$link))
    }else family2 <- family
  }
  if (ncol(y) == 2 & family$family == "binomial") {
    weights <- as.matrix(y[, 1] + y[, 2])
    y <- as.matrix(y[, 1]/weights)
  }
  y <- as.matrix(as.numeric(y))
  offset <- as.vector(model.offset(mf))
  X <- model.matrix(formula, data)
  p <- ncol(X)
  n <- nrow(X)
  if (is.null(offset)){ 
    offs <- matrix(0, n, 1)}else{ offs <- as.matrix(offset)}
  if (is.null(weights)) {
    weights <- matrix(1, n, 1)} else {
      weights <- as.matrix(weights)}
  if (any(weights <= 0)) 
    stop("Only positive weights are allowed!!", call. = FALSE)
  datas <- cbind(X,y, offs, weights)
  
  
  qll <- NULL
  if (family$family == "gaussian") 
    qll <- function(weights, y, mu) -weights * (y - mu)^2/2
  if (family$family == "binomial") 
    qll <- qll <- function(weights, y, mu){
      eps <- eps
      mu <- pmin(pmax(mu, eps), 1-eps)
      weights * (y * log(mu) + (1 - y) * log(1 - mu))
    }
  if (family$family == "poisson") 
    qll <- function(weights, y, mu) weights * (y * log(mu) - mu)
  if (family$family == "Gamma") 
    qll <- function(weights, y, mu) -weights * (y/mu + log(mu))
  if (family$family == "inverse.gaussian") 
    qll <- function(weights, y, mu) weights * (mu - y/2)/mu^2
  if (family$family == "ptfamily") 
    qll <- function(weights, y, mu) weights * (y*log(mu)-mu-log(1-exp(-mu)))
  if (family$family == "Negative Binomial") {
    .Theta <- function() return(.Theta)
    environment(.Theta) <- environment(family$variance)
    .Theta <- .Theta()
    qll <- function(weights, y, mu) weights * (y * log(mu/(mu + .Theta)) + .Theta * 
                                                 log(.Theta/(.Theta + mu)))
  }
  if (family$family == "Tweedie") {
    .Theta <- function() return(p)
    environment(.Theta) <- environment(family$variance)
    .Theta <- .Theta()
    if (.Theta %in% c(0, 1, 2)) {
      if (.Theta == 0) 
        qll <- function(weights, y, mu) -weights * (y - mu)^2/2
      if (.Theta == 1) 
        qll <- function(weights, y, mu) weights * (y * log(mu) - mu)
      if (.Theta == 2) 
        qll <- function(weights, y, mu) -weights * (y/mu + log(mu))
    }
    else qll <- function(weights, y, mu) mu^(-.Theta) * (mu * y/(1 - .Theta) - mu^2/(2 - 
                                                                                       .Theta))
  }
  if (is.null(qll)) 
    qll <- function(weights, y, mu) family$dev.resids(y, mu, weights)
  
  qllp <- function(rho, beta, W, D,n){
    A <- Matrix::Diagonal(n) - rho*W
    Xi <-Matrix::solve(A + Matrix::Diagonal(n, eps),Matrix::Matrix(D[, 1:p]))
    yi <- D[, p + 1]
    ni <- nrow(Xi)
    etai <- Xi %*% beta + D[, p + 2]
    mui <- family$linkinv(etai[,1])
    if(family$family=="binomial"){
      mui <- pmin(pmax(mui,eps),1-eps)
    }
    QL <- sum(qll(weights=D[,p+3], y=D[,p+1], mu=mui))
    return(-QL)
  }
  
  if(family2$family=="ptfamily" &  !all(y>0)){
    stop("Only y>1 are allowed in ptfamily!!", call. = FALSE)
  }
  if (family$family=="poisson"){
  modTemp <- try(sphet::spreg(family$linkfun(y+1)~datas[, 1:p]-1,
                              listw = W, model = "lag"),
                 silent = TRUE)
  if(!inherits(modTemp, "try-error")){
    rho <- modTemp$coefficients[p+1]  
  }else{
    rho <- 0
  }}else if(family$family=="binomial"){
    glm_temp <- glm.fit(x = X,y = y,family = binomial(),
                        weights = weights,offset = offs)
    mu_glm <- glm_temp$fitted.values
    res <- y - mu_glm
    rho <- as.numeric( t(res) %*% (W %*% res) / (t(res) %*% res) )
    rho <- max(min(rho, 0.95), -0.95)
  }else{
    rho <- 0
  }
  A <- Matrix::Diagonal(n) - rho*W
  
  beta_new <- try(glm.fit(y = y, x = Matrix::solve(A + Matrix::Diagonal(n, eps),X), family = family2, 
                          weights = weights, offset = offs), silent = TRUE)
  if(!inherits(beta_new, "try-error")){
    beta_new <- matrix(beta_new$coefficients, nrow=p)
  }else if(!inherits(modTemp, "try-error")){
    beta_new <- matrix(modTemp$coefficients[1:p], nrow=p)
  }else{
    beta_new <- matrix(0, ncol=p)
  }
  
  rho_f <- optim(par=rho,qllp, beta=beta_new, W=W, D=datas,n=n, method="L-BFGS-B", lower=-0.99,upper=0.99)
  rho_new <- rho_f$par
  tolrho <- 1
  niterrho <- 1
 
  Result <- rho_new
  while(tolrho>toler & niterrho < maxit){
    rho <- rho_new
    A = Matrix::Diagonal(n) - rho*W
    
    tol <- 1
    niter <- 1
    
    while(tol > toler & niter < maxit){
      
      beta_old <- beta_new
      
      Xi <- Matrix::solve(A + Matrix::Diagonal(n,eps), X)
      
      eta <- Xi %*% beta_old + offs
      mu  <- family$linkinv(eta[,1])
      
      gprime <- family$mu.eta(eta[,1])
      varmu  <- family$variance(mu)
      
      wgee <- weights * gprime^2 / varmu
      
      H <- Matrix::crossprod(Xi, wgee * Xi)
      
      U <- Matrix::crossprod(Xi, weights* gprime/varmu * (y - mu))
      
      beta_new <- beta_old + Matrix::solve(H, U)
      
      tol <- max(abs((beta_new - beta_old)/(beta_old + 1e-8)))
      
      niter <- niter + 1
    }
    
    rho_f <- optim(par=rho,qllp, beta=beta_new, W=W, D=datas, n=n,
                   method="L-BFGS-B", lower=-0.99,upper=0.99, hessian = TRUE)
    
    rho_new <- rho_f$par
    varrho <- -rho_f$hessian[1,1]
    tolrho <- abs(rho_new-rho)
    
    niterrho <- niterrho + 1
    Result <- c(Result, rho_new)
  }
  #print(rho_f)
  rho <- rho_new
  A <- Matrix::Diagonal(n) - rho * W
  Xi <- Matrix::solve(A + Matrix::Diagonal(n, eps), X)
  eta <- Xi %*% beta_new + offs
  mu  <- family$linkinv(eta[,1])
  
  gprime <- family$mu.eta(eta[,1])
  varmu  <- family$variance(mu)
  
  wgee <- weights * gprime^2 / varmu
  
  # Información de Fisher aproximada (H) y score (U)
  H <- Matrix::crossprod(Xi, wgee * Xi)
  U <- Matrix::crossprod(Xi, weights * gprime / varmu * (y - mu))
  Ui <- Xi * as.vector(weights * gprime / varmu * (y - mu))
  
  # Varianza sandwich directamente
  I0 <- Matrix::solve(H)
  B <- crossprod(Ui)
  vcovs <- I0 %*% B %*% I0

  # Ajuste especial para binomial con weights = 1
  if(family$family == "binomial" & all(weights == 1)){
    vcovs <- I0
  }
  
  rownames(vcovs) <- colnames(X)
  colnames(vcovs) <- colnames(X)
  
  # Cálculo de CIC
  w <- sqrt(weights * gprime^2 / varmu)
  Xw <- sweep(Xi, 1, w, "*")
  phi <- sum((y - mu)^2 / (varmu / weights)) / (n - p)
  CIC <- sum(Matrix::diag(I0 %*% B))
  
  # RJC
  # Log-likelihood (quasi)
  logLik <- -qllp(rho = rho, W = W, D = datas, beta = beta_new, n = n)
  
  RJC <- 2*CIC-2*logLik/phi
  
  # Varianza de rho
  inv_var_rho <- var_rho_inv(A, W, X, beta_new, family,y, weights, phi)
  var_rho <- 1 / inv_var_rho[1]
  var_rho_san <- inv_var_rho[2]
  
  
  estfun <- Matrix::crossprod(Xw, (y - mu))
  rownames(estfun) <- colnames(X)
  colnames(estfun) <- ""
  out_ <- list(coefficients = beta_new, rho=rho, varrho=var_rho, 
               varrhoSan=var_rho_san,    fitted.values = mu, 
               linear.predictors = eta,  arrangedata = datas, vcovs=vcovs,
               prior.weights = weights, y = y, formula = formula, call = match.call(), 
               offset = offs, model = mf, data = data,  
               converged = ifelse(niterrho < maxit, TRUE, FALSE), estfun = estfun, 
               R = vcovs,naive = I0, family = family,  
               phi = phi,CIC = CIC, RJC = RJC, 
               logLik = logLik, deviance = sum(family$dev.resids(y, mu, weights)), 
               df.residual = length(y) - length(beta_new), 
               levels = .getXlevels(attr(mf, "terms"), mf),
               contrasts = attr(X, "contrasts"), 
               start = start, iter = niterrho, linear = TRUE)
  class(out_) <- "GEESAR"
  return(out_)
}



