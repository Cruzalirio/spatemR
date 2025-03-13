GEESAR <- function (formula, family = gaussian(), weights, data, W,
                    corr, start = NULL, 
          toler = 1e-05, maxit = 50, trace = FALSE, ...) {
  mf <- model.frame(formula, data=data)
  y <- as.matrix(model.response(mf))
  if (is(family, "function")) 
    family <- family()
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
    qll <- function(weights, y, mu) weights * (y * log(mu) + (1 - y) * log(1 - mu))
  if (family$family == "poisson") 
    qll <- function(weights, y, mu) weights * (y * log(mu) - mu)
  if (family$family == "Gamma") 
    qll <- function(weights, y, mu) -weights * (y/mu + log(mu))
  if (family$family == "inverse.gaussian") 
    qll <- function(weights, y, mu) weights * (mu - y/2)/mu^2
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
    A = diag(n) - rho*W
    Xi <-solve(A)%*%matrix(D[, 1:p], ncol = p)
    yi <- D[, p + 1]
    ni <- nrow(Xi)
    etai <- Xi %*% beta + D[, p + 2]
    mui <- family$linkinv(etai[,1])
    QL <- sum(qll(weights=D[,p+3], y=D[,p+1], mu=mui))
    return(-QL)
  }
  score <- function(D, W, rho, beta, ni) {
    A = diag(ni) - rho*W
    Xi <-solve(A)%*%matrix(D[, 1:p], ncol = p)
    yi <- D[, p + 1]
    ni <- nrow(Xi)
    etai <- Xi %*% beta + D[, p + 2]
    mui <- family$linkinv(etai[,1])
    wi <- sqrt(family$variance(mui)/D[, p + 3])
    Xiw <- Xi * matrix(family$mu.eta(etai[,1]), nrow(Xi), p)
    Vi <- t(solve(crossprod(A)) * matrix(wi, ni, ni)) * 
      matrix(wi, ni, ni)
    Vi2 <- chol2inv(chol(Vi))
    Xiw2 <- crossprod(Vi2, Xiw)
    cbind(crossprod(Xiw2, (yi - mui)), crossprod(Xiw2,Xiw), 
               crossprod(Xiw2, (tcrossprod(yi - mui)) %*% Xiw2))
  }
  
  rho <- 0
  A = diag(n) - rho*W
  
  beta_new <- try(glm.fit(y = y, x = solve(A)%*%X, family = family2, 
                            weights = weights, offset = offs), silent = TRUE)
  beta_new <- matrix(beta_new$coefficients, nrow=p)
  rho_f <- optim(par=rho,qllp, beta=beta_new, W=W, D=datas,n=n, method="L-BFGS-B", lower=-0.99,upper=0.99)
  rho_new <- rho_f$par
  tolrho <- 1
  niterrho <- 1
  
  # rrr = c()
  # for(rho in seq(-1,0.99,0.01)) rrr = c(rrr, qllp(rho=rho, W=W, D=datas, beta=beta_new, n=n))
  # plot(seq(-1,0.99, 0.01),rrr, type="l")
  # 
  Result <- rho_new
  while(tolrho>toler & niterrho < maxit){
  rho <- rho_new
  A = diag(n) - rho*W
  beta_new <- try(glm.fit(y = y, x = solve(A)%*%X, family = family2, 
                            weights = weights, offset = offs), silent = TRUE)
  beta_new <- matrix(beta_new$coefficients, nrow=p)
  
  niter <- 1
  tol <- 1
  while (tol > toler & niter < maxit) {
    beta_old <- beta_new
    resume2 <- score(D=datas,W=W, rho=rho, beta = beta_old, ni=n)
    kchol <- chol2inv(chol(resume2[1:p, 2:(p + 1)]))
    beta_new <- beta_old + crossprod(kchol, resume2[, 1])
    tol <- max(abs((beta_new - beta_old)/beta_old))
    niter <- niter + 1
  }
  
  rho_f <- optim(par=rho,qllp, beta=beta_new, W=W, D=datas, n=n,
                 method="L-BFGS-B", lower=-0.99,upper=0.99)
  rho_new <- rho_f$par
  tolrho <- abs(rho_new-rho)
  niterrho <- niterrho + 1
  Result = c(Result, rho_new)
  }
  
  #plot(Result, type="l")
  
  rho <- rho_new
  A = diag(n) - rho*W
  rownames(beta_new) <- colnames(X)
  colnames(beta_new) <- ""
  if (niter == maxit) 
    warning("Iteration limit exceeded!!\n", call. = FALSE)
  eta <- solve(A)%*% X %*% beta_new + offs
  mu <- family$linkinv(eta[,1])
  phiis <- diag(solve(crossprod(A)))
  vari <- sqrt(family$variance(mu)/datas[, p + 3]*phiis)
  phi <- sum(((y-mu)/sqrt(vari))^2)/(n-p)
  sera <- try(chol(R), silent = TRUE)
  I0 <- solve(resume2[1:p, 2:(p + 1)])
  I1 <- resume2[1:p, (p + 2):(2 * p + 1)]
  RJC <- crossprod(I0, I1)
  vcovs <- RJC %*% I0
  rownames(vcovs) <- colnames(X)
  colnames(vcovs) <- colnames(X)
  
  w <- sqrt(weights * family2$mu.eta(eta[,1])^2/family2$variance(mu))
  Xw <- matrix(w, nrow(X), ncol(X)) * X
  CIC <- sum(diag((crossprod(Xw)/phi) %*% vcovs))
  RJC <- sqrt((1 - sum(diag(RJC))/(p * phi))^2 + (1 - sum(diag(RJC %*% 
                                                                 RJC))/(p * phi^2))^2)
  logLik <- -qllp(rho=rho, W=W, D=datas, beta=beta_new, n=n)/phi
  estfun <- as.matrix(resume2[, 1])
  rownames(estfun) <- colnames(X)
  colnames(estfun) <- ""
  out_ <- list(coefficients = beta_new, rho=rho,  fitted.values = mu, 
               linear.predictors = eta,  arrangedata = datas, 
               prior.weights = weights, y = y, formula = formula, call = match.call(), 
               offset = offs, model = mf, data = data, 
               score = score, converged = ifelse(niter < maxit, TRUE, FALSE), estfun = estfun, 
               naive = I0, family = family,  phi = phi, phiis = phiis, CIC = CIC, RJC = RJC, 
               logLik = logLik, deviance = sum(family$dev.resids(y, mu, weights)), df.residual = length(y) - length(beta_new), 
               levels = .getXlevels(attr(mf, "terms"), mf),
               contrasts = attr(X, "contrasts"), start = start, iter = niter, linear = TRUE)
  class(out_) <- "glmgee"
  return(out_)
}



