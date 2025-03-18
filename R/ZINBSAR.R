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



estimate_parameters <- function(y, z, x, W, max_iter = 100,
                                tol = 1e-6, epsilon=1e-6, tauReal=NULL) {
  n <- length(y)
  rho <- c(0, 0, 1)
  
  mod0 <- ZIM::zim.fit(y,x,z)
  tau <- 1
  beta <- mod0$para[1:(ncol(x))]
  gamma <- mod0$para[(ncol(x)+1):(ncol(x)+ncol(z))]
  # Initialize matrices
  
  for (iter in 1:max_iter) {
    rho1 <- rho[1]
    rho2 <- rho[2]
    A1 <- diag(n) - rho1 * W
    A2 <- diag(n) - rho2 * W
    # Update p and lambda
    p <- as.vector(1 / (1 + exp(-solve(A1) %*% z %*% gamma)))
    lambda <- as.vector(exp(solve(A2) %*% x %*% beta))
    
    # Compute u
    u <- (1 + ((1 - p) * (1 / (1 + tau * lambda))^(1/tau)) / p)^(-1) * (y == 0)
    
    # Compute gradients
    dp_dgamma <- (solve(A1) %*% z)*matrix(p * (1 - p), nrow=n, ncol=ncol(z), byrow=FALSE)
    dlambda_dbeta <- (solve(A2) %*% x) * matrix(lambda, nrow=n, ncol=ncol(x), byrow=FALSE)
    
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
    # ltheta <- function(theta){
    #   loglik = (1 - u) * dnbinom(y, size = exp(-theta), mu = lambda, log = TRUE)
    #   return(-sum(loglik))
    # }
    # thetaTemp <- optim(0, fn=ltheta, method= "L-BFGS-B", hessian=TRUE)
    # #resTemp = y-lambda
    # if(is.null(tauReal)){
    #   tau <- sum(lambda^2 * (1 - u)^2 * ((y- lambda)^2-lambda))/ sum((1 - u)^2 * lambda^4)
    #   tau = 1/( MASS::theta.ml(y[y>0], lambda[y>0], n=length(y[y>0])))
    # }else{
    #   tau=tauReal
    # }
    # 
    #tau <-  abs(tau)
    
    #tau <- exp(thetaTemp$par)
    # Update rho1 and rho2 simultaneously using Newton-Raphson
    likelihood_grad <- function(rho) {
      rho1 <- rho[1]
      rho2 <- rho[2]
      tau <- rho[3]
      A1 <- diag(n) - rho1 * W
      A2 <- diag(n) - rho2 * W
      p <- as.vector(1 / (1 + exp(-solve(A1) %*% z %*% gamma)))
      #p <- ifelse(p==1, 0.999, p)
      lambda <- as.vector(exp(solve(A2) %*% x %*% beta))
      u <- (1 + ((1 - p) * (1 / (1 + tau * lambda))^(1/tau)) / p)^(-1) * (y == 0)
      l0 <- u*log(p+epsilon)+(1-u)*log(1-p+epsilon)
      l1 <-  (1 - u) * dnbinom(y, size = 1/tau, mu = lambda, log = TRUE)
      log_likelihood <- sum(l0+l1)
      return(-log_likelihood)
    }
    
    parTemp <- optim(par=rho, fn=likelihood_grad, method= "L-BFGS-B", lower=c(-0.99, -0.99, 0.01),
          upper=c(0.99, 0.99, 200), hessian=TRUE)
    
    rhoNew <- parTemp$par[1:2]
    tau <- parTemp$par[3]

    # Convergence check
    if (max(abs(gamma_update), abs(beta_update), abs(parTemp$par- rho)) < tol) {
      break
    }
    rho =parTemp$par
    print(rho)
  }
  
  return(list(gamma = gamma, beta = beta, rho = rhoNew, 
              vcovrho=solve(parTemp$hessian), tau = tau))
}





library(spdep)
library(Matrix)



n = 144
beta = matrix(c(0.1,-0.05,0.05), ncol=1)
gamma = matrix(c(-0.1, -0.05, 0.05), ncol=1)
rho1 = 0
rho2 = 0
data = data.frame(n1=1:n)
# coords
data$lat = rep(1:sqrt(n), sqrt(n))
data$long = sort(rep(1:sqrt(n), sqrt(n)))
# create W matrix
wt = as.matrix(dist(cbind(data$long, data$lat), method = "euclidean", upper=TRUE))
wt1 = ifelse(wt==1, 1, 0)
diag(wt1) = 0

# row standardize
rs = rowSums(wt1)
wt1 = apply(wt1, 2, function(x) x/rs)
lw1 = mat2listw(wt1, style="W")

#wt1 = Matrix(wt1, sparse = TRUE)
inv1 = solve(diag(n)-rho1*wt1)
inv2 = solve(diag(n)-rho2*wt1)
x2 = runif(n)
x3 = runif(n)

X<-cbind(rep(1,n),x2,x3)
Z <- cbind(rep(1,n),x2,x3)

data$eta1 = as.vector(inv1%*%Z%*%gamma)
data$p1 = exp(data$eta1)/(1+exp(data$eta1))

data$nu2 = as.vector(inv2%*%X%*%beta)
data$lambda2 = exp(data$nu2)

tau <- 5

data$y1 = ifelse(runif(n)<data$p1, 0, rnbinom(n=length(data$lambda2), 
                                              mu= data$lambda2, size=1/tau))

y=data$y1
x=X
z=X
W=wt1
epsilon = 0.00001
max_iter = 20
tol = 1e-4


dataM = cbind(data, X)
result <- estimate_parameters(y=data$y1,x=X,z=X,W=wt1, max_iter = 200, tol=1E-4)
result$gamma
result$beta
