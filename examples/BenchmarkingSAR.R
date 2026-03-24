devtools::install_github("Cruzalirio/spatemR")
library(spdep)
library(sphet)
library(spatemR)
library(spatialreg)
library(gamlss)
library(Matrix)
library(tidyverse)
library(gamlss)

set.seed(123)

rho_grid <- c(-0.995,-0.8,-0.5,0,0.5,0.8,0.995)
n_grid   <- c(49,144,400,1600,10000)
R <- 10

b1 <- 1; b2 <- 0.5; b3 <- -0.3
a1 <- 0; a2 <- 0.4; a3 <- -0.2

results <- data.frame()

save_file <- "benchmark_partial.RData"

for(n in n_grid){
  #if(n < 9000) next
  cat("Running n =",n,"\n")
  
  data = data.frame(n1=1:n)
  data$lat = rep(1:sqrt(n), sqrt(n))
  data$long = sort(rep(1:sqrt(n), sqrt(n)))
  
  wt1 = as.matrix(dist(cbind(data$long, data$lat)))
  wt1 = ifelse(wt1 == 1, 1, 0)
  diag(wt1) = 0
  
  rs = rowSums(wt1)
  wt1 = apply(wt1, 2, function(x) x/rs)
  
  lw1 = mat2listw(wt1, style="W")
  
  for(rho1 in rho_grid){
    cat("  rho =", rho1, "\n")
    
    for(r in 1:R){
      
      inv2 = solve(Diagonal(n) - rho1*wt1)
      
      x1 = rnorm(n)
      x2 = rnorm(n, 2, 1)
      x3 = runif(n)
      
      data1 = data.frame(x1,x2,x3)
      
      mu = b1 + b2*x1 + b3*x2
      sigma = sqrt(exp(a1 + a2*x2 + a3*x3))
      
      data1$y1 = as.vector(inv2 %*% rnorm(n, mean=mu, sd=sigma))
      
      # ---------------------------
      # 1) Proposed — SARARgamlss
      # ---------------------------
      t0 <- Sys.time()
      fit <- SARARgamlss(
        y1 ~ x1 + x2,
        sigma.formula = ~ x2 + x3,
        W1 = wt1,
        data = data1,
        type = "SAR",
        tol = 1e-4,
        maxiter = 20
      )
      t1 <- Sys.time()
      results <- rbind(results, data.frame(
        n=n, rho=rho1, rep=r, time=as.numeric(t1-t0),
        model="Proposed"
      ))
      
      # ---------------------------
      # 2) Ro-SAR — sphet::spreg()
      # ---------------------------
      t0 <- Sys.time()
      fit2 <- spreg(y1 ~ x1 + x2, data = data1,
                    listw = lw1, model="sar")
      t1 <- Sys.time()
      results <- rbind(results, data.frame(
        n=n, rho=rho1, rep=r, time=as.numeric(t1-t0),
        model="Ro-SAR"
      ))
      
      # ---------------------------
      # 3) LAGSAR — spdep::lagsarlm()
      # ---------------------------
      t0 <- Sys.time()
      fit3 <- lagsarlm(y1 ~ x1 + x2, data=data1, listw=lw1)
      t1 <- Sys.time()
      results <- rbind(results, data.frame(
        n=n, rho=rho1, rep=r, time=as.numeric(t1-t0),
        model="LAGSAR"
      ))
      
      # ---------------------------
      # 4) GAMLSS — baseline
      # ---------------------------
      t0 <- Sys.time()
      fit4 <- gamlss(y1 ~ x1 + x2, sigma.formula = ~ x2 + x3,
                     data=data1, family=NO)
      t1 <- Sys.time()
      results <- rbind(results, data.frame(
        n=n, rho=rho1, rep=r, time=as.numeric(t1-t0),
        model="GAMLSS"
      ))
      
    }
    
    save(results, file=save_file)
    cat("Saved partial scenario.\n")
  }
}