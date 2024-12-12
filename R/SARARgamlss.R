#' @title Spatial autoregresive models with joint modelling of mean and variance
#' 
#' @description
#' Function to calculate a SARAR, SEM or SAR model with explicit model using 
#' GAMLSS 
#' 


SARARgamlss <- function(formula = formula(data), sigma.formula = ~1,
                        W1= diag(0,nrow(data)), W2= diag(0,nrow(data)),
                        data, tol = 1E-4, maxiter = 20,
                        type = c("SAR", "SARAR", "SEM"),
                        weights = NULL,...){
  if(type=="SAR" & sum(W1)==0){
    print("The SAR Model contains W1=0, is a usual non-spatial GAMLSS")
  }
  if(type=="SAR"){
    W2 = 0*W2
  }
  if(type=="SEM"){
    W1 = 0*W1
  }
  if(type=="SEM" & sum(W2)==0){
    print("The SEM Model contains W2=0, is a usual non-spatial GAMLSS")
  }
  if(type=="SARAR" & sum(W2)==0& sum(W1)==0){
    print("The SAR Model contains W1=0 and W2=0, is a usual non-spatial GAMLSS")
  }
  if(type=="SARAR" & sum(W2)==0){
    print("The SARAR Model assume that W1=W2")
  }
  rho0 <- 0
  lambda <- 0
  n <- nrow(data)
  m0 <- gamlss::gamlss(formula=formula, sigma.formula = sigma.formula, data=data, family = NO())
  Y <- matrix(m0$y, ncol=1)
  namesX <- colnames(model.matrix(m0, what="mu"))
  namesZ <- colnames(model.matrix(m0, what="sigma"))
  var0 <- predict(m0,what="sigma",type="response")^2
  Xbeta <- predict(m0, what = "mu", type="response")
  loglik <- function(rholam, W1, W2, Xbeta,Y, var0){
    AA <- diag(n)-rholam[1]*W1
    BB <- diag(n)-rholam[2]*W2
    VV <- BB%*%(AA%*%Y-Xbeta)/sqrt(var0)
    loglik<- -0.5*sum(log(var0)) + log(det(AA)) +log(det(BB))-
      0.5*sum(VV^2)
    return(-loglik)
  }
  p0 <- optim(par=c(0,0),fn=loglik, method="L-BFGS-B", W1=W1,W2=W2, Xbeta=Xbeta, 
              Y=Y, var0=var0,
              lower =c(-0.999, -0.999),upper=c(0.999, 0.99)  )
  p0 <- p0$par
  tolTemp <- 1
  iter <- 1
  while(tolTemp > tol& iter<maxiter){
    p1 <- p0
    AA <- diag(n)-p1[1]*W1
    BB <- diag(n)-p1[2]*W2
    Ytemp <- as.matrix(BB%*%AA%*%Y)
    Xtemp <- as.matrix(BB%*%model.matrix(m0, what="mu"))
    colnames(Xtemp) <- namesX
    Ztemp <- model.matrix(m0, what="sigma")
    colnames(Ztemp) <- namesZ
    m1 <- gamlss::gamlss(Ytemp~Xtemp-1,
                         ~Ztemp-1)
    var1 <- predict(m1,what="sigma",type="response")^2
    Xbeta <- predict(m1, what = "mu", type="response")
    p0 <- optim(par=c(0,0),fn=loglik, method="L-BFGS-B", W1=W1,W2=W2, Xbeta=Xbeta, 
                Y=Y, var0=var1,
                lower =c(-0.999, -0.999),upper=c(0.999, 0.99)  )
    p0 <- p0$par
    tolTemp <- sum(abs(p1-p0))
    iter <- iter+1
  }
  namesTemp <- c("G.deviance", "residuals","mu.fv","mu.lp",
                 "mu.wv","mu.wt","mu.qr" ,"mu.coefficients",
                 "sigma.fv","sigma.wv","sigma.wt",
                 "sigma.coefficients","P.deviance","aic","sbc") 
  for(names in namesTemp){
    m0[[names]] <- m1[[names]]
  }
  if(type=="SARAR"){
  m0$rho <- p0[1]
  m0$lambda <- p0[2]}else if(type=="SAR"){
    m0$rho <- p0[1]
  }else{
    m0$lambda <- p0[2]
  }
  return(m0)
}
