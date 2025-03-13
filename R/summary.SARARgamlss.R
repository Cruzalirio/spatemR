#' @title Summary for SARARgamlss Model
#' @description This function provides a summary of a fitted SARARgamlss model, 
#' showing the generalized deviance, AIC, SBC, P-deviance, the coefficients for 
#' mean (\( \mu \)) and variance (\( \sigma \)), and spatial parameters (\( \hat{\rho} \) and \( \hat{\lambda} \)).
#' @param object An object of class `SARARgamlss`, typically the result from fitting a SARARgamlss model.
#' @param ... Additional arguments passed to the summary function (not used in this method).
#' @return A printed summary of the SARARgamlss model that includes:
#'   \item{Generalized Deviance}{The generalized deviance of the model.}
#'   \item{AIC}{The Akaike Information Criterion (AIC) of the model.}
#'   \item{SBC}{Schwarz Bayesian Criterion (SBC) for model comparison.}
#'   \item{P-deviance}{The P-deviance statistic of the model.}
#'   \item{Coefficients for Mean (mu)}{The coefficients for the mean (mu) model along with standard errors, t-values, and p-values.}
#'   \item{Coefficients for Variance (sigma)}{The coefficients for the variance (sigma) model along with standard errors, t-values, and p-values.}
#'   \item{Spatial Parameters}{Estimates of the spatial parameters, \(\hat{\rho}\) (spatial lag parameter) and \(\hat{\lambda}\) (spatial error parameter), 
#'   and their respective variances.}
#'   \item{Model Type}{The type of spatial model used (SAR, SARAR, or SEM).}
#'   \item{Iteration Status}{Information on the number of iterations used to fit the model.}
#' @details 
#'   The `summary.SARARgamlss` function provides a detailed output for the fitted `SARARgamlss` model.
#'   It includes a summary of the deviance, AIC, SBC, and P-deviance, along with the coefficients for the mean 
#'   and variance models. Additionally, it reports the spatial parameters \(\hat{\rho}\) and \(\hat{\lambda}\), 
#'   as well as their variances. The spatial model type (SAR, SARAR, or SEM) is also displayed. If the model 
#'   has been fitted successfully, the number of iterations is shown.
#' 
#'   The function is helpful for interpreting the results of a spatial GAMLSS model by displaying the key model
#'   parameters and spatial parameters, which are crucial in understanding the spatial effects in the model.
#' @seealso \link{SARARgamlss}
#' @examples
#' # Load data (example from the spdep package)
#' data(oldcol)
#' 
#' # Create spatial weight matrices W1 and W2
#' W1 <- spdep::nb2mat(COL.nb, style = "W")
#' W2 <- W1  # In this case, assume the same spatial weights for both
#' 
#' # Fit a SARARgamlss model
#' result <- SARARgamlss(formula = CRIME ~ INC + HOVAL, 
#'                       sigma.formula = ~ INC + pb(HOVAL), 
#'                       W1 = W1, W2 = W2, 
#'                       data = COL.OLD, 
#'                       tol = 1E-4, 
#'                       maxiter = 20, 
#'                       type = "SARAR")
#' 
#' # Print the summary
#' summary.SARARgamlss(result)
summary.SARARgamlss <- function(object, ...) {
  if (!inherits(object, "SARARgamlss")) {
    stop("Object must be of class 'SARARgamlss'")
  }
  
  # Print general summary information
  cat("SARARgamlss Model Summary\n")
  cat("---------------------------\n")
  cat("Coefficients for mu (mean):\n")
  print(object$mu.coefficients)
  cat("\n")
  cat("Coefficients for sigma (variance):\n")
  print(object$sigma.coefficients)
  cat("\n")
  
  # Rho and Lambda values
  cat("Estimated rho:", formatC(object$rho, digits = 4, format = "f"), "\n")
  cat("Estimated lambda:", formatC(object$lambda, digits = 4, format = "f"), "\n")
  
  # Deviance, AIC, etc.
  cat("GAMLSS Deviance:", formatC(object$G.deviance, digits = 4), "\n")
  cat("AIC:", formatC(object$aic, digits = 4), "\n")
  cat("SBC:", formatC(object$sbc, digits = 4), "\n")
  cat("Residual degrees of freedom:", object$df.residual, "\n")
}

