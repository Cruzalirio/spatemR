#' Custom Summary Function for SARARgamlss and GEESAR Models
#'
#' This function generates a summary for objects of class 'SARARgamlss' or 'GEESAR'.
#' It combines the summary outputs for both models, including GAMLSS model details, 
#' spatial parameters (rho and lambda), and Wald tests.
#'
#' @param object An object of class 'SARARgamlss' or 'GEESAR'.
#' @param ... Additional arguments (currently unused).
#' @return A list containing the summary for the specified model class.
#' @examples
#' \dontrun{
#' # Example for SARARgamlss model
#' result_sarar <- SARARgamlss(formula = CRIME ~ INC + HOVAL, 
#'                             sigma.formula = ~ INC + pb(HOVAL), 
#'                             W1 = W1, W2 = W2, data = data)
#' summary_custom(result_sarar)
#' 
#' # Example for GEESAR model
#' result_geesar <- GEESAR(formula = CRIME ~ INC + HOVAL, data = data, W = W1)
#' summary_custom(result_geesar)
#' }
#' @export

summary_SAR <- function(object, ...) {
  
  # Verificar si el objeto es de clase 'SARARgamlss'
  if (!is.null(object$gamlss)) {
    
    # Resumen del modelo GAMLSS dentro del objeto SARARgamlss
    sum_gamlss <- summary(object$gamlss, type="qr")
    
    # Extraer parámetros espaciales de SARARgamlss
    rho_hat <- object$spatial$spatial[1]
    lambda_hat <- object$spatial$spatial[2]
    var_rho <- object$spatial$sdspatial[1, 1]
    var_lambda <- object$spatial$sdspatial[2, 2]
    
    # Calcular estadísticas de Wald para rho y lambda
    z_rho <- rho_hat / sqrt(var_rho)
    z_lambda <- lambda_hat / sqrt(var_lambda)
    
    spatial_summary <- data.frame(
      Estimate = c(rho_hat, lambda_hat),
      `Std. Error` = c(sqrt(var_rho), sqrt(var_lambda)),
      `Wald Statistic` = c(z_rho, z_lambda),
      `Pr(>|W|)` = 2 * pnorm(abs(c(z_rho, z_lambda)), lower.tail = FALSE)
    )
    rownames(spatial_summary) <- c("rho", "lambda")
    colnames(spatial_summary) <- c("Estimate", "Std.Error", "z value", "Pr(>|z|)")
    
    # Realizar prueba conjunta de Wald
    waldT <- matrix(c(rho_hat, lambda_hat), ncol = 2) %*% solve(object$spatial$sdspatial) %*% 
      matrix(c(rho_hat, lambda_hat), ncol = 1)
    pvalue <- pchisq(waldT, df = 2, lower.tail = FALSE)
    wald_test <- cbind(Wald = waldT, `Pr(>W)` = pvalue)
    colnames(wald_test) <- c("Wald value", "Pr(>w)")
    rownames(wald_test) <- c("c(rho, lambda)=c(0,0)")
    
    # Crear lista de resultados
    summary_obj <- list(
      gamlss_summary = sum_gamlss,
      spatial_summary = spatial_summary,
      wald_test = wald_test
    )
    
    class(summary_obj) <- "summary.SARARgamlss"
    return(summary_obj)
  }
  
  # Verificar si el objeto es de clase 'GEESAR'
  else if (inherits(object, "GEESAR")) {
    
    beta_hat <- object$coefficients
    vcov_matrix <- object$naive
    std_error <- sqrt(diag(vcov_matrix))
    Z_stat <- beta_hat / std_error
    p_values <- 2 * (1 - pnorm(abs(Z_stat)))
    
    rho_hat <- object$rho
    varrho_hat <- object$varrho  
    
    coef_table <- cbind(
      Estimate = beta_hat,
      Std.Error = std_error,
      Z.value = Z_stat,
      P.value = p_values
    )
    colnames(coef_table) <- c("Estimate", "Std.Error", "z value", "P(>|z|)")
    
    summary_obj <- list(
      coefficients = coef_table,
      rho = rho_hat,
      varrho = varrho_hat,
      dispersion = object$phi,
      df.residual = object$df.residual,
      call = object$call,
      family = object$family$family,
      logLik = object$logLik,
      CIC = object$CIC,
      RJC = object$RJC
    )
    
    class(summary_obj) <- "summary.GEESAR"
    return(summary_obj)
  }
  
  # Si el objeto no pertenece a ninguna de las clases, lanzar error
  else {
    stop("The object must be of class 'SARARgamlss' or 'GEESAR'.")
  }
}

#' Print Method for Summary of SARARgamlss Models
#'
#' This method prints a formatted summary of a 'summary.SARARgamlss' object,
#' including details of the GAMLSS model, spatial parameters (rho and lambda), 
#' and Wald tests.
#'
#' @param x An object of class 'summary.SARARgamlss'.
#' @param ... Additional arguments (currently unused).
#' @export
print.summary.SARARgamlss <- function(x, ...) {
  cat("\nSpatial Parameters Summary:\n")
  print(x$spatial_summary)
  cat("---------------------------------------------------------\n")
  
  # Prueba conjunta de Wald
  if (!is.null(x$wald_test)) {
    cat("\nJoint Wald Test for Spatial Parameters (rho, lambda):\n")
    print(x$wald_test)
  }
  
  cat("=========================================================\n")
}

#' Print Method for Summary of GEESAR Models
#'
#' This method prints a formatted summary of a 'summary.GEESAR' object,
#' including details of the model coefficients, rho, dispersion, and other statistics.
#'
#' @param x An object of class 'summary.GEESAR'.
#' @param ... Additional arguments (currently unused).
#' @export
print.summary.GEESAR <- function(x, ...) {
  cat("\nSummary of GEESAR Model\n")
  cat("=========================================================\n")
  cat("Family:", x$family, "\n")
  cat("Call:", deparse(x$call), "\n")
  cat("---------------------------------------------------------\n")
  cat("Coefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE)
  cat("---------------------------------------------------------\n")
  cat("Estimated rho:", formatC(x$rho, digits = 4, format = "f"), "\n")
  cat("Estimated variance of rho (varrho):", formatC(x$varrho, digits = 4, format = "f"), "\n")
  cat("Estimated dispersion:", formatC(x$dispersion, digits = 4, format = "f"), "\n")
  cat("Log-likelihood:", formatC(x$logLik, digits = 4, format = "f"), "\n")
  cat("CIC:", formatC(x$CIC, digits = 4, format = "f"), "\n")
  cat("RJC:", formatC(x$RJC, digits = 4, format = "f"), "\n")
  cat("Residual degrees of freedom:", x$df.residual, "\n")
  cat("=========================================================\n")
}

