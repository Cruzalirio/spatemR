#' Hurdle Model using GEESAR
#'
#' This function fits a hurdle model using GEESAR, consisting of:
#' (1) A logit model for zero vs. non-zero responses.
#' (2) A truncated Poisson model for positive counts.
#'
#' @param formula A formula specifying the model.
#' @param data The dataset.
#' @param W The spatial weight matrix.
#' @param weights Optional weights.
#' @param toler Convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @param trace Logical. If TRUE, prints progress.
#' @return A list containing the logit and Poisson-truncated models.
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 100
#' x <- rnorm(n)
#' y <- rpois(n, lambda = exp(0.5 * x))
#' y[rbinom(n, 1, 1/(1+exp(-0.5*x)))] <- 0  # Introduce zeros
#' W <- matrix(rbinom(n^2,1,0.2), n, n)  # Example spatial weight matrix
#' diag(W) <- 0
#' rtot <- rowSums(W)
#' W <- W/ifelse(rtot==0, 0.1, rtot)
#' model <- Hurdle_GEESAR(y ~ x, data = data.frame(y, x), W = W)
#' summary_SAR(model$logit_model)
#' summary_SAR(model$poisson_truncated_model)
#'}
#'
#'@export
Hurdle_GEESAR <- function(formula, data, W, weights = NULL, toler = 1e-05, maxit = 200, trace = FALSE) {
  
  # Modelo logit para y == 0 vs y > 0
  data$y_binary <- as.numeric(data[[all.vars(formula)[1]]] > 0)
  
  geesar_logit <- GEESAR(
    formula = update(formula, y_binary ~ .),
    family = binomial(link = "logit"),
    data = data,
    W = W,
    weights = weights,
    toler = toler,
    maxit = maxit,
    trace = trace
  )
  
  # Modelo Poisson truncado para y > 0
  data_positive <- subset(data, data[[all.vars(formula)[1]]] > 0)
  
  geesar_poisson_truncated <- GEESAR(
    formula = formula,
    family = ptfamily(),
    data = data_positive,
    W = W[data[[all.vars(formula)[1]]] > 0, data[[all.vars(formula)[1]]] > 0],
    weights = weights[data[[all.vars(formula)[1]]] > 0],
    toler = toler,
    maxit = maxit,
    trace = trace )
  
  return(list(logit_model = geesar_logit, poisson_truncated_model = geesar_poisson_truncated))
}

