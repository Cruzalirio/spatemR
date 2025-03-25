#' Truncated Poisson Family for GLM
#'
#' This function defines a truncated Poisson family for use in Generalized Linear Models (GLMs), 
#' where zero values are not allowed. It modifies the Poisson likelihood by excluding zero-count observations.
#'
#' @param link Character string or a link-glm object specifying the link function. Accepted values are "log", "identity", and "sqrt".
#' @return An object of class "family" that can be used in \code{glm()}.
#'
#' @examples
#' set.seed(123)
#' y <- rpois(100, lambda = 3)
#' y <- y[y > 0]  # Truncate zeros
#' x <- rnorm(length(y))
#' model <- glm(y ~ x, family = ptfamily())
#' summary(model)
#'

ptfamily <- function (link = "log") 
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  family <- "ptfamily"
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s", 
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu * (1 - exp(-mu)) / (1 - exp(-mu) - mu * exp(-mu))
  validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
  dev.resids <- function(y, mu, wt) {
    -2 * wt * (y * log(mu) - mu - lfactorial(y) - log(1 - exp(-mu)))
  }
  aic <- function(y, n, mu, wt, dev) -2 * sum(dpois(y, lambda = mu, log = TRUE) - log(1 - exp(-mu)))
  initialize <- expression({
    if (any(y < 1)) stop("Values less than 1 not allowed for the 'Poisson truncated' family")
    n <- rep.int(1, nobs)
    mustart <- y + 0.1
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) 
      warning("ignoring prior weights")
    ftd <- fitted(object)
    rpois(nsim * length(ftd), ftd)
  }
  structure(list(family = family, link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, 
                 dispersion = 1), class = "family")
}

