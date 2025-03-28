% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.GEESAR.R
\name{summary_SAR}
\alias{summary_SAR}
\title{Custom Summary Function for SARARgamlss and GEESAR Models}
\usage{
summary_SAR(object, ...)
}
\arguments{
\item{object}{An object of class 'SARARgamlss' or 'GEESAR'.}

\item{...}{Additional arguments (currently unused).}
}
\value{
A list containing the summary for the specified model class.
}
\description{
This function generates a summary for objects of class 'SARARgamlss' or 'GEESAR'.
It combines the summary outputs for both models, including GAMLSS model details, 
spatial parameters (rho and lambda), and Wald tests.
}
\examples{
\dontrun{
# Example for SARARgamlss model
result_sarar <- SARARgamlss(formula = CRIME ~ INC + HOVAL, 
                            sigma.formula = ~ INC + pb(HOVAL), 
                            W1 = W1, W2 = W2, data = data)
summary_SAR(result_sarar)

# Example for GEESAR model
result_geesar <- GEESAR(formula = CRIME ~ INC + HOVAL, data = data, W = W1)
summary_SAR(result_geesar)
}
}
