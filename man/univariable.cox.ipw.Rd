% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Regression.Cox.IPW.R
\name{univariable.cox.ipw}
\alias{univariable.cox.ipw}
\title{Create IPW Univariate Cox PH Regression Table#'}
\usage{
univariable.cox.ipw(
  time,
  status,
  covariate,
  label,
  type,
  reflevel = NULL,
  data,
  weights = NULL
)
}
\arguments{
\item{time}{a variable for time to event}

\item{status}{a variable for event (or censoring) indicator (1 indicates event, 0 is censored)}

\item{covariate}{a vector of variables included in the univariable analysis}

\item{label}{a vector of labels shown in the table}

\item{type}{a vector of types either categorical or numerical (\code{"cat"} or \code{"num"})}

\item{reflevel}{a vector of reference levels for categorical variables. NA is used for numerical variable.}

\item{data}{data frame dataset}

\item{weights}{a variable for inverse probability weights}
}
\value{
matrix with estimates of coefficient, HR, and P-value using IPW
}
\description{
The function creates a result table for IPW univariate Cox PH regression.
}
