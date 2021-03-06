% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.r
\name{make.adeba}
\alias{make.adeba}
\title{Initialize a density estimate}
\usage{
make.adeba(data, range = c(Inf, Inf), alpha = NULL, beta = 0.5, pilot,
  transform = TRUE, na.rm = FALSE, parallel = FALSE,
  log_prior = uniform_log_prior)
}
\arguments{
\item{data}{Dataset.}

\item{range}{Domain of each variable (column) of the dataset.
The idea is to allow for bounded kernels in the future, but at the moment
this argument has no effect.}

\item{alpha}{Alpha values controlling the global bandwidth scaling.
This should be set to \code{NULL}, implying automatic estimation from the
data, unless you know what you are doing.}

\item{beta}{See \code{\link{adeba}}.}

\item{pilot}{Pilot function values for the first iteration. Leave unspecified
to use a constant pilot.}

\item{transform}{If estimating multivariate densities, this argument controls
whether the dataset should be pre-rotated and scaled before estimation.
If the covariance matrix of \code{data} is large this is a very good idea,
but it takes a litte extra time to do.}

\item{na.rm}{Whether to remove observations with missing values (\code{TRUE})
or throw an error (\code{FALSE}).}

\item{parallel}{See \code{\link{adeba}}.}

\item{log_prior}{A function that calculates a log prior from a data frame
  with parameters. See \code{\link{log_prior}} for options and details.

  \emph{NOTE:}\cr
  Unless you know ADEBAs internals well you should probably not
  touch this argument. It was only exposed to facilitate analyses on the
  prior's importance, or rather lack of it, which were added as a supplement
  to the original publication.}
}
\value{
An unfitted density estimate that it to be passed to
  \code{\link{iterate}}.
  The estimate consist of the following components:
  \describe{
    \item{\code{iterations}}{Number of iterations calculated.}
    \item{\code{parameters}}{All the parameters and posterior values.}
    \item{\code{bandwidths}}{Bandwidths corresponding to the parameters.
      These are pre-calculated to make the results easier for the user to
      digest and manipulate, and since are often needed multiple times.}
    \item{\code{distance}}{Distance matrix of the original data set.}
    \item{\code{constant}}{Logical vector marking columns in the data set
      without any variation. These are excluded from the analysis.}
    \item{\code{transform}}{To be able to use spherical kernels on
      multivariate data sets with variables of very different variances or
      high linear dependence, a PCA-based transformation is applied to the
      data prior to computation. This element contains that function.}
    \item{\code{retransform}}{Function for converting transformed data back
      to the original domain. Needed by \code{\link{radeba}}.}
  }
}
\description{
Initialize a density estimate
}
\author{
Christofer \enc{Bäcklin}{Backlin}
}
