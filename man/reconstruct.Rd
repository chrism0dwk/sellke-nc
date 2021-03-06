% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/reconstruct.R
\name{reconstruct}
\alias{reconstruct}
\title{Recalculate infection times using Sellke construction.}
\usage{
reconstruct(beta, epidemic, condition.on.I = F)
}
\arguments{
\item{beta}{The new infection rate}

\item{epidemic}{an epidemic simulated using \code{sellkeSim}}

\item{condition.on.I}{Should we condition on only known infections?}
}
\value{
A new epidemic, with rescaled infection times.
}
\description{
\code{reconstruct} returns a rescaled \code{epidemic}, conditioned on
removal times, Sellke thresholds, and a new \code{beta}.
}
\details{
This function takes a Sellke-defined epidemic, and rescales the infection times based on a
new transmission rate \code{beta}.  If the rescaling results in either a discontinuous
epidemic, or a violation of I < R, then the original epidemic is returned with a warning.
}

