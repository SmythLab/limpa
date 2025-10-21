\name{observedMomentsCN}
\alias{observedMomentsCN}
\title{Observed Distribution Moments from Complete Normal Model}

\description{
Mean and standard-deviation of the observed data distribution under the complete normal model.
}

\usage{
observedMomentsCN(mean.comp=6, sd.comp=1, dpc=c(-4,0.7))
}

\arguments{
  \item{mean.comp}{mean of complete normal distribution.}
  \item{sd.comp}{standard deviation of complete normal distribution.}
  \item{dpc}{numeric vector of length 2 giving the DPC intercept and slope.}
}

\details{
Under the complete normal model, calculate the mean and standard deviation of the observed data distribution.
}

\value{
A list with compoenents
  \item{mean.obs}{mean of observed data distribution.}
  \item{sd.obs}{standard deviation of observed data distribution.}
  \item{prob.obs}{unconditional probability that values are observed.}
}

\examples{
observedMomentsCN(mean.comp=6, sd.comp=2)
}

\concept{Missing value model}
