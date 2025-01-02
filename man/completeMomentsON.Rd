\name{completeMomentsON}
\alias{completeMomentsON}
\title{Complete Distribution Moments from Observed Normal Model}

\description{
Mean and standard-deviation of the complete data distribution under the observed normal model.
}

\usage{
completeMomentsON(mean.obs=6, sd.obs=1, dpc=c(-4,0.7))
}

\arguments{
  \item{mean.obs}{mean of observed normal distribution.}
  \item{sd.obs}{standard deviation of observed normal distribution.}
  \item{dpc}{numeric vector of length 2 giving the DPC intercept and slope.}
}

\details{
Under the observed normal model, calculate the mean and standard deviation of the complete data distribution that would have occurred if the missing value mechanism hadn't operated.
}

\value{
A list with compoenents
  \item{mean.comp}{mean of complete data distribution.}
  \item{sd.comp}{standard deviation of complete data distribution.}
  \item{prob.obs}{unconditional probability that values are observed.}
}

\examples{
completeMomentsON(mean.obs=6, sd.obs=2)
}
