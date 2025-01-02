\name{simCompleteDataON}
\alias{simCompleteDataCN}
\alias{simCompleteDataON}
\title{Simulate Complete Data From Complete or Observed Normal Models}

\description{
Simulate a vector complete data together with the associated missing value events, under two different models.
}

\usage{
simCompleteDataCN(n, mean.comp=6, sd.comp=1, dpc=c(-4,0.7))
simCompleteDataON(n, mean.obs=6, sd.obs=1, dpc=c(-4,0.7))
}

\arguments{
  \item{n}{number of values to simulate.}
  \item{mean.comp}{mean of complete normal distribution.}
  \item{sd.comp}{standard deviation of complete normal distribution.}
  \item{mean.obs}{mean of observed normal distribution.}
  \item{sd.obs}{standard deviation of observed normal distribution.}
  \item{dpc}{numeric vector of length 2 giving the DPC intercept and slope.}
}

\details{
These functions simulate a vector of complete log2-expression data and identify which will be observed and which will be missing.
The complete values themselves are all non-missing, but some will be undetected in a hypothetical real dataset.
\code{simCompleteDataCN} simulates data according to the complete normal model (CN), while \code{simCompleteDataON} simulates data according to the observed normal model (ON).

These functions can be used to explore the differences between the complete and observed normal models.
Under the CN model, the complete values (including both observed and unobserved) are exactly normally distributed, while the subset that are observed are only approximately normal.
Under the ON model, the opposite is true. 
The observed values are exactly normal while the complete values are only approximately normal.
}

\value{
A list with compoenents
  \item{y.complete}{vector of complete values.}
  \item{is.missing}{vector of \code{TRUE} or \code{FALSE} values indicating whether each value will be missing.}
  \item{prob.missing}{conditional probability given \code{y.complete} that each value will be missing.}
}

\examples{
# Complete values are only approximately normal under the ON model.
out <- simCompleteDataON(100, mean.obs=6, sd.obs=1)
mean(out$prob.missing)
qqnorm(out$y.complete)
qqline(out$y.complete)
}
