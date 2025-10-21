\name{ztbinom}
\alias{ZeroTruncatedBinomial}
\alias{dztbinom}
\alias{pztbinom}
\title{Zero-Truncated Binomial Distribution}

\description{
Density and distribution function for the zero-truncated binomial distribution, using the same arguments as for the R stats binomial distribution functions.
}

\usage{
dztbinom(x, size, prob, log = FALSE, logit.p = FALSE)
pztbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
}

\arguments{
  \item{x,q}{vector of quantiles.}
  \item{size}{number of trials (zero or more).}
  \item{prob}{probability of success on each trial.}
  \item{log}{logical; if \code{TRUE}, the log-density is returned.}
  \item{logit.p}{logical; if \code{TRUE}, success probabilities \code{prob} are on the logit scale.}
  \item{lower.tail}{logical; if \code{TRUE}, probabilities are P(X<q) otherwise P(X>q).}
  \item{log.p}{logical; if \code{TRUE}, tail probabilities are on the log-scale.}
}

\details{
\code{dztbinom} and \code{pztbinom} perform similarly to the R stats functions \code{dbinom} and \code{pbinom} except for the zero-truncation.
}

\value{
Output values give density (\code{dztbinom}) or cumulative probability (\code{pztbinom}) values
for the zero-truncated binomial distribution with parameters \code{size} and \code{prob}.

Output is a vector of length equal to the maximum length of any of the arguments \code{x}, \code{q}, \code{size} or \code{prob}.
If the first argument is the longest, then all the attributes of the input argument are preserved on output, for example, a matrix \code{x} will give a matrix on output.
Elements of input vectors that are missing will cause the corresponding elements of the result to be missing, as will non-positive values for \code{size} or negative values for \code{prob}.
}

\examples{
# Compare to binomial
x <- 1:3
dztbinom(x, size=3, prob=0.5)
dbinom(x, size=3, prob=0.5)
pztbinom(x, size=3, prob=0.5)
pbinom(x, size=3, prob=0.5)
}

\concept{Zero truncated binomial distribution}
