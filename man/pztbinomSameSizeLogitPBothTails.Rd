\name{pztbinomSameSizeLogitPBothTails}
\alias{pztbinomSameSizeLogitPBothTails}

\title{Zero-Truncated Binomial Distribution P-Values}

\description{
Left and right p-values for zero-truncated binomial deviates.
}

\usage{
pztbinomSameSizeLogitPBothTails(q, size, logit.prob)
}

\arguments{
  \item{q}{vector of quantiles.}
  \item{size}{number of trials (one or more). Should be a single value.}
  \item{logit.prob}{success probabilities on logit scale.}
}

\details{
This function returns left and right p-values for zero-truncated binomial deviates.
The \code{size} argument is assumed to be the same for all deviates.
This function is slow but accurate, even for very small or very large success probabilities.
}

\value{
A list with components:
  \item{left.p.value}{left tail p-values.}
  \item{right.p.value}{matrix of complete log2-expression values without NAs.}
  \item{total.prob}{data.frame with columns \code{Protein} and \code{DE.Status} giving protein ID and true DE status.}
}

\examples{
x <- 1:3
pztbinomSameSizeLogitPBothTails(x, size=3, logit.prob=0)
}

\keyword{distribution}
