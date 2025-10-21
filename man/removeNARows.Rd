\name{removeNARows}
\alias{removeNARows}
\alias{removeNARows.default}
\alias{removeNARows.EList}
\title{Remove Entirely NA Rows from Matrix or EList}

\description{
Remove rows from a matrix that have fewer than a user-specified minimum number of non-missing observations.
}

\usage{
\method{removeNARows}{default}(y, nobs.min = 1, ...)
}

\arguments{
  \item{y}{a matrix or an EList object.}
  \item{nobs.min}{minimum number of non-missing observations for rows to be kept.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
Produces a new matrix keeping only those rows that have at least the specified number of non-missing values.
}

\value{
A matrix or EList the same as \code{y} but with entirely or mostly missing rows removed.
}

\examples{
y <- matrix(rnorm(25),5,5)
y[y < -0.5] <- NA
removeNARows(y)
}

\concept{Filter features}

