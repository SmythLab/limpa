\name{estimateDPCIntercept}
\alias{estimateDPCIntercept}
\title{Estimate DPC Intercept}

\description{
Estimate the DPC intercept given a value for the slope.
}

\usage{
estimateDPCIntercept(y, dpc.slope = 0.8, trace = FALSE)
}

\arguments{
  \item{y}{numeric matrix of log2-intensities, or any data object than can be coerced to a matrix.
    Includes NAs.
    Rows correspond to peptide precursors and columns to samples.}
  \item{dpc.slope}{DPC slope.}
  \item{trace}{if \code{TRUE}, then progress information will be printed from each glm iteration.}
}

\details{
Estimates the intercept coefficient of the detection probability curve (DPC) by using \code{imputeByExpTilt} to impute complete data, then fitting a binomial glm model with the slope as an offset vector.
If the dataset is large, then similar \code{y} values are aggregated before fitting the glm.
}

\value{
A single numeric value giving the intercept.
}

\seealso{
\code{\link{imputeByExpTilt}}.
}

\examples{
y <- simProteinDataSet(n.peptides=100, n.groups=1, dpc.slope=0.7)
estimateDPCIntercept(y, dpc.slope=0.7)
}
