\name{imputeByExpTilt}
\alias{imputeByExpTilt}
\alias{imputeByExpTilt.default}
\alias{imputeByExpTilt.EListRaw}
\alias{imputeByExpTilt.EList}
\alias{expTiltByRows}
\alias{expTiltByColumns}
\title{Impute Missing Values by Exponential Tilting}
\description{
Impute missing values in a log-expression matrix by applying exponential tilting to rows, columns or both.
}

\usage{
\method{imputeByExpTilt}{default}(y, dpc.slope = 0.7, prior.logfc = NULL, by = "both", \dots)
expTiltByRows(y, dpc.slope = 0.7, sigma.obs = NULL)
expTiltByColumns(y, dpc.slope = 0.7)
}

\arguments{
  \item{y}{an \code{EList} object or a numeric matrix of log-expression values. Columns are samples and rows are peptides or proteins.
  For \code{expTiltByRows} or \code{expTiltByRows}, should be a numeric matrix.}
  \item{dpc.slope}{slope of detection probability curve.}
  \item{prior.logfc, sigma.obs}{simple standard deviation to be expected between observed values for the same peptide or protein. Can a single value or vector of length \code{nrow(y)}. By default is estimated from the data.}
  \item{by}{character value. Should imputation by rows (\code{"rows"}), by columns (\code{"columns"}) or both (\code{"both"})?}
  \item{\dots}{other arguments are not used.}
}

\details{
Implements exponential tilting strategy outlined by Li & Smyth (2023).
The imputed values are the expected values of the missing value distribution.

The strategy can be applied to rows or columns.
If \code{by="both"}, the imputated values are an average of the row and column imputations, weighted inversely by the prediction variances.
}

\value{
An object of the same class as \code{y} but with NAs imputed.
}

\references{
Li M, Smyth GK (2023).
Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics.
\emph{Bioinformatics} 39(5), btad200.
\doi{10.1093/bioinformatics/btad200}
}

\examples{
y <- matrix(rnorm(25),5,5)
y[1,1] <- NA
imputeByExpTilt(y)
}

\concept{Quantification}
