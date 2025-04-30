\name{dpcCN}
\alias{dpcCN}
\title{Detection Probability Curve Assuming Complete Normal Model}

\description{
Detection probability curve for label-free shotgun proteomics data
assuming a complete normal model for the peptide intensities.
}

\usage{
dpcCN(y, dpc.start= c(-4,0.7), iterations = 3, verbose = TRUE)
}

\arguments{
  \item{y}{numeric matrix of log2-intensities.
    Rows correspond to peptide precursors and columns to samples.}
  \item{dpc.start}{numeric vector of length 2 giving starting estimates for the DPC intercept and slope.}
  \item{iterations}{number of outer iterations.}
  \item{verbose}{if \code{TRUE}, then progress information will be printed from each iteration.}
}

\details{
Estimate the detection probability curve (DPC) for label-free shotgun proteomics data
by maximum posterior assuming that the complete log-intensities are normally distributed (the "complete normal" model).
The complete log-intensities are the values that would have been observed if the missing value mechanism had not operated.

The algorithm uses an alternating iteration (Smyth, 1996), alternately estimating the row-wise means and standard deviations (mu and sigma) for fixed DPC and estimating the DPC for fixed mu and sigma.
The inner estimations use the BFGS algorithm implemented in the \code{optim} function.
Three outer iterations are usually sufficient.

\code{\link{dpc}} estimates the DPC by a different method, described in Li & Smyth (2023), based on exponential tilting and assuming that only the observed values are normally distributed (the "observed normal" model).
}

\note{
This function may underestimate the DPC slope if entirely missing peptides are omitted and the proportion of peptides that are entirely missing by chance is not small.

\code{dpcCN} can take several minutes on large datasets so, by default, progress information is turned on with \code{verbose=TRUE}.
The function will run quietly if \code{verbose=FALSE} is set.
}

\value{
A list with components
  \item{dpc}{numeric vector of length 2 giving estimated DPC coefficients.}
  \item{mu}{numeric vector of length \code{nrow(y)} giving estimated complete data row-wise means.}
  \item{sigma}{numeric vector of length \code{nrow(y)} giving estimated complete data row-wise standard deviations.}
}

\references{
Li M, Smyth GK (2023).
Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics.
\emph{Bioinformatics} 39(5), btad200.
\url{10.1093/bioinformatics/btad200}

Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}

Li M, Cobbold SA, Smyth GK (2025).
Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values.
\emph{bioRxiv} 2025/651125.
\doi{10.1101/2025.04.28.651125}

Smyth GK (1996).
Partitioned algorithms for maximum likelihood and other non-linear estimation.
\emph{Statistics and Computing} 6, 201-216.
\doi{10.1007/BF00140865}
\url{https://gksmyth.github.io/pubs/partitio.pdf}
}

\seealso{
\code{\link{dpc}}.
}

\examples{
y <- simProteinDataSet(n.peptides=100, n.groups=1)
out <- dpcCN(y)
out$dpc
}
