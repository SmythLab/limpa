\name{dpc}
\alias{dpc}
\title{Detection Probability Curve Assuming Observed Normal Model}

\description{
Detection probability curve for label free shotgun proteomics data
assuming observed normal log-intensities.
}

\usage{
dpc(y, maxit = 100, eps = 1e-04, b1.upper = 1)
}

\arguments{
  \item{y}{
    numeric matrix of log2-transformed intensities.
    Rows correspond to peptide precursors and columns to samples.
    Any object such as an EList that can be coerced to a matrix is also acceptable.
  }
  \item{maxit}{maximum number of iterations.}
  \item{eps}{convergence tolerance.}
  \item{b1.upper}{upper bound for the DPC slope.}
}

\details{
Estimate the detection probability curve (DPC) for label-free shotgun proteomics data using the method described by Li & Smyth (2023).
This function assumes that the observed log-intensities are normally distributed (the "observed normal" model),
and uses exponential tilting to reformulate the DPC in terms of observed statistics instead of in terms of unobserved quantities.
}

\value{
A list with components
  \item{dpc}{estimated DPC coefficients.}
  \item{history}{iteration history.}
  \item{dpc.start}{initial values estimated for the DPC coefficients.}
  \item{prop.detected}{proportion of observed values for each row.}
  \item{mu.prior}{prior value for row-wise means for observed values.}
  \item{n.prior}{precision of prior for row-wise means, expressed as effective number of observations.}
  \item{s2.prior}{prior value for row-wise variances for observed values.}
  \item{df.prior}{precision of prior for row-wise variances, expressed as prior degrees of freedom.}
  \item{mu.obs}{posterior row-wise means for observed values.}
  \item{s2.obs}{posterior row-wise variances for observed values.}
  \item{mu.mis}{posterior row-wise means for values that are missing.}
}

\references{
Li M, Smyth GK (2023).
Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics.
\emph{Bioinformatics} 39(5), btad200.
\url{10.1093/bioinformatics/btad200}
}

\seealso{
\code{\link{dpcCN}}
}

\examples{
y <- simProteinDataSet(n.peptides=100, n.groups=1)
out <- dpc(y)
out$dpc
}
