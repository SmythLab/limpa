\name{dpcDE}
\alias{dpcDE}
\title{Fit Linear Model With Precision Weights}
\description{
Fit linear models and make precision weights from the DPC-Quant standard errors.
}

\usage{
dpcDE(y, design, plot=TRUE, ...)
}

\arguments{
  \item{y}{protein-level EList produced by dpcQuant().}
  \item{design}{design matrix.}
  \item{plot}{should the variance trend be plotted?}
  \item{\dots}{other arguments are passed to \code{voomaLmFitWithImputation}.}
}

\details{
Calls \code{\link{voomaLmFitWithImputation}} to compute vooma precision weights from the DPC-Quant standard errors stored in \code{y} and to use those weights to fit protein-wise linear models.
Any \code{voomaLmFit} functionality can be used, giving access to optional empirical sample weights or random blocks.
}

\value{
An \code{MArrayLM} object suitable for analysis in limma.
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}

Li M, Cobbold SA, Smyth GK (2025).
Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values.
\emph{bioRxiv} 2025/651125.
\doi{10.1101/2025.04.28.651125}
}

\seealso{
  \code{\link{voomaLmFitWithImputation}}. Also \code{voomaLmFit} in the limma package.
}

\examples{
y.peptide <- simProteinDataSet()
y.protein <- dpcQuant(y.peptide, "Protein", dpc=c(-4,0.7))
Group <- factor(y.peptide$targets$Group)
design <- model.matrix(~Group)
fit <- dpcDE(y.protein, design)
}
