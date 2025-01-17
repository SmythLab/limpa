\name{plotProtein}
\alias{plotProtein}
\title{Plot protein summary with error bars by DPC-Quant}

\description{
Plot the log-intensity of a protein summarized by DPC-Quant for each sample with error bars. 
}

\usage{
plotProtein(y, protein, col = "black", cex = 2, lwd = 2, ...)
}

\arguments{
  \item{y}{protein-level EList produced by \code{dpcQuant()}.}
  \item{protein}{
    A vector of length 1. Can be the name of the protein or the numeric index that locates
    the protein to plot from rows of \code{y}. 
  }
  \item{col}{Color for the points and error bars. }
  \item{cex}{Size for the points. }
  \item{lwd}{Line width for the error bars. }
  \item{\dots}{other arguments are passed to \code{plot()}.}
}

\details{
Plot the sample-wise protein quantification results from \code{dpcQuant()} for each protein. 
The error bars (standard errors) indicate the quantification uncertainty associated with each estimate. 
Typically within a dataset, the larger the error bar is, the more missing values there are in the 
precursor/peptide-level data for that protein. 
}

\value{
A plot is created on the current graphics device.
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}
}

\examples{
y.peptide <- simProteinDataSet()
y.protein <- dpcQuant(y.peptide, "Protein", dpc=c(-4,0.7))
plotProtein(y.protein, protein = "Protein01", col = rep(c("blue", "red"), each = 5))
}