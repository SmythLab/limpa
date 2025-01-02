\name{dpcQuantHyperparam}
\alias{dpcQuantHyperparam}
\title{Estimate Hyperparameters for DPC-Quant}

\description{
Estimate hyperparameters for the DPC-based protein quantification method (DPC-Quant). 
}

\usage{
dpcQuantHyperparam(y, protein.id, dpc.slope = 0.7, sd.quantile.for.logFC = 0.9, robust = FALSE, ...)
}

\arguments{
  \item{y}{
    a numeric matrix of peptide-level log2-expression values.
    Columns are samples and rows are peptides or precursors.
  }
  \item{protein.id}{
    a character vector of length \code{nrow(y)} giving protein IDs.
  }
  \item{dpc.slope}{
    slope of the DPC.
  }
  \item{sd.quantile.for.logFC}{
    a number between 0 and 1.
    The quantile of the precursor-level variances to represent the typical between-sample variation. 
  }
  \item{robust}{
    should robust empirical Bayes moderation be applied to the protein standard deviations?
    \code{robust=TRUE} will cause very large standard deviations to be squeezed less strongly towards the prior value.
  }
  \item{\dots}{
    other arguments are passed to \code{imputeByExpTilt}.
  }
}

\details{
Estimates and returns the empirical Bayes hyperparameters required for DPC-Quant protein quantification. 
This function is called by the \code{dpcQuant} function. 
}

\value{
A list with components
  \item{prior.mean}{
    mean of the global prior distribution for protein log-expression values.
    Represents the typical average log-expression of a protein.
  }
  \item{prior.sd}{
    standard deviation of the global prior distribution for protein log-expression values.
    Represents the standard deviation of average log-expression across proteins.
  }
  \item{prior.logFC}{
    standard deviation to be expected between log-expression values for the same protein across conditions.
  }
  \item{sigma}{
    protein standard deviations from additive model fitted to peptide log expression values.
    Numeric vector of same length as \code{unique(protein.id))}.
  }
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}
}

\seealso{
\code{\link{dpcQuant}}, \code{\link{imputeByExpTilt}}
}
