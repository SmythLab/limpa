\name{peptides2Proteins}
\alias{peptides2Proteins}
\title{DPC-Quant for Many Proteins}
\description{
Quantify protein expression values by the DPC-Quant method.
}

\usage{
peptides2Proteins(y, protein.id, sigma = 0.5, dpc = c(-4, 0.7),
        prior.mean = 6, prior.sd = 10, prior.logFC = 2,
        standard.errors = FALSE, newton.polish = FALSE, verbose = FALSE, chunk = 1000L)
}

\arguments{
  \item{y}{
    a numeric matrix of log-expression values.
    Columns are samples and rows are peptides or precursors.
  }
  \item{protein.id}{
    protein IDs.
    Character vector of length \code{nrow(y)}.
  }
  \item{sigma}{
    standard deviations of peptide-level expression values.
    Numeric vector of same length as \code{unique(protein.id))}.
   }
  \item{dpc}{
    numeric vector giving intercept and slope of detection probability curve (DPC).
  }
  \item{prior.mean}{
    mean of the global prior distribution for protein log-expression values.
    Represents the typical average log-expression of a protein.
  }
  \item{prior.sd}{
    standard deviation of the global prior distribution for protein log-expression values.
    Represents the standard deviation of average log-expression across proteins.
  }
  \item{prior.logFC}{
    standard deviation to be expected between log-expression values for the same protein.
  }
  \item{standard.errors}{
    logical, should standard errors for the protein expression values be returned?
  }
  \item{newton.polish}{
    logical.
    If \code{TRUE} then one Newton iteration will be done to refine the optimization after the BFGS algorithm has finished.
    Ignored if \code{standard.errors=FALSE}.
  }
  \item{verbose}{
    should progress information be output?
    If \code{TRUE}, then progress information is output every \code{chunk} proteins.
  }
  \item{chunk}{
    When \code{verbose=TRUE}, how often to output progress information.
    By default, reports every 1000 proteins.
  }
}

\details{
Implements the DPC-Quant method, which returns maximum posterior estimates for protein expression values.
}

\value{
An \code{EList} object with a row for each protein.
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

\examples{
y.peptide <- simProteinDataSet(8,n.groups=1,samples.per.group=4,prop.missing=0.2)
y.protein <- peptides2Proteins(y.peptide$E, y.peptide$genes$Protein)
}

\concept{Quantification}
