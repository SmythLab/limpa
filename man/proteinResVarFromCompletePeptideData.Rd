\name{proteinResVarFromCompletePeptideData}
\alias{proteinResVarFromCompletePeptideData}
\title{Protein Residual Variances From Complete Peptide Data}

\description{
Get protein-wise residual variances by fitting a two-way additive model to the complete (imputed) peptide data for each protein.
}

\usage{
proteinResVarFromCompletePeptideData(y, protein.id, reorder=FALSE)
}

\arguments{
  \item{y}{
    a numeric matrix of complete peptide log2-expression values without NAs.
    Columns are samples and rows are peptides or precursors.
  }    
  \item{protein.id}{
    a character vector of length \code{nrow(y)} giving protein IDs.
  }
  \item{reorder}{
    does the data need to sorted into protein order?
    If \code{TRUE}, then the rows of \code{y} will be sorted so that peptides for the same protein are in consecutive rows.
    If \code{FALSE}, the rows are assumed to be already sorted.
  }
}

\details{
This function operates on complete data after imputation of missing values, and is used to get the \code{sigma} hyperparameters required by \code{peptides2Proteins} and \code{dpcQuant}.
The function fits an additive linear model (~ sample + peptide) to the peptide data for each protein and returns the residual variances.
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
\code{\link{dpcQuant}}, \code{\link{peptides2Proteins}}
}

\examples{
y <- simProteinDataSet(8, n.groups=1, samples.per.group=4, prop.missing=0)
proteinResVarFromCompletePeptideData(y$E, y$genes$Protein)
}

\concept{Quantification}
