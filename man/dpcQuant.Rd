\name{dpcQuant}
\alias{dpcQuant}
\alias{dpcQuant.default}
\alias{dpcQuant.EList}
\title{Quantify Proteins}

\description{
Use the DPC to quantify protein expression values.
}

\usage{
\method{dpcQuant}{EList}(y, protein.id = "Protein.Group", dpc = NULL, verbose = TRUE, ...)
}

\arguments{
  \item{y}{
    a numeric matrix or EList of peptide-level log2-expression values.
    Columns are samples and rows are peptides or precursors.
  }
  \item{protein.id}{
    protein IDs.
    Either an annotation column name (if \code{y} is an EList) or a character vector of length \code{nrow(y)}.
  }
  \item{dpc}{
    numeric vector giving intercept and slope of DPC.
    Alternatively the output objects from \code{dpc} or \code{dpcCN} are also acceptable.
  }
  \item{verbose}{
    should progress information be output?
    If \code{TRUE}, then progress information is output every 1000 proteins.
  }
  \item{\dots}{
    other arguments are passed to \code{dpcQuantHyparam}.
  }
}

\details{
Implements the DPC-Quant method, which quantifies protein log2-expression values from peptide data.
The method represents missing values probabilistically using the PDC and returns maximum posterior estimates for all the protein log2-expression values, so that there are no missing values in the final summary.

The \code{dpc} function should be run before running this function.
The output from this function can be input to \code{dpcDE}.
}

\value{
An \code{EList} object with a row for each protein, with the following extra components:
  \item{other$n.observations}{matrix giving the number of missing non-missing peptide observations supporting each protein expression value.}
  \item{other$standard.error}{matrix giving the standard error of each protein expression value.}
}

\note{
\code{dpcQuant} can take several minutes on large datasets so, by default, progress information is turned on with \code{verbose=TRUE}.
The function will run quietly if \code{verbose=FALSE} is set.
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}
}

\seealso{
  \code{\link{dpc}}, \code{\link{dpcQuantHyperparam}}, \code{\link{dpcDE}}, \code{\link{EList-class}}.
}

\examples{
y.peptide <- simProteinDataSet(n.groups=1,samples.per.group=4,prop.missing=0.2)
y.protein <- dpcQuant(y.peptide, "Protein", dpc=c(-4,0.7))
}
