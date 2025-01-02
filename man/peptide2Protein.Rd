\name{peptides2ProteinBFGS}
\alias{peptides2ProteinBFGS}
\alias{peptides2ProteinNewton}
\alias{peptides2ProteinWithoutNAs}
\title{DPC-Quant for One Protein}
\description{
Convert a matrix of peptide log-expression values for one protein to protein-level expression values by the DPC-Quant method.
}

\usage{
peptides2ProteinBFGS(y, sigma = 0.5, weights = NULL, dpc = c(-4, 0.7),
        prior.mean = 6, prior.sd = 10, prior.logFC = 2,
        standard.errors = TRUE, newton.polish = TRUE, start = NULL)
peptides2ProteinNewton(y, sigma = 0.5, weights = NULL, dpc = c(-4, 0.7),
        prior.mean = 6, prior.sd = 10, prior.logFC = 2,
        standard.errors = TRUE, tol=1e-6, maxit=10, start = NULL, verbose = FALSE)
peptides2ProteinWithoutNAs(y, sigma = 0.5, weights = NULL, dpc = c(-4, 0.7),
        prior.mean = 6, prior.sd = 10, prior.logFC = 2)
}

\arguments{
  \item{y}{
    a numeric matrix of log-expression values.
    Columns are samples and rows are peptides or precursors.
    Typically contains NAs, but NAs are not allowed for \code{peptides2ProteinWithoutNAs}.
  }
  \item{sigma}{
    standard deviation of peptide-level expression values after allowing for peptide and sample baseline differences.
  }
  \item{weights}{
    numeric matrix of same size as \code{y} containing positive precision weights.
    The precision of the log-expression values is summarized by \code{sigma/sqrt(weights)}.
  }
  \item{dpc}{
    numeric vector giving intercept and slope of the detection probability curve (DPC).
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
  \item{start}{
    numeric vector of starting values for the linear model coefficients.
    Of length \code{ncol(y)+nrow(y)-1)}.
  }
  \item{tol}{
    stopping criterion tolerance for Newton's method, to be achieved by the average local slope statistic.
  }
  \item{maxit}{
    maximum number of iterations for Newton's method.
  }
  \item{verbose}{
    logical.
    If \code{TRUE}, progress will be output at each iteration.
  }
}

\details{
Implements the DPC-Quant method, which returns maximum posterior estimates for protein expression values.

\code{peptides2ProteinBFGS} maximizes the posterior using the BFGS algorithm with analytic first derivatives.
The standard errors are computed from analytic second derivatives.

\code{peptides2ProteinNewton} maximizes the posterior using Newton's method.
}

\value{
\code{peptides2ProteinBFGS} and \code{peptides2ProteinNewton} return a list with components.
  \item{protein.expression}{numeric vector giving the estimated protein log-expression value for each sample.}
  \item{standard.error}{numeric vector giving standard errors for the protein log-expression values.}
  \item{value}{the minimized objective function, minus twice the log-posterior distribution.}
  
\code{peptides2ProteinWithoutNAs} returns a numeric vector of protein expression values.
}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
\url{http://hdl.handle.net/11343/351600}

Smyth GK (2005).
Optimization and nonlinear equations.
In: \emph{Encyclopedia of Biostatistics Second Edition}, Volume 6, P. Armitage and T. Colton (eds.), Wiley, London, pages 3857-3863.
\url{https://gksmyth.github.io/pubs/OptimNonlinEqnPreprint.pdf}
}

\examples{
y <- matrix(rnorm(12),3,4)
y[1:2,1] <- NA
y[1,2] <- NA
peptides2ProteinBFGS(y)
}
