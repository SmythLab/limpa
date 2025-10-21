\name{plotPeptides}
\alias{plotPeptides}
\alias{plotPeptides.default}
\alias{plotPeptides.EList}

\title{Plot Peptide Log-Intensities for One Protein}

\description{
Plot the peptide-level log-intensities for a specified protein.
}

\usage{
\method{plotPeptides}{default}(y, cex=1.5, lwd=1.5, col="blue", step.down=0.5, ...)
\method{plotPeptides}{EList}(y, index, cex=1.5, lwd=1.5, col="blue", step.down=0.5, ...)
}

\arguments{
  \item{y}{
    a numeric matrix or EList of peptide-level log2-intensity values.
    Columns are samples and rows are peptides or precursors.
  }
  \item{index}{
    index such that \code{y[index,]} will select rows for one protein.
  }
  \item{cex}{
    plot symbol size.
  }
  \item{lwd}{
    line width.
  }
  \item{col}{
    line color.
  }
  \item{step.down}{
    missing values are plotted this amount below the minimum observed value for each peptide.
  }
  \item{\dots}{
    other arguments are passed to \code{plot}.
  }
}

\details{
Plots the log-intensities for the peptides or precursors for one protein, connecting the points belonging to the same peptide.
}

\value{
A plot is made on the current graphics device.
Closed dots correspond to observed values and open dots to missing values.
The matrix of plotting points is also invisibly returned.
}

\seealso{
  \code{\link{plotProtein}}.
}

\examples{
y.peptide <- simProteinDataSet(20,n.groups=1,samples.per.group=4,peptides.per.protein=4)
plotPeptides(y.peptide, 9:12)
}

\concept{Plots}
