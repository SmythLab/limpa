\name{simProteinDataSet}
\alias{simProteinDataSet}
\title{Simulate Peptide Data with NAs By Complete Normal Model}
\description{
Simulate peptide-level log2-expression values from a mass spectrometry experiment.
}
\usage{
simProteinDataSet(n.peptides = 100, n.groups = 2, samples.per.group = 5, peptides.per.protein = 4,
     mu.range = c(2,10), sigma = 0.4, prop.de = 0.2, fc = 2,
     dpc.intercept = NULL, dpc.slope = 0.7, prop.missing = 0.4)
}
\arguments{
  \item{n.peptides}{number of peptides (rows of output).}
  \item{n.groups}{number of experimental groups (conditions).}
  \item{samples.per.group}{number of samples per group.}
  \item{peptides.per.protein}{number of peptides per protein.}
  \item{mu.range}{range of log2-expression values, in terms of expected value per peptide.}
  \item{sigma}{standard deviation of log2-expression values for each peptide in each group.}
  \item{prop.de}{proportion of differentially expressed proteins.}
  \item{fc}{true fold-change for differentially expressed proteins.}
  \item{dpc.intercept}{intercept of detection probability curve. Usually determined from \code{dpc.slope} and \code{prop.missing}.}
  \item{dpc.slope}{slope of detection probability curve.}
  \item{prop.missing}{proportion of missing values (at average log2-expression). Ignored if \code{dpc.intercept} is not \code{NULL}.}
}

\details{
Simulate peptide-level log2-expression values (log2-intensities) from a mass spectrometry experiment.
Values are generated and missing values assigned according to the complete normal model.

Each group of successive peptides is assumed to belong to one protein.
If the protein is differentially expressed (DE), then each peptide belonging to that protein is also DE with the same fold-change.

If \code{dpc.intercept} is not specified, then it is chosen to ensure that the proportion of missing values is equal to \code{prop.missing} at the average log2-expression value.

The simulated data is stored in an EList object, the standard limma package data class for log-expression values.
Peptides are ordered by average expected expression level.
Some of the more lowly expressed peptides may be entirely NA, depending on the argument settings.
}

\value{
EList containing simulated log2-expression values with \code{n.peptides} rows and \code{n.groups * n.samples.per.group} columns.
The EList contains the following components:
  \item{E}{matrix of peptide log2-expression values with NAs.}
  \item{other$E.complete}{matrix of complete log2-expression values without NAs.}
  \item{genes}{data.frame with columns \code{Protein} and \code{DE.Status} giving protein ID and true DE status.}
  \item{targets}{data.frame with column \code{Group} giving group identity for each sample.}
}

\examples{
y <- simProteinDataSet(n.peptides=10, n.groups=1)
show(y)
}

