\name{readSpectronaut}
\alias{readSpectronaut}
\title{Read Peptide-Precursor Intensities From Spectronaut Output}
\description{
Read Spectronaut Reports.tsv file into EList object.
}

\usage{
readSpectronaut(file = "Report.tsv", path = NULL, sep = "\t", log = TRUE,
                run.column = "R.Raw File Name", qty.column = "EG.TotalQuantity",
                q.columns = c("EG.Qvalue", "PG.Qvalue"), q.cutoffs = 0.01)
}

\arguments{
  \item{file}{the name of the file from which the data are to be read.}
  \item{path}{character string giving the directory containing the file. Defaults to the current working directory.}
  \item{sep}{the field separator character}
  \item{log}{logical. If \code{TRUE} then intensities will be returned on the log2 scale, otherwise unlogged with zeros.}
  \item{run.column}{column containing run}
  \item{qty.column}{column containing intensities}
  \item{q.columns}{column headings in the Spectronaut output containing Q-values for peptide identification. Character vector.}
  \item{q.cutoffs}{cutoffs to apply to the Q-value columns. Only peptides with values below the cutoffs will be retained. Numeric vector of same length as \code{q.columns}.}
}

\details{
Spectronaut (https://biognosys.com/software/spectronaut/) writes a file in long (data.frame) format, typically called \code{Report.tsv}, containing normalized intensities for peptide precursors.
\code{readSpectronaut} reads this file and produces an object in limma EList or EListRaw format.
}

\value{
If \code{log=FALSE}, an EListRaw object containing precursor-level unlogged intensities with zeros and protein annotation.
If \code{log=TRUE}, an EList object containing precursor-level log2 intensities with NAs and protein annotation.
Rows are peptide-precursors and columns are samples.
Peptide precursor and protein annotation is stored in the `genes` output component.
}

\examples{
\dontrun{
y <- readSpectronaut()
dpcfit <- dpc(y)
}
}

\keyword{reading data}
