\name{readMaxQuant}
\alias{readMaxQuant}
\title{Read Peptide-Precursor Intensities From MaxQuant Output}
\description{
Read MaxQuant peptide output into EList object.
}

\usage{
readMaxQuant(
    file = "peptides.txt", 
    path = NULL, sep = "\t", log = TRUE,
    peptide.column = "Sequence",
    qty.column = NULL,
    qty.column.key = "Intensity ", 
    q.columns = c("PEP"), q.cutoffs = 0.01,
    extra.columns = c("Proteins", "Leading razor protein", 
    "Gene names", "Unique (Groups)", "Protein group IDs")
)
}

\arguments{
  \item{file}{the name of the file from which the data are to be read.}
  \item{path}{character string giving the directory containing the file. Defaults to the current working directory.}
  \item{sep}{the field separator character}
  \item{log}{logical. If \code{TRUE} then intensities will be returned on the log2 scale, otherwise unlogged with zeros.}
  \item{peptide.column}{column containing peptide IDs. String of length 1L.}
  \item{qty.column}{columns containing intensities.}
  \item{qty.column.key}{character string of the key that identify columns containing intensities. If not NULL, \code{qty.column} will be overwritten.}
  \item{q.columns}{column headings in the output containing q-values for peptide identification. Character vector.}
  \item{q.cutoffs}{cutoffs to apply to the q-value columns. Only peptides with values below the cutoffs will be retained. Numeric vector of same length as \code{q.columns}.}
  \item{extra.columns}{extra columns that are appended to the precursor annotation matrix.}
}

\details{
\code{readMaxQuant} reads standard peptide.txt output from MaxQuant and produces an object in limma EList or EListRaw format.
}

\value{
If \code{log=FALSE}, an EListRaw object containing precursor-level unlogged intensities with zeros and protein annotation.
If \code{log=TRUE}, an EList object containing precursor-level log2 intensities with NAs and protein annotation.
Rows are peptide-precursors and columns are samples.
Peptide precursor and protein annotation is stored in the \code{genes} output component.
}

\examples{
\dontrun{
y <- readMaxQuant("peptides.txt")
dpcfit <- dpc(y)
}
}

\concept{Reading data}
