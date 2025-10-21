\name{readFragPipe}
\alias{readFragPipe}
\title{Read Peptide-Precursor Intensities From FragPipe Output}
\description{
Read FragPipe combined peptide, modified peptide or ion output into EList object.
}

\usage{
readFragPipe(
    file = "combined_ion.tsv", path = NULL, sep = "\t", log = TRUE,
    peptide.column = c("Modified Sequence", "Charge"),
    qty.column = NULL,
    qty.column.key = " Intensity", 
    extra.columns = c("Protein", "Protein ID", "Gene", "Protein Description", "Mapped Proteins"),
    match.type.key = NULL, 
    maxlfq = FALSE
)
}

\arguments{
  \item{file}{the name of the file from which the data are to be read, typically should be one of "combined_ion.tsv", "combined_peptide.tsv", or "combined_modified_peptide.tsv" output by FragPipe.}
  \item{path}{character string giving the directory containing the file. Defaults to the current working directory.}
  \item{sep}{the field separator character}
  \item{log}{logical. If \code{TRUE} then intensities will be returned on the log2 scale, otherwise unlogged with zeros.}
  \item{peptide.column}{column containing peptide IDs. Character vector. If length > 1, these columns will be concatenated to make a unique precursor ID.}
  \item{qty.column}{columns containing intensities.}
  \item{qty.column.key}{character string of the key that identify columns containing intensities. If not NULL, \code{qty.column} will be overwritten.}
  \item{extra.columns}{extra columns that are appended to the precursor annotation matrix.}
  \item{match.type.key}{character string of the key that identify columns containing match type annotation. }
  \item{maxlfq}{logical. If \code{TRUE}, columns containing MaxLFQ intensities will be output as \code{E} in the output object. Only set to \code{TRUE} when \code{qty.column.key} identifies the MaxLFQ columns.}
}

\details{
FragPipe (Yu et al 2023) writes a file in wide format called \code{combine_peptide.tsv} containing normalized intensities for peptide precursors from all experimental samples.
\code{readFragPipe} reads this file and produces an EList or EListRaw object.
}

\value{
If \code{log=FALSE}, an EListRaw object containing precursor-level unlogged intensities with zeros and protein annotation.
If \code{log=TRUE}, an EList object containing precursor-level log2 intensities with NAs and protein annotation.
Rows are peptide-precursors and columns are samples.
Peptide precursor and protein annotation is stored in the \code{genes} output component.
}

\references{
Yu F, Teo GC, Kong AT, Fr\"olich K, Li GX, Demichev V, Nesvizhskii AI (2023).
Analysis of DIA proteomics data using MSFragger-DIA and FragPipe computational platform.
\emph{Nature Communications} 14, 4154.
\doi{10.1038/s41467-023-39869-5}
}

\examples{
\dontrun{
ypep <- readFragPipe()
dpcfit <- dpc(ypep)
yprot <- dpcQuant(ypep, dpc=dpcfit)
}
}

\concept{Reading data}
