\name{readDIANN}
\alias{readDIANN}
\title{Read Peptide-Precursor Intensities From DIA-NN Output}
\description{
Read DIA-NN Reports.tsv file into EList object.
}

\usage{
readDIANN(file = "Report.tsv", path = NULL, format = "parquet", sep = "\t", log = TRUE,
  run.column = "Run",
  precursor.column = "Precursor.Id",
  qty.column = "Precursor.Normalised",
  q.columns = c("Global.Q.Value", "Lib.Q.Value"),
  q.cutoffs = c(0.01, 0.01), 
  extra.columns = c("Protein.Group", "Protein.Names", "Genes", "Proteotypic"))
}

\arguments{
  \item{file}{the name of the file from which the data are to be read. Or it can also be the data.frame read from the report in the long format, where each row is an observation.}
  \item{path}{character string giving the directory containing the file. Defaults to the current working directory.}
  \item{format}{character string giving the format of the file. Possible values are "tsv" or "parquet". Default is "parquet".}
  \item{sep}{the field separator character}
  \item{log}{logical. If \code{TRUE} then intensities will be returned on the log2 scale, otherwise unlogged with zeros.}
  \item{run.column}{character string giving the name of the column containing run names.}
  \item{precursor.column}{character string giving the name of the column containing precursor IDs.}
  \item{qty.column}{character string giving the name of the column containing precursor intensities.}
  \item{q.columns}{column headings in the DIA-NN output containing Q-values for peptide identification. Character vector.}
  \item{q.cutoffs}{cutoffs to apply to the Q-value columns. Only peptides with values below the cutoffs will be retained. Numeric vector of same length as \code{q.columns}.}
  \item{extra.columns}{extra columns that are appended to the precursor annotation matrix.}
}

\details{
DIA-NN (Demichev et al 2020) writes a file in long (data.frame) format, typically called \code{Report.tsv}, containing normalized intensities for peptide precursors.
\code{readDIANN} reads this file and produces an EList or EListRaw object.

From version 2.0, DIA-NN returns the main report in Apache Parquet format (\url{https://github.com/vdemichev/DiaNN/releases}). 
\code{readDIANN} can read the Parquet file directly or, alternatively, one can read the Parquet file into a data.frame, and use \code{readDIANN} to process the long-format data.frame into a limma EList or EListRaw object. 
}

\value{
If \code{log=FALSE}, an EListRaw object containing precursor-level unlogged intensities with zeros and protein annotation.
If \code{log=TRUE}, an EList object containing precursor-level log2 intensities with NAs and protein annotation.
Rows are peptide-precursors and columns are samples.
Peptide precursor and protein annotation is stored in the \code{genes} output component.
}

\references{
Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M (2020).
DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput.
\emph{Nature Methods} 17(1), 41-44.
}

\examples{
\dontrun{
ypep <- readDIAN()
ypep <- filterCompoundProteins(ypep)
ypep <- filterNonProteotypicPeptides(ypep)
dpcfit <- dpc(ypep)
yprot <- dpcQuant(ypep, dpc=dpcfit)
}
}

\concept{Reading data}
