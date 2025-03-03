\name{readDIANN}
\alias{readDIANN}
\title{Read Peptide-Precursor Intensities From DIA-NN Output}
\description{
Read DIA-NN Reports.tsv file into EList object.
}

\usage{
readDIANN(file = "Report.tsv", path = NULL, format="tsv", sep = "\t", log = TRUE,
        q.columns = c("Global.Q.Value", "Lib.Q.Value"), q.cutoffs = c(0.01, 0.01))
}

\arguments{
  \item{file}{the name of the file from which the data are to be read.}
  \item{path}{character string giving the directory containing the file. Defaults to the current working directory.}
  \item{format}{character string giving the format of the file. Possible values are "tsv" or "parquet". Default is "tsv". }
  \item{sep}{the field separator character}
  \item{log}{logical. If \code{TRUE} then intensities will be returned on the log2 scale, otherwise unlogged with zeros.}
  \item{q.columns}{column headings in the DIA-NN output containing Q-values for peptide identification. Character vector.}
  \item{q.cutoffs}{cutoffs to apply to the Q-value columns. Only peptides with values below the cutoffs will be retained. Numeric vector of same length as \code{q.columns}.}
}

\details{
DIA-NN (Demichev et al 2020) writes a file in long (data.frame) format, typically called \code{Report.tsv}, containing normalized intensities for peptide precursors.
\code{readDIANN} reads this file and produces an object in limma EList or EListRaw format.
}

\value{
If \code{log=FALSE}, an EListRaw object containing precursor-level unlogged intensities with zeros and protein annotation.
If \code{log=TRUE}, an EList object containing precursor-level log2 intensities with NAs and protein annotation.
Rows are peptide-precursors and columns are samples.
Peptide precursor and protein annotation is stored in the `genes` output component.
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

\keyword{reading data}
