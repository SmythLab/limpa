\name{filterCompoundProteins}
\alias{filterCompoundProteins}
\alias{filterNonProteotypicPeptides}
\alias{filterSingletonPeptides}
\alias{filterCompoundProteins.default}
\alias{filterNonProteotypicPeptides.default}
\alias{filterSingletonPeptides.default}
\alias{filterCompoundProteins.EList}
\alias{filterNonProteotypicPeptides.EList}
\alias{filterSingletonPeptides.EList}
\alias{filterCompoundProteins.EListRaw}
\alias{filterNonProteotypicPeptides.EListRaw}
\alias{filterSingletonPeptides.EListRaw}

\title{Filtering Based On Protein Annotation}

\description{
Filter peptides or proteins from the dataset based on uniqueness of annotation.
}

\usage{
\method{filterCompoundProteins}{default}(y, protein.group, \dots)
\method{filterCompoundProteins}{EList}(y, protein.group="Protein.Group", \dots)
\method{filterSingletonPeptides}{default}(y, protein.group, min.n.peptides = 2, \dots)
\method{filterSingletonPeptides}{EList}(y, protein.group="Protein.Group", min.n.peptides = 2, \dots)
\method{filterNonProteotypicPeptides}{default}(y, proteotypic, \dots)
\method{filterNonProteotypicPeptides}{EList}(y, proteotypic="Proteotypic", \dots)
}

\arguments{
  \item{y}{a matrix, \code{EList} object or \code{EListRaw} object containing log2-expression values.}
  \item{protein.group}{protein group for each row of \code{y}.
    Can be either a character vector of length \code{nrow(y)} or the name of an annotation column.}
  \item{proteotypic}{indicates whether each peptide is proteotypic (detectable and unique to one protein).
    Should contain 0/1 or TRUE/FALSE values.
    Can be either a vector of length \code{nrow(y)} or the name of an annotation column.}
  \item{min.n.peptides}{minimum number of peptides required in a protein.}
  \item{\dots}{other arguments are not currently used.}
}

\details{
Filter peptide or proteins from the dataset based on uniqueness of annotation.
\code{filterCompoundProteins} removes compound protein groups consisting of multiple proteins separated by ";" delimiters.
\code{filterSingletonPeptides} removes proteins with only one peptide.
\code{filterNonProteotypicPeptides} removes peptides that belong to more than one protein, using the "Proteotypic" annotation column that is returned by DIA-NN and other proteomics quantification software.
}

\value{
An object the same as \code{y} but with non-compliant rows removed.
}

\seealso{
\code{\link{readDIANN}}
}

\concept{Filter features}
