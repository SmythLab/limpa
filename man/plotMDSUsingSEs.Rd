\title{Multidimensional Scaling Plot of Gene Expression Profiles, Using Standard Errors}
\name{plotMDSUsingSEs}
\alias{plotMDSUsingSEs}

\description{
Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical z-statistic of differences between the samples.
}

\usage{
plotMDSUsingSEs(y, top = 500, labels = NULL, pch = NULL, cex = 1, dim.plot = c(1,2),
     xlab = NULL, ylab = NULL, plot = TRUE, var.explained = TRUE, \dots)
}

\arguments{
  \item{y}{\code{EList} produced by \code{dpcQuant} or \code{dpcImpute}.}
  \item{top}{number of top genes used to calculate pairwise distances.}
  \item{labels}{character vector of sample names or labels. Defaults to \code{colnames(x)}.}
  \item{pch}{plotting symbol or symbols. See \code{\link{points}} for possible values. Ignored if \code{labels} is non-\code{NULL}.}
  \item{cex}{numeric vector of plot symbol expansions.}
  \item{dim.plot}{integer vector of length two specifying which principal components should be plotted.}
  \item{xlab}{title for the x-axis.}
  \item{ylab}{title for the y-axis.}
  \item{plot}{logical. If \code{TRUE} then a plot is created on the current graphics device.}
  \item{var.explained}{logical. If \code{TRUE} then the percentage variation explained is included in the axis labels.}
  \item{\dots}{any other arguments are passed to \code{plot}, and also to \code{text} (if \code{pch} is \code{NULL}).}
}

\details{
This function uses multidimensional scaling (MDS) to produce a principal coordinate (PCoA) plot showing the relationships between the expression profiles represented by the columns of \code{x}.
Distances on the plot represent the \emph{leading z-statistic}.
The leading log-fold-change between a pair of samples is defined as the root-mean-square average of the \code{top} largest z-statistics between those two samples.

If \code{pch=NULL}, then each sample is represented by a text label, defaulting to the column names of \code{x}.
If \code{pch} is not \code{NULL}, then plotting symbols are used.

See \code{\link[graphics]{text}} for possible values for \code{col} and \code{cex}.
}

\value{
If \code{plot=TRUE} or if \code{x} is an object of class \code{"MDS"}, then a plot is created on the current graphics device.

An object of class \code{"MDS"} is also invisibly returned.
This is a list containing the following components:
\item{eigen.values}{eigen values}
\item{eigen.vectors}{eigen vectors}
\item{var.explained}{proportion of variance explained by each dimension}
\item{distance.matrix.squared}{numeric matrix of squared pairwise distances between columns of \code{x}}
\item{dim.plot}{dimensions plotted}
\item{x}{x-xordinates of plotted points}
\item{y}{y-cordinates of plotted points}
\item{gene.selection}{gene selection method}
}

\author{Gordon Smyth}

\references{
Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, and Smyth GK (2015).
limma powers differential expression analyses for RNA-sequencing and microarray studies.
\emph{Nucleic Acids Research} 43, e47.
\url{http://nar.oxfordjournals.org/content/43/7/e47}
}

\seealso{
\code{plotMDS} in the limma package.
}

\examples{
# See dpcQuant()
}

\concept{plots}
