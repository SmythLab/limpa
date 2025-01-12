\name{voomaLmFitWithImputation}
\alias{voomaLmFitWithImputation}

\title{Apply vooma-lmFit Pipeline With Automatic Estimation of Sample Weights and Block Correlation}

\description{
Estimate the variance trend, use it to compute observational weights and use the weights to a fit a linear model.
Includes automatic estimation of sample weights and block correlation.
Equivalent to calling vooma(), arrayWeights(), duplicateCorrelation() and lmFit() iteratively.
}

\usage{
voomaLmFitWithImputation(y, design = NULL, prior.weights = NULL, imputed = NULL,
           block = NULL,
           sample.weights = FALSE, var.design = NULL, var.group = NULL, prior.n = 10,
           predictor = NULL, span = NULL, legacy.span = FALSE,
           plot = FALSE, save.plot = FALSE, keep.EList = TRUE)
}

\arguments{
  \item{y}{
    a numeric \code{matrix}, \code{EList} object, or any object containing log-expression data that can be coerced to a matrix.
    Rows correspond to genes and columns to samples.
  }
  \item{design}{
    design matrix with rows corresponding to samples and columns to coefficients to be estimated.
    Defaults to the unit vector meaning that samples are treated as replicates.
  }
  \item{prior.weights}{
    prior weights.
    Can be a numeric matrix of individual weights of same dimensions as the \code{counts},
    or a numeric vector of sample weights with length equal to \code{ncol(counts)},
    or a numeric vector of gene weights with length equal to \code{nrow(counts)}.
  }
  \item{imputed}{
    logical matrix of the same size as \code{y} indicating whether each observation was entirely imputed.
  }
  \item{block}{
    vector or factor specifying a blocking variable on the arrays.
    Has length equal to \code{ncol(y)}.
  }
  \item{sample.weights}{
    logical value. If \code{TRUE} then empirical sample quality weights will be estimated.
  }
  \item{var.design}{
    design matrix for predicting the sample variances.
    Defaults to the sample-specific model whereby each sample has a different variance.
  }
  \item{var.group}{
    vector or factor indicating groups to have different sample weights.
    This is another way to specify \code{var.design} for groupwise sample weights.
  }
  \item{prior.n}{
    prior number of genes for squeezing the weights towards equality.
    Larger values squeeze the sample weights more strongly towards equality.
  }
  \item{predictor}{
    precision predictor.
    Either a column vector of length \code{nrow(y)} or a numeric matrix of the same dimensions as \code{y} that predicts the precision of each log-expression value.
    Is used as a second covariate together with the log-intensities to predict the variances and produce the final precision weights.
  }
  \item{span}{
    width of the smoothing window, as a proportion of the data set.
    Defaults to a value between 0.3 and 1 that depends the number of genes (\code{nrow(y)}).
    Equal to 1 if the number of genes is less than or equal to 50, then decreases slowly to 0.3 if the number of genes is very large.
  }
  \item{legacy.span}{
    logical.
    If \code{TRUE}, then the original default setting will be used for \code{span}, which is slightly smaller than the new default.
  }
  \item{plot}{
    logical.
    If \code{TRUE}, a plot of the mean-variance trend is displayed.
  }
  \item{save.plot}{
    logical, should the coordinates and line of the plot be saved in the output?
  }
  \item{keep.EList}{
    logical. If \code{TRUE}, then the \code{EList} object containing log-expression values and observation weights will be saved in the component \code{EList} of the output object.
  }
}

\details{
This function is analogous to \code{voomLmFit} in the edgeR package but for microarray-like data with continuous log-expression values.
The function is equivalent to calling \code{vooma()} followed by \code{lmFit()}, optionally with \code{arrayWeights()} and \code{duplicateCorrelation()} as well to estimate sample weights and block correlation.
The function finishes with \code{lmFit()} and returns a fitted model object.

Like \code{vooma}, \code{voomaLmFitWithImputation} estimates the mean-variance relationship in the data and uses it to compute appropriate precision weights for each observation.
The mean-variance trend is estimated from gene-level data but is extrapolated back to individual observations to obtain a precision weight (inverse variance) for each observation.
The weights are then used by \code{lmFit()} to adjust for heteroscedasticity.

Like \code{voomLmFit}, which corrects for loss of residual degrees of freedom due to entirely zero counts in a group (Lun & Smyth 2017), \code{voomaLmFitWithImputation} corrects for loss of residual degrees of freedom due to entirely imputed values in a group (L.
This adjustment prevents from the residual standard deviations from being underestimated due to zero variance between identical imputed values in a group.

If \code{span=NULL}, then an optimal span value is estimated depending on \code{nrow(y)}.
The span is chosen by \code{chooseLowessSpan} with \code{n=nrow(y)}, \code{small.n=50}, \code{min.span=0.3} and \code{power=1/3}.
If \code{legacy.span=TRUE}, then the \code{chooseLowessSpan} arguments are reset to \code{small.n=10}, \code{min.span=0.3} and \code{power=0.5} to match the settings used by \code{vooma} in limma version 3.59.1 and earlier.

If \code{predictor} is not \code{NULL}, then the variance trend is modeled as a function of both the mean log-expression and the \code{predictor} using a multiple linear regression with the two predictors.
In this case, the \code{predictor} is assumed to be some prior predictor of the precision or standard deviation of each log-expression value.
Any \code{predictor} that is correlated with the precision of each observation should give good results.
This ability to model the variance trend using two covariates (mean log-expression and the predictor covariate) was described for the first time by Li (2024).

Sample weights will be estimated using \code{arrayWeights} if \code{sample.weights = TRUE} or if either \code{var.design} or \code{var.group} are non-NULL.
An intra-block correlation will be estimated using \code{duplicateCorrelation} if \code{block} is non-NULL.
In either case, the whole estimation pipeline will be repeated twice to update the sample weights and/or block correlation.
}

\value{
An MArrayLM object containing linear model fits for each row of data.
If sample weights are estimated, then the output object will include a \code{targets} data.frame component with the sample weights as a column with heading \code{"sample.weights"}.

If \code{save.plot=TRUE} then the output object will include components \code{voom.xy} and \code{voom.line}.
\code{voom.xy} contains the x and y coordinates of the points in the vooma variance-trend plot and \code{voom.line} contains the estimated trend line.

If \code{keep.EList=TRUE}, then the output includes component \code{EList} with sub-components \code{Elist$E} and \code{EList$weights}.
If \code{y} was an EList object, then the output \code{EList} preserves all the components of \code{y} and adds the weights.
}

\author{Mengbo Li and Gordon Smyth}

\references{
Li M (2024).
Linear Models and Empirical Bayes Methods for Mass Spectrometry-based Proteomics Data.
PhD Thesis, University of Melbourne.
http://hdl.handle.net/11343/351600

Lun ATL, Smyth GK (2017).
No counts, no variance: allowing for loss of degrees of freedom when assessing biological variability from RNA-seq data.
\emph{Statistical Applications in Genetics and Molecular Biology} 16(2), 83-93.
\doi{10.1515/sagmb-2017-0010}
}

\seealso{
\code{\link{vooma}}, \code{\link{lmFit}}, \code{voomLmFit} (in the edgeR package).
}

\examples{
# Example with a precision predictor
group <- gl(2,4)
design <- model.matrix(~group)
y <- matrix(rnorm(500*8),500,8)
u <- matrix(runif(length(y)),500,8)
yu <- y*u
fit <- voomaLmFitWithImputation(yu,design,plot=TRUE,predictor=u)

# Reproducing vooma plot from output object
fit <- voomaLmFitWithImputation(yu,design,predictor=u,save.plot=TRUE)
do.call(plot,fit$voom.xy)
do.call(lines,fit$voom.line)
}
