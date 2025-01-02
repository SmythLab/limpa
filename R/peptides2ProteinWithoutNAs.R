peptides2ProteinWithoutNAs <- function(y, sigma=0.5, weights=NULL, dpc=c(-4,0.7), prior.mean=6, prior.sd=10, prior.logFC=2)
# Summarize peptide to protein log-expression for one protein when there are no NAs.
# This function is not particularly fast as currently written; peptides2ProteinWithDerivs() is faster.
# Created 15 July 2023. Last modified 15 July 2023.
{
# Check y
  if(anyNA(y)) stop("NAs not allowed")

# DPC parameters
  dpc.intercept <- dpc[1]
  dpc.slope <- dpc[2]

# Additive model for peptides and samples
  npeptides <- nrow(y)
  nsamples <- ncol(y)
  i1n <- seq_len(nsamples)
  Samples <- factor(rep.int(seq_len(nsamples),rep.int(npeptides,nsamples)))
  if(npeptides > 1L) {
    Peptides <- factor(rep.int(seq_len(npeptides),nsamples))
    contrasts(Peptides) <- contr.sum(npeptides)
    X <- model.matrix(~0+Samples+Peptides)
  } else {
    X <- model.matrix(~0+Samples)
  }

# Check weights
  if(is.null(weights)) {
    weights <- matrix(1,npeptides,nsamples)
  } else {
    if(!identical(dim(y),dim(weights))) stop("weights must be conformal with y")
  }

# Convert y to vector
  SampleNames <- colnames(y)
  if(is.null(SampleNames)) SampleNames <- paste0("Sample", formatC(i1n,format="d",flag="0",width=1L+floor(log10(nsamples))))
  y <- as.vector(y)
  w <- as.vector(weights)/sigma^2

# Prior covariance matrix
  s2_1 <- prior.sd^2 - prior.logFC^2 / nsamples
  s2_2 <- prior.logFC^2
  prior.Sigma <- matrix(s2_1,nsamples,nsamples)
  diag(prior.Sigma) <- s2_1 + s2_2
  prior.Sigma.inverse <- solve(prior.Sigma)

#  a <- 1/prior.logFC^2
#  b <- 1/(nsamples*prior.sd)^2 - a/nsamples
#  prior.Sigma.inverse <- matrix(b,nsamples,nsamples)
#  diag(prior.Sigma.inverse) <- diag(prior.Sigma.inverse) + a

# Ridge regression
  XtX <- crossprod(X,w*X)
  XtX[i1n,i1n] <- XtX[i1n,i1n] + prior.Sigma.inverse
  Xty <- crossprod(X,w*y)
  Xty[i1n] <- Xty[i1n] + rowSums(prior.Sigma.inverse) * prior.mean
  drop(solve(XtX, Xty))[i1n]
}
