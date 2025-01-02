peptides2ProteinNewton <- function(y, sigma=0.5, weights=NULL, dpc=c(-4,0.7), prior.mean=6, prior.sd=10, prior.logFC=2, standard.errors=TRUE, tol=1e-6, maxit=10, start=NULL, verbose=FALSE)
# Summarize peptide to protein log-expression for one protein
# assuming a complete normal model, additive linear model,
# logit-linear detection probabilities and
# and a normal prior distribution for protein log-expression.
# Uses Newton's method.
# Input:
# y: numeric matrix of log-intensities, rows are peptides, columns are samples
# sigma: complete data residual standard deviation for the additive linear model
# Created 18 May 2024. Last modified 1 Jan 2025.
{
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
  N <- nrow(X)
  nbeta <- ncol(X)

# Check weights
  if(is.null(weights)) {
    weights <- matrix(1,npeptides,nsamples)
  } else {
    if(!identical(dim(y),dim(weights))) stop("weights must be conformal with y")
    if(min(weights) < 1e-15) stop("weights must be positive")
  }

# Starting values
  if(is.null(start)) {
    beta <- rep_len(0,nbeta)
    yimp <- suppressWarnings(imputeByExpTilt(y))
    beta[i1n] <- colMeans(yimp)
    if(npeptides > 1L) {
      b <- rowMeans(yimp)
      b <- b-mean(b)
      beta[(nsamples+1):nbeta] <- b[-npeptides]
    }
  } else {
    if(!identical(length(start),nbeta)) stop("`start` is not correct length")
    beta <- start
  }

# Convert y to vector
  SampleNames <- colnames(y)
  if(is.null(SampleNames)) SampleNames <- paste0("Sample", formatC(i1n,format="d",flag="0",width=1L+floor(log10(nsamples))))
  y <- as.vector(y)
  sigma <- sigma/sqrt(as.vector(weights))

# Subset observed and missing values
  mis <- is.na(y)
  obs <- !mis
  mis <- which(mis)
  obs <- which(obs)
  nobs <- length(obs)
  nmis <- length(mis)
  yo <- y[obs]
  sigmao <- sigma[obs]
  sigmam <- sigma[mis]
  Xo <- X[obs,,drop=FALSE]
  Xm <- X[mis,,drop=FALSE]

# Gauss quadrature
  if(nmis) {
    gq <- gauss.quad.prob(16,dist="normal")
    z <- matrix(gq$nodes,nmis,16,byrow=TRUE)
  }

  Q <- Inf
  for (iter in seq_len(maxit)) {

#   Minus twice log-posterior
    mu <- X %*% beta
#   Observed part
    if(length(obs)) {
      muo <- mu[obs]
      logposto <- sum(((yo-muo)/sigmao)^2)
    } else {
      logposto <- 0
    }
#   Missing part
    if(length(mis)) {
      mum <- mu[mis]
      eta <- dpc.intercept + dpc.slope * (mum + sigmam * z)
      G <- plogis(eta,lower.tail=FALSE) %*% gq$weights
      logpostm <- -2*sum(log(G))
    } else {
      logpostm <- 0
    }
#   Prior part
    ProteinExpr <- beta[seq_len(nsamples)]
    ProteinExpr.Mean <- mean(ProteinExpr)
    logpostp <- (ProteinExpr.Mean - prior.mean)^2 / prior.sd^2
    logpostp <- logpostp + sum( (ProteinExpr - ProteinExpr.Mean)^2 ) / prior.logFC^2
#   Total
    Qold <- Q
    Q <- logposto + logpostm + logpostp
    if(verbose) message("Q ",Q)
    if(Q > Qold) warning("Divergence, likelihood increased")

#   Derivative
    DQ.Dmu <- rep_len(0,N)
#   Observed part
    if(length(obs)) {
      muo <- mu[obs]
      DQ.Dmu[obs] <- -(yo-muo)/sigmao^2
    }
#   Missing part
    if(length(mis)) {
      Gdot <- -dpc.slope * dlogis(eta) %*% gq$weights
      DQ.Dmu[mis] <- -Gdot/G
    }
    DQ <- crossprod(X,DQ.Dmu)
#   Prior part
    ProteinExpr <- beta[i1n]
    ProteinExpr.Mean <- mean(ProteinExpr)
    DQ.Prior <- (ProteinExpr.Mean-prior.mean)/nsamples/prior.sd^2 + (ProteinExpr-ProteinExpr.Mean)/prior.logFC^2
    DQ[i1n] <- DQ[i1n] + DQ.Prior
#   Total
    DQ <- as.vector(2*DQ)

#   Second derivative
    D2Q.Dmu <- rep_len(0,N)
    D2Q.Dmu[obs] <- 1/sigma[obs]^2
    if(length(mis)) {
      Gddot <- dpc.slope^2/4 * (tanh(eta/2) / cosh(eta/2)^2) %*% gq$weights
      D2Q.Dmu[mis] <- (Gdot*Gdot - G*Gddot)/G/G
    }
    D2Q <- crossprod(X,D2Q.Dmu*X)
    D2Q.Prior <- 1/nsamples^2/prior.sd^2 - 1/nsamples/prior.logFC^2 + diag(nsamples)/prior.logFC^2
    D2Q[i1n,i1n] <- D2Q[i1n,i1n] + D2Q.Prior
    D2Q <- 2*D2Q

#   R <- chol(D2Q)
#   delta <- backsolve(R,backsolve(R,DQ,transpose=TRUE),transpose=FALSE)
    delta <- solve(D2Q,DQ)
    beta <- beta - delta
    eps <- crossprod(DQ,delta)/N
    if(verbose) message("eps ", eps)
    if(eps < tol) break
  }

# Standard errors
  if(standard.errors)
    standard.error <- sqrt(2*diag(chol2inv(chol(D2Q)))[i1n])
  else
    standard.error <- NULL

# Output
  list(protein.expression=beta[i1n],value=Q,deriv=DQ,hessian=D2Q,standard.error=standard.error,iter=iter,eps=eps)
}
