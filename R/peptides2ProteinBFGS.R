peptides2ProteinBFGS <- function(y, sigma=0.5, weights=NULL, dpc=c(-4,0.7), prior.mean=6, prior.sd=10, prior.logFC=2, standard.errors=TRUE, newton.polish=TRUE, start=NULL)
# Summarize peptide to protein log-expression for one protein
# assuming a complete normal model, additive linear model,
# logit-linear detection probabilities and
# and a normal prior distribution for protein log-expression.
# Uses BFGS algorithm with derivatives and improved calculation.
# Computes standard errors from analytic Hessian matrix.
# Input:
# y: numeric matrix of log-intensities, rows are peptides, columns are samples
# sigma: complete data residual standard deviation for the additive linear model
# Created 5 June 2023. Last modified 25 Jun 2025.
{
# log-intensities
  npeptides <- nrow(y)
  nsamples <- ncol(y)

# Check weights
  if(is.null(weights)) {
    weights <- matrix(1,npeptides,nsamples)
  } else {
    if(!identical(dim(y),dim(weights))) stop("weights must be conformal with y")
    if(min(weights) < 1e-15) stop("weights must be positive")
  }

# Consolidate entirely NA rows
  AllNA <- which(rowSums(is.na(y))==nsamples)
  if(length(AllNA) > 1L) {
    y <- y[-AllNA[-1],,drop=FALSE]
    weights[AllNA[1],] <- colSums(weights[AllNA,])
    weights <- weights[-AllNA[-1],,drop=FALSE]
    if(!is.null(start)) start <- start[-(nsamples-1+AllNA[-1])]
    npeptides <- npeptides +1L - length(AllNA)
  }

# DPC parameters
  dpc.intercept <- dpc[1]
  dpc.slope <- dpc[2]

# Additive model for peptides and samples
  i1n <- seq_len(nsamples)
  Samples <- factor(rep.int(i1n,rep.int(npeptides,nsamples)))
  if(npeptides > 1L) {
    Peptides <- factor(rep.int(seq_len(npeptides),nsamples))
    contrasts(Peptides) <- contr.sum(npeptides)
    X <- model.matrix(~0+Samples+Peptides)
  } else {
    X <- model.matrix(~0+Samples)
  }
  N <- nrow(X)
  nbeta <- ncol(X)

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
  yo <- y[obs]
  sigmao <- sigma[obs]
  sigmam <- sigma[mis]
  probmis <- rep_len(0,length(mis))
  dprobmis <- probmis
  Xo <- X[obs,,drop=FALSE]
  Xm <- X[mis,,drop=FALSE]

# Gauss quadrature
  if(length(mis)) {
    gq <- gauss.quad.prob(16,dist="normal")
    z <- matrix(gq$nodes,length(mis),16,byrow=TRUE)
  }

# Minus twice log-posterior
  minusLogPosterior <- function(beta) {
    mu <- X %*% beta
#   Observed part
    if(length(obs)) {
      logposto <- sum(((yo-mu[obs])/sigmao)^2)
    } else {
      logposto <- 0
    }
#   Missing part
    if(length(mis)) {
      eta <- dpc.intercept + dpc.slope * (mu[mis] + sigmam * z)
      G <- plogis(eta,lower.tail=FALSE) %*% gq$weights
      logpostm <- -2*sum(log(G))
    } else {
      logpostm <- 0
    }
#   Prior part
    ProteinExpr <- beta[i1n]
    ProteinExpr.Mean <- mean(ProteinExpr)
    logpostp <- (ProteinExpr.Mean - prior.mean)^2 / prior.sd^2
    logpostp <- logpostp + sum( (ProteinExpr - ProteinExpr.Mean)^2 ) / prior.logFC^2
#   Return
    logposto + logpostm + logpostp
  }

# Derivative of minus twice log-posterior
  dMinusLogPosterior <- function(beta) {
    mu <- X %*% beta
    N <- nrow(X)
    DQ.Dmu <- rep_len(0,N)
#   Observed part
    DQ.Dmu[obs] <- -(yo-mu[obs])/sigmao^2
#   Missing part
    if(length(mis)) {
      mum <- mu[mis]
      eta <- dpc.intercept + dpc.slope * (mum + sigmam * z)
      G <- plogis(eta,lower.tail=FALSE) %*% gq$weights
      Gdot <- -dpc.slope * dlogis(eta) %*% gq$weights
      DQ.Dmu[mis] <- -Gdot/G
    }
    DQ <- crossprod(X,DQ.Dmu)
#   Prior part
    ProteinExpr <- beta[i1n]
    ProteinExpr.Mean <- mean(ProteinExpr)
    DQ.Prior <- (ProteinExpr.Mean-prior.mean)/nsamples/prior.sd^2 + (ProteinExpr-ProteinExpr.Mean)/prior.logFC^2
    DQ[i1n] <- DQ[i1n] + DQ.Prior
#   Return
    2*DQ
  }

# Minimize minus twice log-posterior
  out <- optim(beta,fn=minusLogPosterior,gr=dMinusLogPosterior,method="BFGS",hessian=FALSE)

# Output protein expression values
  ProteinExpression <- out$par[i1n]
  names(ProteinExpression) <- SampleNames
  out$protein.expression <- ProteinExpression

# Output derivatives and standard errors
  if(standard.errors) {
    ProteinExpr.Mean <- mean(ProteinExpression)
    beta <- out$par
    mu <- X %*% beta
    DQ.Dmu <- D2Q.Dmu <- rep_len(0,N)
    DQ.Dmu[obs] <- -(yo-mu[obs])/sigmao^2
    D2Q.Dmu[obs] <- 1/sigmao^2
    if(length(mis)) {
      eta <- dpc.intercept + dpc.slope * (mu[mis] + sigmam * z)
      G <- plogis(eta,lower.tail=FALSE) %*% gq$weights
      Gdot <- -dpc.slope * dlogis(eta) %*% gq$weights
      Gddot <- dpc.slope^2/4 * (tanh(eta/2) / cosh(eta/2)^2) %*% gq$weights
      DQ.Dmu[mis] <- -Gdot/G
      D2Q.Dmu[mis] <- (Gdot*Gdot - G*Gddot)/G/G
    }
    DQ <- crossprod(X,DQ.Dmu)
    DQ.Prior <- (ProteinExpr.Mean - prior.mean)/nsamples/prior.sd^2 + (ProteinExpression-ProteinExpr.Mean)/prior.logFC^2
    DQ[i1n] <- DQ[i1n] + DQ.Prior
    DQ <- 2*DQ
    out$deriv <- as.vector(DQ)
    D2Q <- crossprod(X,D2Q.Dmu*X)
    D2Q.Prior <- 1/nsamples^2/prior.sd^2 - 1/nsamples/prior.logFC^2 + diag(nsamples)/prior.logFC^2
    D2Q[i1n,i1n] <- D2Q[i1n,i1n] + D2Q.Prior
    D2Q <- 2*D2Q
    out$hessian <- D2Q

#   Invert D2Q with error catching
    tryCatch(Qinv <- chol2inv(chol(D2Q)), error=function(e) {
      warning("Hessian not invertible at first try.")
      eps <- max(diag(D2Q)) * 1e-8
      eps <- rep(c(0,eps),c(nsamples,npeptides-1))
      diag(D2Q) <- diag(D2Q) + eps
      Qinv <- chol2inv(chol(D2Q))
    })

#   Use D2Q and DQ to undertake one Newton iteration for almost no cost
    if(newton.polish) {
      beta <- beta - Qinv %*% DQ
      out$protein.expression <- beta[i1n]
      names(out$protein.expression) <- SampleNames
      out$standard.error <- sqrt(2*diag(Qinv)[i1n])
      out$par.bgfs <- out$par
      out$par <- beta
      out$value.bfgs <- out$value
      out$value <- minusLogPosterior(beta)
    } else {
      out$standard.error <- sqrt(2*diag(Qinv)[i1n])
    }
  }

  out
}
