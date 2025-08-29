dpcCN <- function(y, dpc.slope.start=0.7, dpc.start=NULL, iterations=2L, subset=1000L, verbose=TRUE)
# MLE for DPC curve assuming complete normal model.
# Created 14 Dec 2024. Last modified 23 Jun 2025.
{
# Check y
  y <- as.matrix(y)
  npeptides <- nrow(y)
  nsamples <- ncol(y)
  if(npeptides < 3) stop("Too few rows of data")
  if(nsamples < 2) stop("Too few samples")

# Subset large datasets
  if(npeptides > subset) {
    set.seed(20250620)
    invisible(runif(100))
    i <- sample.int(npeptides, subset)
    y <- y[i,]
    npeptides <- subset
  }

# DPC
  if(is.null(dpc.start)) {
    beta1 <- dpc.slope.start
    beta0 <- estimateDPCIntercept(y,dpc.slope=beta1)
  } else {
    beta0 <- dpc.start[1]
    beta1 <- dpc.start[2]
  }

# Remove rows that are entirely missing
#  y <- y[rowSums(is.na(y)) < nsamples,,drop=FALSE]
#  npeptides <- nrow(y)

# Gauss quadrature
  gq <- gauss.quad.prob(16,dist="normal")

# Minus twice log-likelihood for beta0 and beta1
  minusLogLik <- function(par,y,mu,sigma,gq) {
    beta0 <- par[1]
    beta1 <- exp(par[2]/10)
    Q <- 0
    for (i in seq_len(npeptides)) {
      yi <- y[i,]
      ismis <- is.na(yi)
      isobs <- !ismis
      nmis <- sum(ismis)
      eta <- beta0+beta1*yi[isobs]
      Q <- Q + sum(plogis(eta,log.p=TRUE))
      if(nmis) {
        eta <- beta0+beta1*(mu[i]+sigma[i]*gq$nodes)
        ProbMis <- sum(gq$weights * plogis(eta,lower.tail=FALSE))
        Q <- Q + nmis*log(ProbMis)
      }
    }
    -2*Q
  }

# Deriv for beta0 and beta1
  dMinusLogLik <- function(par,y,mu,sigma,gq) {
    beta0 <- par[1]
    beta1 <- exp(par[2]/10)
    DQ.beta0 <- DQ.beta1 <- 0
    for (i in seq_len(npeptides)) {
      yi <- y[i,]
      ismis <- is.na(yi)
      isobs <- !ismis
      nmis <- sum(ismis)
      eta <- beta0+beta1*yi[isobs]
      ProbObs <- plogis(eta)
      DProbObs.beta0 <- dlogis(eta)
      DProbObs.beta1 <- yi[isobs]*dlogis(eta)
      DQ.beta0 <- DQ.beta0 + sum(DProbObs.beta0/ProbObs)
      DQ.beta1 <- DQ.beta1 + sum(DProbObs.beta1/ProbObs)
      if(nmis) {
        z <- mu[i]+sigma[i]*gq$nodes
        eta <- beta0+beta1*z
        ProbMis <- sum(gq$weights * plogis(eta,lower.tail=FALSE))
        DProbMis.beta0 <- -sum(gq$weights * dlogis(eta))
        DProbMis.beta1 <- -sum(gq$weights * z * dlogis(eta))
        DQ.beta0 <- DQ.beta0 + nmis*DProbMis.beta0/ProbMis
        DQ.beta1 <- DQ.beta1 + nmis*DProbMis.beta1/ProbMis
      }
    }
    -2*c(DQ.beta0,DQ.beta1*beta1/10)
  }

# Starting estimates for mu and sigma for each peptide
  y.impute <- imputeByExpTilt(y,dpc.slope=beta1)
  mu <- rowMeans(y.impute)
  sigma <- sqrt(rowSums((y.impute-mu)^2)/(nsamples-1))

# Prior distribution
  mu.mean <- mean(mu)
  mu.sd <- sd(mu)
  logsigma <- log(sigma+0.001)
  logsigma.mean <- mean(logsigma)
  logsigma.sd <- sd(logsigma)

# Main iteration
  for (iter in seq_len(iterations)) {
    if(verbose) message("Iter ",iter,": ",appendLF=FALSE)

#   Estimate mu and sigma for each peptide
    for (i in seq_len(npeptides)) {
#     Maximize posterior
      out <- .estimateMuSigmaForOnePeptide(y[i,],dpc=c(beta0,beta1),mu.start=mu[i],sigma.start=sigma[i],mu.mean=mu.mean,mu.sd=mu.sd,logsigma.mean=logsigma.mean,logsigma.sd=logsigma.sd)
      mu[i] <- out$mu
      sigma[i] <- out$sigma
    }

#   BFGS for beta0 and beta1
    par <- c(beta0,log(beta1)*10)
    out <- optim(par,fn=minusLogLik,gr=dMinusLogLik,y=y,mu=mu,sigma=sigma,gq=gq,method="BFGS")
    beta0 <- out$par[1]
    beta1 <- exp(out$par[2]/10)
    if(verbose) message("dpc = ",format(c(beta0,beta1)))
  }

# Output
  out <- list()
  out$dpc <- c(beta0, beta1)
  names(out$dpc) <- c("beta0", "beta1")
  out$mu <- mu
  out$sigma <- sigma
  out$n.detected <- rowSums(!is.na(y))
  out$nsamples <- nsamples
  out
}

.estimateMuSigmaByRow <- function(y, dpc=c(-4,0.7), mu.start=NULL, sigma.start=NULL, mu.mean=6, mu.sd=5, logsigma.mean=0, logsigma.sd=3, gq=gauss.quad.prob(16,dist="normal"))
# Maximum posterior estimates for mu and sigma for each row of a matrix
# given DPC and assuming complete normal model.
# Created 22 Dec 2024. Last modified 31 Dec 2024.
{
  y <- as.matrix(y)
  npeptides <- nrow(y)
  if(is.null(mu.start)) mu.start <- rep_len(mu.mean,npeptides)
  if(is.null(sigma.start)) mu.start <- rep_len(exp(logsigma.mean),npeptides)
  for (i in seq_len(nrow(y))) {
    out <- .estimateMuSigmaForOnePeptide(y[i,],dpc=dpc,mu.start=mu[i],sigma.start=sigma[i],mu.mean=mu.mean,mu.sd=mu.sd,logsigma.mean=logsigma.mean,logsigma.sd=logsigma.sd,gq=gq)
    mu[i] <- out$mu
    sigma[i] <- out$sigma
  }
  list(mu=mu,sigma=sigma)
}

.estimateMuSigmaForOnePeptide <- function(y, dpc=c(-4,0.7), mu.start=NULL, sigma.start=NULL, mu.mean=6, mu.sd=5, logsigma.mean=0, logsigma.sd=3, gq=gauss.quad.prob(16,dist="normal"))
# Maximum posterior estimates for mu and sigma for one peptide
# given DPC and assuming complete normal model.

# Input:
# y: vector log-intensities with NAs
# dpc: detection probability curve intercept and slope
# mu.start: starting value
# sigma.start: starting value
# mu.mean: prior mean for mu
# mu.sd: prior sd for mu
# logsigma.mean: prior mean for logsigma
# logsigma.sdn: prior sd for logsigma

# Created 14 Dec 2024. Last modified 23 Dec 2024.
{
# Observed and missing values
  mis <- is.na(y)
  obs <- !mis
  mis <- which(mis)
  obs <- which(obs)
  nmis <- length(mis)
  nobs <- length(obs)

# DPC
  dpc.intercept <- dpc[1]
  dpc.slope <- dpc[2]

# Starting values
  if(is.null(mu.start)) mu.start <- mu.mean
  if(is.null(sigma.start)) sigma.start <- exp(logsigma.mean)

# Minus twice the log-posterior function
  yo <- y[obs]
  minusLogPosterior <- function(par) {
    mu <- par[1]
    logsigma <- par[2]
    sigma <- exp(logsigma)
#   Prior part
    Q <- ((mu-mu.mean)/mu.sd)^2 + ((logsigma-logsigma.mean)/logsigma.sd)^2
#   Observed part
#   Using REML (nobs-1) instead of ML nobs
    if(nobs) Q <- Q + 2*(nobs-1)*logsigma + sum((yo-mu)^2)/sigma^2
#   Missing part
    if(nmis) {
      eta <- dpc.intercept + dpc.slope * (mu + sigma * gq$nodes)
      G <- sum(gq$weights * plogis(eta,lower.tail=FALSE))
      Q <- Q - 2*nmis*log(G)
    }
#   Total
    Q
  }

# Derivatives function
  Deriv <- c(0,0)
  dMinusLogPosterior <- function(par) {
    mu <- par[1]
    logsigma <- par[2]
    sigma <- exp(logsigma)
#   Prior part
    DQ.mu <- 2*(mu-mu.mean)/mu.sd^2
    DQ.logsigma <- 2*(logsigma-logsigma.mean)/logsigma.sd^2
#   Observed part
    if(nobs) {
      DQ.mu <- DQ.mu - 2*sum(yo-mu)/sigma^2
      DQ.logsigma <- DQ.logsigma + 2*(nobs-1) - 2*sum((yo-mu)^2)/sigma^2
    }
#   Missing part
    if(nmis) {
      eta <- dpc.intercept + dpc.slope * (mu + sigma * gq$nodes)
      G <- sum(gq$weights * plogis(eta,lower.tail=FALSE))
      DQ.mu <- DQ.mu + 2*nmis * dpc.slope * sum(gq$weights*dlogis(eta)) / G
      DQ.logsigma <- DQ.logsigma + 2*nmis * sigma * dpc.slope * sum(gq$weights*gq$nodes*dlogis(eta)) / G
    }
#   Total
    Deriv[1] <- DQ.mu
    Deriv[2] <- DQ.logsigma
    Deriv
  }

# Maximize log-likelihood
  par <- c(mu.start,log(sigma.start))
  out <- optim(par,fn=minusLogPosterior,gr=dMinusLogPosterior,method="BFGS")

# Output
  out$mu <- out$par[1]
  out$sigma <- exp(out$par[2])
  out$deriv <- dMinusLogPosterior(out$par)
  out$model <- "CN"
  out
}
