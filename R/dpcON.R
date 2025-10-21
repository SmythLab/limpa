dpcON <- function(y, dpc.start=NULL, dpc.slope.start=0.7, robust=FALSE, verbose=FALSE)
# Estimate detection probability curve (DPC) assuming ON model.
# Created 22 Jun 2025. Last modified 12 Oct 2025.
{
# Check y
  y <- as.matrix(y)
  nsamples <- ncol(y)
  n.detected <- nsamples - rowSums(is.na(y))

# Remove entirely NA rows
  if(min(n.detected) < 0.5) {
    if(verbose) cat(sum(n.detected < 0.5), "peptides are completely missing in all samples.\n")
    y <- y[n.detected > 0.5,,drop=FALSE]
    n.detected <- n.detected[n.detected > 0.5]
  }
  npeptides <- nrow(y)

# Estimate hyperparamters
  hp <- .dpcHyperparam(y)
  mu.obs <- hp$mu_obs.post
  s2.obs <- hp$s2_obs.post
  mu.obs.mean <- mean(mu.obs)
  mu.obs.cor <- mu.obs - mu.obs.mean

# Negative log-likelihyood
  negLL <- function(par) {
    b0m <- par[1]
    b1 <- par[2]
    eta <- b0m + b1*mu.obs.cor - 0.5*b1*b1*s2.obs
    -sum(dztbinom(n.detected,size=nsamples,prob=eta,log=TRUE,logit.p=TRUE))
  }

# Starting values
  if(is.null(dpc.start)) {
#   Empirical logist transformation
    b1 <- dpc.slope.start
    LogitProb <- log( (n.detected + 0.5) / (nsamples - n.detected + 0.5) )
    w <- 1/( 1/(n.detected+0.5) + 1/(nsamples-n.detected+0.5) )
    X <- cbind(1,mu.obs-b1*s2.obs)
    fit <- lm.wfit(X,LogitProb,w)
    dpc.start <- as.vector(fit$coef)
  }
  b0 <- dpc.start[1]
  b1 <- dpc.start[2]
  if(verbose) message(" Start dpc = ",format(c(b0,b1)))

# Maximize likelihood by BFGS
  b0m <- b0 + b1*mu.obs.mean
  par <- c(b0m,b1)
  out <- optim(par,fn=negLL,method="BFGS")
  if(verbose) {
    b0m <- out$par[1]
    b1 <- out$par[2]
    b0 <- b0m - b1*mu.obs.mean
    message("   MLE dpc = ",format(c(b0,b1)))
  }

# Robust
  if(robust) {
    b0m <- out$par[1]
    b1 <- out$par[2]
    eta <- b0m + b1*mu.obs.cor - 0.5*b1*b1*s2.obs
    PZT <- pztbinomSameSizeLogitPBothTails(n.detected, nsamples, logit.prob=eta)
    p <- 2*pmin(PZT$left.p.value,PZT$right.p.value)
    FDR <- p.adjust(p,method="BH")
#   if(verbose) cat("Downweighting",sum(FDR < 0.05),"outliers\n")
#   FDR.left <- p.adjust(PZT$left.p.value)
#   if(min(FDR.left) < 0.05) {
#     i <- which(FDR.left < 0.05)
#     df.total <- (hp$df.prior + (n.detected-1))[i]
#     s2.obs[i] <- s2.obs[i] * df.total / qchisq(0.05,df=df.total)
#     if(verbose) cat("Floating",length(i),"observed variances\n")
#   }
    negLLWt <- function(par) {
      b0m <- par[1]
      b1 <- par[2]
      eta <- b0m + b1*mu.obs.cor - 0.5*b1*b1*s2.obs
      -sum(FDR * dztbinom(n.detected,size=nsamples,prob=eta,log=TRUE,logit.p=TRUE))
    }
    out <- optim(par,fn=negLLWt,method="BFGS")
    if(verbose) {
      b0m <- out$par[1]
      b1 <- out$par[2]
      b0 <- b0m - b1*mu.obs.mean
      message("Robust dpc = ",format(c(b0,b1)))
    }
  } else {
    FDR <- NULL
  }

# Output
  b0m <- out$par[1]
  b1 <- out$par[2]
  b0 <- b0m - b1*mu.obs.mean
  mu.mis <- mu.obs - b1 * s2.obs
  list(dpc = c(beta0=b0, beta1=b1),
       dpc.start = dpc.start, 
       n.detected = n.detected,
       nsamples = nsamples,
       mu.prior = hp$mu0,
       n.prior = hp$n0,
       df.prior = hp$df.prior,
       s2.prior = hp$s2.prior,
       mu.obs = mu.obs,
       mu.mis = mu.mis, 
       s2.obs = s2.obs,
       neg.loglik = out$value,
       outlier.FDR = FDR,
       model = "ON")
}
