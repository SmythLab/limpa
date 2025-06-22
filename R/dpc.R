dpc <- function(y, maxit = 100, eps = 1e-4, b1.upper = 1)
# Estimate detection probability curve (DPC) assuming ON model.
# Mengbo Li and Gordon Smyth
# Created 16 May 2022 as part of proDP package.
# Migrated to limpa package 10 Sept 2024. Last modified 31 Dec 2024.
{

  y <- as.matrix(y)
  narrays <- ncol(y)
  n.detected <- rowSums(!is.na(y))

  # Check input
  if (identical(min(n.detected),0)) {
    message(sum(n.detected == 0), " peptides are completely missing in all samples.")
    y <- y[n.detected > 0.5, ]
    n.detected <- n.detected[n.detected > 0.5]
  }

  dp <- n.detected / narrays
  wt <- rep_len(narrays, nrow(y))

  # Get hyperparamters
  hp <- .dpcHyperparam(y)
  mu_obs <- hp$mu_obs.post
  s2_obs <- hp$s2_obs.post

  # Get start values for betas assuming standard binomial
  glmfit0 <- glm(dp ~ mu_obs,
                 weights = wt,
                 family = binomial)
  betaStart <- coef(glmfit0)
  names(betaStart) <- c("beta0", "beta1")

  # Get an initial logit-linear fit under zero-truncated binomial
  fit0 <- .logitZTBinom(
    dp = dp,
    X = matrix(mu_obs, ncol = 1),
    wt = wt,
    beta0 = betaStart,
    b0.upper = 0,
    b1.upper = b1.upper
  )
  betaStart <- fit0$params
  mu_mis <- mu_obs - betaStart[2]*s2_obs

  betas <- betaStart
  betas.hist <- matrix(betas, nrow = 1)
  negLL <- .dpc.negLL(betas, dp, wt, mu_obs, mu_mis)
  negLL.hist <- negLL
  for (i in seq_len(maxit)) {
#### Can this call to optim() be replaced by fitZTLogit? Are the bounds required?
    ztbinomFit <- stats::optim(betas,
                               .dpc.negLL,
                               dp = dp,
                               wt = wt,
                               mu_obs = mu_obs,
                               mu_mis = mu_mis,
                               method = "L-BFGS-B",
                               lower = c(-Inf, 0),
                               upper = c(0, b1.upper))
    newBetas <- ztbinomFit$par
    mu_mis <- mu_obs - newBetas[2]*s2_obs
    newNegLL <- .dpc.negLL(newBetas, dp, wt, mu_obs, mu_mis)
    if (negLL - newNegLL < eps) break
    betas.hist <- rbind(betas.hist, newBetas)
    negLL.hist <- c(negLL.hist, newNegLL)
    betas <- newBetas
    negLL <- newNegLL
  }

  info <- cbind(betas.hist, negLL.hist)
  colnames(info) <- c("beta0", "beta1", "neg.ZBLL")
  rownames(info) <- paste0("i=", 0:(nrow(info)-1))

  list(dpc = betas,
       history = info,
       dpc.start = betaStart,
       n.detected = n.detected,
       nsamples = narrays,
       mu.prior = hp$mu0,
       n.prior = hp$n0,
       df.prior = hp$df.prior,
       s2.prior = hp$s2.prior,
       mu.obs = mu_obs,
       s2.obs = s2_obs,
       mu.mis = mu_mis)
}




.dpcHyperparam <- function(y)
# Obtain hyperparameters and empirical Bayes moderated mu_obs and s2_obs values for fitting the DPC.
# Migrated from the protDP package on 10 Sept 2024. Last modified 31 Dec 2024.
{
  # Observed sample means and variances
  nmis <- rowSums(is.na(y))
  nobs <- ncol(y) - nmis
  mu <- rowMeans(y, na.rm = TRUE)
  df.residual <- nobs - 1L
  s2 <- rowSums((y - mu)^2, na.rm = TRUE) / df.residual

  # Estimate variance hyperparameters
  sv <- squeezeVar(s2, df = df.residual)
  df.prior <- sv$df.prior
  df.total <- df.prior + df.residual
  s2.prior <- sv$var.prior
  s2.post <- sv$var.post

  # Estimate mu0
  mu0 <- mean(mu)

  # Estimate n0
  n <- nobs
  modt <- (mu - mu0) / sqrt(s2.post/n)
  n0 <- .n0EstimateFromModT(tstat = modt,
                            df = df.total,
                            n = n,
                            maxit = 20L,
                            eps = 1e-5,
                            trace = FALSE)

  # posterior mean and variance
  mu_obs.post <- (n*mu + n0*mu0) / (n + n0)
  if (is.infinite(df.prior)) {
    s2_obs.post <- sv$var.post
  } else {
    s2_obs.post <- (n*n0*(mu-mu0)^2/(n+n0) + (n-1)*s2 +
                      df.prior*s2.prior) /(n+df.prior)
    s2_obs.post[is.na(s2_obs.post)] <- s2.prior
  }

  list(mu0 = mu0,
       n0 = n0,
       df.prior = df.prior,
       s2.prior = s2.prior,
       mu_obs.post = mu_obs.post,
       s2_obs.post = s2_obs.post)
}



.n0EstimateFromModT <- function(tstat,
                               df,
                               n,
                               maxit = 20L,
                               eps = 1e-5,
                               trace = FALSE)
# Moment estimation of the effective sample size hyperparameter n0.
# tstat is assumed to be sqrt((n+n0)/n0) * T where T is t-distributed on df degrees of freedom.
# Created 4 Aug 2022 as part of protDP package.
# Migrated to limpa 10 Sept 2024. Last modified 31 Dec 2024.
{
  tstat <- as.numeric(tstat)
  ExpectedLogF <- mean(logmdigamma(df/2) - logmdigamma(1/2))
  RHS <- 2*mean(log(abs(tstat))) - ExpectedLogF

  # Monotonically convergent Newton iteration for v0 = 1/n0
  # Solve mean(log(1+n*v0)) = RHS
  if(RHS <= 0) return(Inf)
  v0 <- 0
  iter <- 0L
  repeat {
    iter <- iter+1L
    Diff <- RHS - mean(log(1+n*v0))
    Deriv <- mean( n/(1+n*v0) )
    Step <- Diff/Deriv
    v0 <- v0 + Step
    if(trace) message("iter=",iter, " v0=",v0," Step=",Step)
    if(Step < eps) break
    if(iter >= maxit) {
      warning("Iteration limit reached")
      break
    }
  }
  1/v0
}


.logitZTBinom <- function(dp, X, wt, beta0, b0.upper = 0, b1.upper = Inf)
# Fit an empirical logit spline assuming zero-truncated binomial distribution
#### This function and then next can probably be replaced by fitZTLogit(), which should be faster but doesn't place a bound on b1.
{
  df <- length(beta0) - 1
  params <- beta0
  params.hist <- matrix(params, nrow = 1)
  negLL <- .logitZTBinom.negLL(params, dp, wt, X)
  negLL.hist <- negLL

  # The order of parameters goes alpha, b0, b1, and etc.
  lower.bounds <- c(-Inf, rep(0, df))
  if (df > 0) {
    upper.bounds <- c(b0.upper, b1.upper, rep(Inf, df-1))
  } else {
    upper.bounds <- c(Inf)
  }
  ztbinomFit <- stats::optim(params,
                             .logitZTBinom.negLL,
                             dp = dp, wt = wt, X = X,
                             method = "L-BFGS-B",
                             lower = lower.bounds,
                             upper = upper.bounds)
  newParams <- ztbinomFit$par
  newnegLL <- .logitZTBinom.negLL(newParams, dp, wt, X)
  params.hist <- rbind(params.hist, newParams)
  negLL.hist <- c(negLL.hist, newnegLL)

  # Clean up results
  info <- cbind(params.hist, negLL.hist)
  colnames(info) <- c(names(beta0), "neg.LL")
  rownames(info) <- paste("iter", 0:(nrow(info)-1), sep = " ")

  list(params = newParams, info = info)

}




.logitZTBinom.negLL <- function(params, dp, wt, X)
# Negative log-likelihood under zero-truncated binomial distribution to fit an empirical logit spline
# This is the objective function for logitZTBinom().
# Created 11 Sep 2024. Last modified 21 Jun 2025.
{
  df <- length(params) - 1
  if (df > 0) X <- cbind(1, X)
  eta <- colSums(t(X) * params)
  -sum(dztbinom(x = dp*wt, size = wt, prob = eta, log = TRUE, logit.p=TRUE))
}




.dpc.negLL <- function(params, dp, wt, mu_obs, mu_mis)
# Negative log-likelihood under zero-truncated binomial distribution to fit DPC
# This is the objective function for dpc().
# Created 11 Sep 2024. Last modified 21 Jun 2025.
{
  b0 <- params[1]
  b1 <- params[2]
  eta <- b0 + 0.5*b1*(mu_obs + mu_mis)
  -sum(dztbinom(x = dp*wt, size = wt, prob = eta, log = TRUE, logit.p=TRUE))
}
