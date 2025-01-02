fitZTLogit <- function(n.successes, n.trials, X=NULL, capped=FALSE, beta.start=NULL, alpha.start=0.95)
# Fit capped logit regression to zero-truncated binomial data
# Created 22 Aug 2023. Last modified 21 Dec 2024.
{
# Check X matrix
  if(is.null(X)) X <- matrix(1,length(n.successes),1)
  nbeta <- ncol(X)
  ntheta <- nbeta + capped

# Minus twice log-likelihood function
  if(capped) {
    M2LL <- function(theta) {
      logitalpha <- theta[1L]
      beta <- theta[-1L]
      eta <- X %*% beta
      p <- plogis(logitalpha) * plogis(eta)
      -2*sum(dztbinom(n.successes,n.trials,prob=p,log=TRUE))
    }
  } else {
    M2LL <- function(theta) {
      eta <- X %*% theta
      p <- plogis(eta)
      -2*sum(dztbinom(n.successes,n.trials,prob=p,log=TRUE))
    }
  }

# Derivative of minus twice log-likelihood
  if(capped) {
    derivM2LL <- function(theta) {
      logitalpha <- theta[1]
      beta <- theta[-1]
      eta <- X %*% beta
      p <- plogis(logitalpha) * plogis(eta)
      deriv.p <- (n.successes - n.trials*p) / (p*(1-p))
      deriv.p <- deriv.p - dbeta(p,1,n.trials) / pbeta(p,1,n.trials)
      deriv.logitalpha <- sum(deriv.p * dlogis(logitalpha) * plogis(eta))
      deriv.beta <- t(X) %*% (deriv.p * plogis(logitalpha) * dlogis(eta))
      -2*c(deriv.logitalpha,deriv.beta)
    }
  } else {
    derivM2LL <- function(theta) {
      eta <- X %*% theta
      p <- plogis(eta)
      deriv.p <- (n.successes - n.trials*p) / (p*(1-p))
      deriv.p <- deriv.p - dbeta(p,1,n.trials) / pbeta(p,1,n.trials)
      -2 * (t(X) %*% (deriv.p * dlogis(eta)))
    }
  }

# Starting values
  if(is.null(beta.start)) {
    eta <- log((n.successes+0.5)/(n.trials-n.successes+0.5))
    beta.start <- lm.fit(X,eta)$coef
  }
  if(capped) {
    logitalpha.start <- qlogis(alpha.start)
    theta <- c(logitalpha.start,beta.start)
  } else {
    theta <- beta.start
  }

# Optimization
  bfgs <- optim(theta,fn=M2LL,gr=derivM2LL,method="BFGS")
  
# Output
  out <- list()
  if(capped) {
    out$beta <- bfgs$par[-1]
    out$alpha <- as.vector(plogis(bfgs$par[1]))
    out$p <- out$alpha * plogis(X %*% out$beta)
  } else{
    out$beta <- bfgs$par
    out$alpha <- 1
    out$p <- plogis(X %*% out$beta)
  }
  if(is.null(names(out$beta))) names(out$beta) <- paste0("b",0:(nbeta-1L))
  out$deviance <- bfgs$value
  out$calls <- bfgs$counts
  out
}
