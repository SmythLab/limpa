simCompleteDataCN <- function(n, mean.comp=6, sd.comp=1, dpc=c(-4,0.7))
# Simulate a vector of complete data from complete normal (CN) model.
# Created 25 Dec 2024.
{
  mu.comp <- rep_len(mean.comp, length.out=n)
  sigma.comp <- rep_len(sd.comp, length.out=n)
  y.complete <- rnorm(n, mean=mu.comp, sd=sigma.comp)
  prob.mis <- plogis(dpc[1] + dpc[2] * y.complete, lower.tail=FALSE)
  is.mis <- as.logical(rbinom(n, prob=prob.mis, size=rep_len(1L,n)))
  list(y.complete=y.complete, is.missing=is.mis, prob.missing=prob.mis)
}

simCompleteDataON <- function(n, mean.obs=6, sd.obs=1, dpc=c(-4,0.7))
# Simulate a vector of complete data from observed normal (ON) model.
# Created 25 Dec 2024.
{
  mu.obs <- rep_len(mean.obs, length.out=n)
  sigma.obs <- rep_len(sd.obs, length.out=n)
  mu.mis <- mu.obs - dpc[2] * sigma.obs^2
  sigma.mis <- sigma.obs
  y.obs <- rnorm(n, mean=mu.obs, sd=sigma.obs)
  y.mis <- rnorm(n, mean=mu.mis, sd=sigma.mis)
  prob.mis <- plogis(dpc[1] + dpc[2]*mu.obs - 0.5*(dpc[2]*sigma.obs)^2, lower.tail=FALSE)
  is.mis <- as.logical(rbinom(n, prob=prob.mis, size=rep_len(1L,n)))
  y.complete <- y.obs
  y.complete[is.mis] <- y.mis[is.mis]
  list(y.complete=y.complete, is.missing=is.mis, prob.missing=prob.mis)
}
