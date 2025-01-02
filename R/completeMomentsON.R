observedMomentsCN <- function(mean.comp=6, sd.comp=1, dpc=c(-4,0.7))
# Simulate a vector of complete data from complete normal (CN) model.
# Created 25 Dec 2024.
{
  n <- max(length(mean.comp),length(sd.comp))
  mu.comp <- rep_len(mean.comp, length.out=n)
  sigma.comp <- rep_len(sd.comp, length.out=n)
  gq <- gauss.quad.prob(16, dist="normal")
  y <- matrix(gq$nodes,n,16,byrow=TRUE)
  y <- mu.comp + sigma.comp * y
  prob.obs <- drop(plogis(dpc[1] + dpc[2]*y) %*% gq$weights)
  mu.obs <- drop((plogis(dpc[1] + dpc[2]*y) * y) %*% gq$weights) / prob.obs
  sigma2.obs <- drop((plogis(dpc[1] + dpc[2]*y) * (y-mu.obs)^2) %*% gq$weights) / prob.obs
  list(mean.obs=mu.obs, sd.obs=sqrt(sigma2.obs), prob.obs=prob.obs)
}

completeMomentsON <- function(mean.obs=6, sd.obs=1, dpc=c(-4,0.7))
# Mean and variance of complete data given observed normal (ON) model.
# Created 25 Dec 2024.
{
  n <- max(length(mean.obs),length(sd.obs))
  mu.obs <- rep_len(mean.obs, length.out=n)
  sigma.obs <- rep_len(sd.obs, length.out=n)
  mu.mis <- mu.obs - dpc[2] * sigma.obs^2
  prob.obs <- plogis(dpc[1] + dpc[2]*mu.obs - 0.5*(dpc[2]*sigma.obs)^2)
  prob.mis <- 1 - prob.obs
  mu.comp <- mu.obs - prob.mis * dpc[2] * sigma.obs^2
  sigma2.comp <- sigma.obs^2 + prob.obs * (mu.obs - mu.comp)^2 + prob.mis * (mu.mis - mu.comp)^2
  list(mean.comp=mu.comp, sd.comp=sqrt(sigma2.comp), prob.obs=prob.obs)
}
