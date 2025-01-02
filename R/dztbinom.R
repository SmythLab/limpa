dztbinom <- function(x, size, prob, log = FALSE)
# Density for zero-truncated binomial distribution
# Created 12 Aug 2023.
{
# Log-probabilities for zero-truncated binomial
  d <- dbinom(x=x,size=size,prob=prob,log=TRUE) - pbinom(0.5,size=size,prob=prob,log.p=TRUE,lower.tail=FALSE)

# Make sure that zero values get zero density
  if(min(x,na.rm=TRUE) < 1) d[x==0] <- -Inf

# Return result
  if(log) d else exp(d)
}

pztbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE)
# Cumulative probablities for zero-truncated binomial distribution
# Created 12 Aug 2023. Last modified 23 Aug 2023.
{
  if(lower.tail) {
    p <- pbinom(q=q,size=size,prob=prob,lower.tail=TRUE,log.p=FALSE)
    p <- p - pbinom(0.5,size=size,prob=prob,lower.tail=TRUE,log.p=FALSE)
    p <- p / pbinom(0.5,size=size,prob=prob,lower.tail=FALSE,log.p=FALSE)
    if(log.p) p <- log(p)
  } else {
    p <- pbinom(q=q,size=size,prob=prob,lower.tail=FALSE,log.p=TRUE)
    p <- p - pbinom(0.5,size=size,prob=prob,lower.tail=FALSE,log.p=TRUE)
    if(!log.p) p <- exp(p)
  }

# Return result
  p
}
