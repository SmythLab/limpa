dztbinom <- function(x, size, prob, log=FALSE, logit.p=FALSE)
# Density for zero-truncated binomial distribution.
# If logit.p=TRUE, then prob is on the logit scale.
# Created 12 Aug 2023. Last modified 21 Jun 2025.
{
# Compute probabilities on the log-scale
  if(logit.p) {
#   Ensure all arguments have same length
    n <- max(length(x),length(size),length(prob))
    if(length(x) < n) x <- rep_len(x,n)
    size <- rep_len(size,n)
    prob <- rep_len(prob,n)
    pos <- prob > 0
    neg <- !pos
    pos <- which(pos)
    neg <- which(neg)
    d <- x
    if(length(neg)) {
      p <- plogis(prob[neg],lower.tail=TRUE)
      s <- size[neg]
      d[neg] <- dbinom(x=x[neg],size=s,prob=p,log=TRUE) - pbinom(0.5,size=s,prob=p,log.p=TRUE,lower.tail=FALSE)
    }
#   For large probabilities, convert to opposite tail
#   to avoid subtractive cancellation for probabilities near 1.
    if(length(pos)) {
      p <- plogis(prob[pos],lower.tail=FALSE)
      s <- size[pos]
      d[pos] <- dbinom(x=s-x[pos],size=s,prob=p,log=TRUE) - pbinom(s-0.5,size=s,prob=p,log.p=TRUE,lower.tail=TRUE)
    }
  } else {
    d <- dbinom(x=x,size=size,prob=prob,log=TRUE) - pbinom(0.5,size=size,prob=prob,log.p=TRUE,lower.tail=FALSE)
  }

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
