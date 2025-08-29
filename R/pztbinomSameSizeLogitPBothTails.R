pztbinomSameSizeLogitPBothTails <- function(q, size, logit.prob)
# Tail probablities for zero-truncated binomial distribution
# assuming same `size` for all values and success probability on logit scale.
# Computes right and left tail p-values.
# Created 24 Jun 2025. Last modified 26 Jun 2025.
{
# Check input
  n <- length(q)
  if(length(size) > 1L) stop("Need same `size` for all values")

# Matrix of probabilities
  q1 <- 1:size
  lch <- lchoose(size,q1)
  logp <- plogis(logit.prob,lower.tail=TRUE,log.p=TRUE)
  log1p <- plogis(logit.prob,lower.tail=FALSE,log.p=TRUE)
  d <- matrix(lch,size,n) + matrix(q1,size,n)*matrix(logp,size,n,byrow=TRUE)+matrix(size-q1,size,n)*matrix(log1p,size,n,byrow=TRUE)
  d <- exp(d)

# Zero-truncated total probability
  total.prob <- colSums(d)

# Right and left p-values, allowing for fractional q values
  left.p.value <- right.p.value <- q
  left.p.value[] <- right.p.value[] <- 0
  qfloor <- pmin(floor(q),size)
  qceiling <- pmax(ceiling(q),1)
  ileft <- iright <- 1:n
  if(min(qfloor) < 1) ileft <- ileft[qfloor >= 1]
  if(max(qceiling) < 1) iright <- iright[qceiling <= size]
  for (j in ileft) left.p.value[j] <- sum(d[1:qfloor[j],j])
  for (j in iright) right.p.value[j] <- sum(d[qceiling[j]:size,j])
  left.p.value <- left.p.value / total.prob
  right.p.value <- right.p.value / total.prob

# Output
  list(left.p.value=left.p.value,right.p.value=right.p.value,total.prob=total.prob)
}
