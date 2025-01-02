imputeByExpTilt <- function(y, ...)
  UseMethod("imputeByExpTilt")

imputeByExpTilt.default <- function(y, dpc.slope=0.7, prior.logfc=NULL, by="both", ...)
# Impute NAs using dpc and exponential tilting
# Created 21 July 2023. Last modified 23 May 2024.
{
# Check input y
  y <- as.matrix(y)

# Missing values
  imis <- which(is.na(y))

# If fewer than 2 observations overall, just return zeros
  global.nobs <- length(y) - length(imis)
  if(global.nobs < 1.5) {
    y[imis] <- 0
    warning("Less than 2 observations, exponential tilting not possible, replacing NAs with zeros")
    return(y)
  }

  by <- match.arg(by, c("rows","columns","both"))
  NRows <- nrow(y)
  NSamples <- ncol(y)
  if(identical(NRows,1L) && identical(by,"both")) {
    by <- "rows"
    if(identical(NSamples,1L)) return(y)
  }

# Impute by rows
  if(identical(by,"rows")) {
    br <- expTiltByRows(y, dpc.slope=dpc.slope, sigma.obs=prior.logfc)
    yt <- t(y)
    yt[is.na(yt)] <- rep.int(br$mu.missing,br$n.missing)
    y <- t(yt)
  }

# Impute by columns
  if(identical(by,"columns")) {
    bc <- expTiltByColumns(y, dpc.slope=dpc.slope)
    y[imis] <- rep.int(bc$mu.missing,bc$n.missing)
  }

# Impute by rows and columns
  if(identical(by,"both")) {
    br <- expTiltByRows(y, dpc.slope=dpc.slope, sigma.obs=prior.logfc)
    bc <- expTiltByColumns(y, dpc.slope=dpc.slope)

#   Use global exponential tilting for NA row or column values
    if(anyNA(br$mu.missing) || anyNA(bc$mu.missing)) {
      global.nobs <- length(y) - length(imis)
      global.degfree <- global.nobs - 1L
      global.mu.obs <- mean(y,na.rm=TRUE)
      global.sigma2.obs <- var(as.vector(y),na.rm=TRUE)
      global.mu.missing <- global.mu.obs - dpc.slope * global.sigma2.obs
      global.prediction.variance <- (global.nobs+1)/global.nobs * global.sigma2.obs + 2/global.degfree * (dpc.slope * global.sigma2.obs)^2
      i <- is.na(br$mu.missing)
      br$mu.missing[i] <- global.mu.missing
      br$prediction.variance[i] <- global.prediction.variance
      i <- is.na(bc$mu.missing)
      bc$mu.missing[i] <- global.mu.missing
      bc$prediction.variance[i] <- global.prediction.variance
    }
    
#   Weighted average of the row and column mu_mis
    i <- .row(c(NRows,NSamples))[is.na(y)]
    j <- .col(c(NRows,NSamples))[is.na(y)]
    mumisrw <- br$mu.missing / br$prediction.variance
    mumiscw <- bc$mu.missing / bc$prediction.variance
    mumis <- (mumisrw[i] + mumiscw[j]) / (1/br$prediction.variance[i] + 1/bc$prediction.variance[j])
    y[imis] <- mumis
  }

  y
}

imputeByExpTilt.EList <- function(y, ...)
# Impute NAs using dpc and exponential tilting for EList object
# Created 11 Sep 2023. Last modified 12 Sep 2023.
{
  y$E <- imputeByExpTilt(y$E, ...)
  y
}

imputeByExpTilt.EListRaw <- function(y, ...)
# Impute NAs using dpc and exponential tilting for EListRaw object
# Created 11 Sep 2023. Last modified 19 Sep 2023.
{
# Convert to EList object
  y$E[y$E < 1e-8] <- NA
  y$E <- log2(y$E)
  y <- new("EList",unclass(y))

  y$E <- imputeByExpTilt(y$E, ...)
  y
}

expTiltByRows <- function(y, dpc.slope=0.7, sigma.obs=NULL)
# Exponential tilting applied to rows
# Created 30 Oct 2023. Last modified 23 May 2024.
{
# Check input
  NRows <- nrow(y)
  NSamples <- ncol(y)

# Missing values
  nmis <- rowSums(is.na(y))
  nobs <- NSamples - nmis

# Row means
  m <- rowMeans(y,na.rm=TRUE)

# Standard deviation between samples
  if(is.null(sigma.obs)) {
    if(max(nobs) < 1.5) {
      sigma2.obs <- NA
    } else {
      degfree <- nobs-1L
      s2 <- rowSums((y-m)^2,na.rm=TRUE)
      i <- nobs > 1.5
      s2 <- s2[i]/degfree[i]
      sigma2.obs <- quantile(s2,0.9,na.rm=TRUE)
    }
  } else {
    sigma2.obs <- sigma.obs^2
  }
  degfree <- max(nobs)-1L

# Rowwise values for mu_mis
  mumis <- m - dpc.slope * sigma2.obs

# Prediction variance
  predvar <- (nobs+1)/nobs * sigma2.obs + 2/degfree * (dpc.slope * sigma2.obs)^2

  list(mu.missing=mumis,prediction.variance=predvar,n.missing=nmis)
}

expTiltByColumns <- function(y, dpc.slope=0.7)
# Exponential tilting applied to columns
# Created 30 Oct 2023. Last modified 23 May 2024.
{
# Check input
  NRows <- nrow(y)
  NSamples <- ncol(y)

# Missing values
  nmis <- colSums(is.na(y))
  nobs <- NRows - nmis

# Column means
  m <- colMeans(y,na.rm=TRUE)

# Standard deviation between rows
  degfree <- nobs-1L
  sigma2.obs <- rowSums((t(y)-m)^2,na.rm=TRUE)/degfree

# Columnwise values for mu_mis
  mumis <- m - dpc.slope * sigma2.obs

# Prediction variance
  predvar <- (nobs+1)/nobs * sigma2.obs + 2/degfree * (dpc.slope * sigma2.obs)^2

  list(mu.missing=mumis,prediction.variance=predvar,n.missing=nmis)
}
