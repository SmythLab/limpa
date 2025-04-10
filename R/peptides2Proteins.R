peptides2Proteins <- function(y, protein.id, sigma=0.5, dpc=c(-4,0.7), prior.mean=6, prior.sd=10, prior.logFC=2, standard.errors=FALSE, newton.polish=FALSE, verbose=FALSE, chunk=1000L)
# Summarize peptide to protein log-expression for many proteins.
# Created 10 July 2023. Last modified 9 March 2025.
{
# Check y
  y <- as.matrix(y)
  NRows <- nrow(y)
  nsamples <- ncol(y)

# Check protein.id
  if(anyNA(protein.id)) stop("NAs not allowed in protein.id")
  if(!identical(NRows,length(protein.id))) stop("length of protein.id must agree with row dimension of y")

# Check that peptides are in protein order
  if( !identical(protein.id,sort(protein.id)) ) {
    o <- order(protein.id)
    protein.id <- protein.id[o]
    y <- y[o,,drop=FALSE]
  }

# Count peptides per protein
  d <- which(!duplicated(protein.id))
  nproteins <- length(d)
  protein.id.unique <- protein.id[d]
  d <- c(d,NRows+1)
  i <- seq_len(nproteins)
  npeptides <- d[i+1L]-d[i]

# Check sigma
  nsigma <- length(sigma)
  if(identical(nsigma,1L)) {
    sigma <- rep_len(sigma,nproteins)
  } else {
    if(!identical(nsigma,nproteins)) stop("Length of sigma must equal number of unique proteins")
  }

# Check chunk
  if(verbose) chunk <- as.integer(chunk)

# Output matrices
  z <- matrix(0,nproteins,nsamples)
  rownames(z) <- protein.id.unique
  colnames(z) <- colnames(y)
  nobs <- matrix(0L,nproteins,nsamples)
  dimnames(nobs) <- dimnames(z)
  if(standard.errors) stderr <- z

# Simple imputation to get starting values
  yimp <- imputeByExpTilt(y,dpc.slope=dpc[2],prior.logfc=prior.logFC)

# Summarize peptide to protein expression
  last.row <- 0
  for (i.protein in seq_len(nproteins)) {
    npeptidesi <- npeptides[i.protein]
    i <- last.row + seq_len(npeptidesi)
    yi <- y[i,,drop=FALSE]
    nobsi <- colSums(!is.na(yi))
#   Starting values
    yimpi <- yimp[i,,drop=FALSE]
    if(npeptidesi > 1L) {
      nbeta <- nsamples + npeptidesi - 1L
      beta <- rep_len(0,nbeta)
      beta[1:nsamples] <- colMeans(yimpi)
      b <- rowMeans(yimpi)
      b <- b-mean(b)
      beta[(nsamples+1):nbeta] <- b[-npeptidesi]
    } else {
      beta <- drop(yimpi)
    }
    IsImp <- which(nobsi==0)
    if(length(IsImp) > 1L) beta[IsImp] <- mean(beta[IsImp])
#   Maximize posterior
    out <- peptides2ProteinBFGS(yi, sigma=sigma[i.protein], dpc=dpc,
      prior.mean=prior.mean, prior.sd=prior.sd, prior.logFC=prior.logFC,
      standard.errors=standard.errors, newton.polish=newton.polish, start=beta)
#   Assemble output
    z[i.protein,] <- out$protein.expression
    nobs[i.protein,] <- nobsi
    if(standard.errors) stderr[i.protein,] <- out$standard.error
    last.row <- last.row + npeptides[i.protein]
    if(verbose) {
      if(identical(i.protein %% chunk,0L)) message("Proteins: ", i.protein, " Peptides: ", last.row)
    }
  }
  if(verbose) {
    message("Proteins: ", i.protein, " Peptides: ", last.row)
  }

# Collect output
  out <- list(E=z)
  out$genes <- data.frame(NPeptides=npeptides,PropObs=rowMeans(nobs)/npeptides)
  out$other$n.observations <- nobs
  if(standard.errors) out$other$standard.error <- stderr
  new("EList",out)
}
