peptides2Proteins <- function(y, protein.id, sigma=0.5, dpc=c(-4,0.7), prior.mean=6, prior.sd=10, prior.logFC=2, standard.errors=FALSE, newton.polish=FALSE, verbose=FALSE, chunk=1000L)
# Summarize peptide to protein log-expression for many proteins.
# Created 10 July 2023. Last modified 29 Dec 2024.
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

# Output matrices
  z <- matrix(0,nproteins,nsamples)
  rownames(z) <- protein.id.unique
  colnames(z) <- colnames(y)
  nobs <- matrix(0L,nproteins,nsamples)
  dimnames(nobs) <- dimnames(z)
  if(standard.errors) stderr <- z

# Summarize peptide to protein expression
  last.row <- 0
  for (i.protein in seq_len(nproteins)) {
    i <- last.row + seq_len(npeptides[i.protein])
    yi <- y[i,,drop=FALSE]
    out <- peptides2ProteinBFGS(yi, sigma=sigma[i.protein], dpc=dpc,
      prior.mean=prior.mean, prior.sd=prior.sd, prior.logFC=prior.logFC,
      standard.errors=standard.errors, newton.polish=newton.polish)
    z[i.protein,] <- out$protein.expression
    nobs[i.protein,] <- colSums(!is.na(yi))
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
