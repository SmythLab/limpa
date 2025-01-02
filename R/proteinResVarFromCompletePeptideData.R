proteinResVarFromCompletePeptideData <- function(y, protein.id, reorder=FALSE)
# Get protein-wise residual variances from complete peptide data
# by fitting an additive model for each protein.
# Created 06 Nov 2023. Last modified 29 Dec 2024.
{
# Check y
  y <- as.matrix(y)
  if(anyNA(y)) stop("y should be complete data without NAs")
  NRows <- nrow(y)
  nsamples <- ncol(y)

# Check protein.id
  if(anyNA(protein.id)) stop("NAs not allowed in protein.id")
  if(!identical(NRows,length(protein.id))) stop("length of protein.id must agree with row dimension of y")

# Reorder by protein.id if necessary
# Otherwise, assume to be already ordered
  if(reorder) {
    o <- order(protein.id)
    protein.id <- protein.id[o]
    y <- y[o,,drop=FALSE]
  }

# Count peptides per protein
  npeptides <- drop(rowsum(rep_len(1L,NRows),protein.id,reorder=FALSE))
  nproteins <- length(npeptides)
  df.residual <- (nsamples-1)*(npeptides-1)

# Sweep out row means
  y <- y - rowMeans(y)

# Sweep out protein-wise column means
  ProtColMeans <- rowsum(y,protein.id,reorder=FALSE) / npeptides
  i.protein <- rep(seq_len(nproteins),npeptides)
  y <- y-ProtColMeans[i.protein,]

# Protein-wise residual variances
  s2 <- rowSums(rowsum(y^2,protein.id,reorder=FALSE)) / df.residual

  list(s2=s2, df.residual=df.residual)
}
