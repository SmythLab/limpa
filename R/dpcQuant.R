dpcQuant <- function(y, ...)
  UseMethod("dpcQuant")

dpcQuant.default <- function(y, protein.id, dpc=NULL, dpc.slope=0.8, verbose=TRUE, chunk=1000L, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 31 Dec 2024. Last modified 23 Mar 2025.
{
# Check y
  y <- as.matrix(y)

# Check protein.id
  if(missing(protein.id)) stop("Need protein ID")
  if(!identical(nrow(y),length(protein.id))) stop("length(protein.id) must match nrows(y)")

# Construct EList and pass to EList method
  z <- new("EList",list(E=y,genes=data.frame(Protein=protein.id)))
  dpcQuant(z,protein.id=protein.id,dpc=dpc,dpc.slope=dpc.slope,verbose=verbose,chunk=chunk,...)
}

dpcQuant.EList <- function(y, protein.id="Protein.Group", dpc=NULL, dpc.slope=0.8, verbose=TRUE, chunk=1000L, ...)
# Use the DPC to quantify protein expression values by maximum posterior.
# Created 27 Dec 2024. Last modified 19 Jun 2025.
{
# Check dpc
  if(is.list(dpc)) dpc <- dpc$dpc
  if(is.null(dpc)) {
    beta0 <- estimateDPCIntercept(y, dpc.slope=dpc.slope)
    if(verbose) message("DPC intercept estimated as ", formatC(beta0))
    dpc <- c(beta0=beta0, beta1=dpc.slope)
  }

# Get vector of protein IDs
# protein.id can be either an annotation column name or a vector of IDs.
  protein.id <- as.character(protein.id)
  if(identical(length(protein.id),1L)) {
    ColName <- protein.id
    protein.id <- y$genes[[ColName]]
    if(is.null(protein.id)) stop("Column \"",ColName,"\" not found in y$genes")
  } else {
    if(!identical(nrow(y),length(protein.id))) stop("length(protein.id) must match nrows(y)")
  }

# If all proteins have just one peptide, call dpcImpute instead
  a <- anyDuplicated(protein.id)
  if(identical(a,0L)) {
    message("All proteins have exactly one peptide: calling dpcImpute() instead")
    return(dpcImpute(y,dpc,verbose=verbose,chunk=chunk,...))
  }

# Sort peptides in protein order
  o <- order(protein.id)
  protein.id <- protein.id[o]
  y <- y[o,]

# Estimate Bayes hyperparameters
  if(verbose) message("Estimating hyperparameters ...")
  h <- dpcQuantHyperparam(y, protein.id=protein.id, dpc.slope=dpc[2], ...)

# Summarize peptides to proteins by maximum posterior estimation
  if(verbose) message("Quantifying proteins ...")
  y.protein <- peptides2Proteins(y,
                          protein.id = protein.id,
                          dpc = dpc,
                          sigma = h$sigma,
                          prior.mean = h$prior.mean,
                          prior.sd = h$prior.sd,
                          prior.logFC = h$prior.logFC,
                          standard.errors = TRUE,
                          verbose = verbose,
                          chunk = chunk)

# Add back original original annotation
  d <- !duplicated(protein.id)
  if(!is.null(y$genes)) {
    # remove duplicated columns
    extra_anno <- which(!(colnames(y$genes) %in% colnames(y.protein$genes)))
    genes <- y$genes[, extra_anno, drop = FALSE]
    # remove precursor-level columns
    dup_protein.id <- duplicated(protein.id)
    kp_cols <- vapply(seq_len(ncol(genes)), function(ii) {
      dup_ii <- duplicated(genes[, ii])
      all(dup_ii[dup_protein.id])
    }, logical(1L))
    genes <- genes[d, kp_cols, drop = FALSE]
    row.names(genes) <- row.names(y.protein)
    y.protein$genes <- data.frame(genes,y.protein$genes)
  }
  y.protein$targets <- y$targets
  y.protein$dpc <- dpc
  y.protein
}
