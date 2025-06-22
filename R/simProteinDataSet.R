simProteinDataSet <- function(n.peptides=100, n.groups=2, samples.per.group=5, peptides.per.protein=4, mu.range=c(2,10), sigma=0.4, prop.de=0.2, fc=2, dpc.intercept=NULL, dpc.slope=0.7, prop.missing=0.4)
# Simulate peptide data with missing values.
# Created 20 Dec 2024. Last modified 26 May 2025.
{
# Setup peptide-wise mu and sigma
  mu <- seq(from=mu.range[1], to=mu.range[2], length.out=n.peptides)
  sigma <- rep(0.6, n.peptides)

# Choose dpc intercept to give roughly the required proportion of NAs
  if(is.null(dpc.intercept)) {
    mean.mu <- mean(mu.range)
    eta.prop.missing <- qlogis(prop.missing, lower.tail=FALSE)
    dpc.intercept <- eta.prop.missing - dpc.slope * mean.mu
#   Two Newton iterations polish
    gq <- gauss.quad.prob(16,l=mu.range[1],u=mu.range[2])
    eta <- dpc.intercept + dpc.slope * gq$nodes
    dpc.intercept <- dpc.intercept - (sum(gq$weights*plogis(eta)) - (1-prop.missing) ) / sum(gq$weights*dlogis(eta))
    eta <- dpc.intercept + dpc.slope * gq$nodes
    dpc.intercept <- dpc.intercept - (sum(gq$weights*plogis(eta)) - (1-prop.missing) ) / sum(gq$weights*dlogis(eta))
  }

# Generate y
  y <- rnorm(n.peptides*n.groups*samples.per.group, mean=mu, sd=sigma)
  dim(y) <- c(n.peptides, n.groups*samples.per.group)

# Peptide names
  w <- floor(log10(n.peptides)+1)
  peptide.i <- seq_len(n.peptides)
  rownames(y) <- paste0("Peptide",formatC(peptide.i,width=w,flag=0))

# Sample names
  colnames(y) <- paste0("S",seq_len(ncol(y)))

# Protein names
  n.proteins <- ceiling(n.peptides / peptides.per.protein)
  protein.i <- rep(seq_len(n.proteins),each=peptides.per.protein,length.out=n.peptides)
  w <- floor(log10(n.proteins)+1)
  Protein <- paste0("Protein",formatC(protein.i,width=w,flag=0))

# Add DE with 2-fold change between first two groups
# All peptides in a DE protein are DE
  DE.Status <- rep_len("NotDE",n.peptides)
  if(n.groups > 1) {
    log2FC <- log2(fc)/2
    n.DE.prot <- round(prop.de*n.proteins)
    if(n.DE.prot > 0.5) {
      is.DE.prot <- sample(seq_len(n.proteins), n.DE.prot)
      n.Up.prot <- ceiling(n.DE.prot/2)
      is.Up.prot <- is.DE.prot[seq_len(n.Up.prot)]
      is.Up.pep <- peptide.i[protein.i %in% is.Up.prot]
      DE.Status[is.Up.pep] <- "Up"
      j1 <- seq_len(samples.per.group)
      j2 <- j1 + samples.per.group
      y[is.Up.pep,j1] <- y[is.Up.pep,j1] - log2FC
      y[is.Up.pep,j2] <- y[is.Up.pep,j2] + log2FC
      n.Dn.prot <- n.DE.prot - n.Up.prot
      if(n.Dn.prot > 0.5) {
        is.Dn.prot <- is.DE.prot[(n.Up.prot+1):n.DE.prot]
        is.Dn.pep <- peptide.i[protein.i %in% is.Dn.prot]
        DE.Status[is.Dn.pep] <- "Down"
        y[is.Dn.pep,j1] <- y[is.Dn.pep,j1] + log2FC
        y[is.Dn.pep,j2] <- y[is.Dn.pep,j2] - log2FC
      }
    }
  }

# Add NAs
  probmis <- plogis(dpc.intercept + y*dpc.slope, lower.tail=FALSE)
  ismis <- rbinom(length(y), prob=probmis, size=1)
  y.complete <- y
  y[ismis==1] <- NA

# EList
  Genes <- data.frame(Protein=Protein,DEStatus=DE.Status)
  row.names(Genes) <- rownames(y)
  Group <- rep(seq_len(n.groups),each=samples.per.group)
  targets <- data.frame(Group=Group)
  row.names(targets) <- colnames(y)
  new("EList", list(E=y, genes=Genes, targets=targets, other=list(E.complete=y.complete)))
}
