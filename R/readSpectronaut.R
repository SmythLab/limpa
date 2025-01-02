readSpectronaut <- function(file="Report.tsv", path=NULL, sep="\t", log=TRUE, run.column = "R.Raw File Name", qty.column = "EG.TotalQuantity", q.columns = c("EG.Qvalue", "PG.Qvalue"), q.cutoffs = 0.01)
# Read normal table output from Spectronaut. 
# Gordon Smyth and Mengbo Li
# Created 18 December 2023. Last modified 23 November 2024.
{
  # Read Spectronaut report file
  if (!is.null(path)) file <- file.path(path, file)
  Select <- c(run.column, "PG.ProteinAccessions", "PG.UniProtIds", "PG.Genes", "EG.PrecursorId", qty.column, "EG.IsImputed", q.columns)
  Report <- fread(file, sep = "\t", select = Select)
  colnames(Report)[1] <- "Run"
  colnames(Report)[6] <- "EG.TotalQuantity"

  # Filter by imputed
  Report <- Report[Report$EG.IsImputed == FALSE, ]

  # Filter by q-values
  if (length(q.columns) > 0L) {
    q.columns <- q.columns[q.columns %in% colnames(Report)]
    if (!identical(length(q.cutoffs), length(q.columns))) {
      q.cutoffs <- rep(q.cutoffs[1], length(q.columns))
      message("Length of q-value columns does not match with length of q-value cutoffs. Use q.cutoffs[1] for all columns.")
    }
    kp <- rep(TRUE, nrow(Report))
    for (qcol in seq_along(q.columns)) {
      kp <- kp & Report[[q.columns[qcol]]] <= q.cutoffs[qcol]
    }
    Report <- Report[kp, ]
  }

  # Convert intensities to wide format
  Samples <- unique(Report$Run)
  Precursors <- unique(Report$EG.PrecursorId)
  y <- matrix(0, length(Precursors), length(Samples))
  mSample <- match(Report$Run, Samples)
  mPrecursor <- match(Report$EG.PrecursorId, Precursors)
  i <- mPrecursor + (mSample - 1L) * length(Precursors)
  y[i] <- Report$EG.TotalQuantity
  colnames(y) <- Samples
  rownames(y) <- Precursors

  # Precursor annotation in wide format
  d <- duplicated(Report$EG.PrecursorId)
  Genes <- data.frame(Report[!d, 2:5])
  row.names(Genes) <- Precursors

  # Output either unlogged EListRaw (with zeros) or logged Elist (with NAs)
  if (log) {
    y[y <= 1L] <- NA
  # Remove rows that are missing in all samples
    ind <- rowMeans(is.na(y)) < 1
    y <- y[ind, , drop = FALSE]
    Genes <- Genes[ind, ]
  # Log
    y <- log2(y)
    new("EList", list(E = y, genes = Genes))
  } else {
    new("EListRaw", list(E = y, genes = Genes))
  }
}
