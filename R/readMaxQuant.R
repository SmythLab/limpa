readMaxQuant <- function(
    file="peptides.txt", path=NULL, sep="\t", log=TRUE,
    peptide.column = "Sequence",
    qty.column = NULL,
    qty.column.key = "Intensity ", 
    q.columns = c("PEP"), q.cutoffs = 0.01,
    extra.columns = c("Proteins", "Leading razor protein", "Gene names", "Unique (Groups)", "Protein group IDs")
)
  # Read peptides.txt file from MaxQuant. 
  # Gordon Smyth and Mengbo Li
  # Created 24 June 2025. Last modified 15 October 2025.
{
  # Read peptides.txt file
  if (!is.null(path)) file <- file.path(path, file)
  
  # Read column names
  all.columns <- suppressWarnings(fread(file, sep = sep, nrows = 0L, showProgress = FALSE, data.table=FALSE))
  all.columns <- colnames(all.columns)
  
  # get qty.columns by key
  if (!is.null(qty.column.key)) {
    message("Searching for column names matching to `qty.column.key`.")
    qty.column <- all.columns[grepl(qty.column.key, all.columns)]
  }
  
  Select <- c(peptide.column, qty.column, q.columns, extra.columns)
  if (any(!(Select %in% all.columns))) {
    no.in.Select <- setdiff(Select, all.columns)
    message(paste("Columns", paste(no.in.Select, collapse = ","), "not in data!", sep = " "))
    message("Reading the rest of the columns only.")
  }
  Select <- intersect(Select, all.columns)
  extra.columns <- intersect(extra.columns, all.columns)
  
  # Read in all required columns
  Report <- fread(file, sep = sep, select = Select, 
  data.table = FALSE, showProgress = FALSE, integer64="double")
  colnames(Report)[which(colnames(Report) == peptide.column)] <- "Peptide.Sequence"
  
  # Filter by q-values
  if (length(q.columns) > 0L) {
    q.columns <- q.columns[q.columns %in% colnames(Report)]
    if (!identical(length(q.cutoffs), length(q.columns))) {
      q.cutoffs <- rep(q.cutoffs[1], length(q.columns))
      message("Length of q-value columns does not match with length of q-value cutoffs.
      Use q.cutoffs[1] for all columns.")
    }
    kp <- rep(TRUE, nrow(Report))
    for (qcol in seq_along(q.columns)) {
      kp <- kp & Report[[q.columns[qcol]]] <= q.cutoffs[qcol]
    }
    Report <- Report[kp, ]
  }

  # Get the intensity matrix
  y <- as.matrix(Report[, qty.column, drop = FALSE])
#  y <- lapply(y, function(col)
#    vapply(col, as.numeric, numeric(1))
#  )
#  y <- as.data.table(y)
  row.names(y) <- Report$Peptide.Sequence
  if (!is.null(qty.column.key)) {
    colnames(y) <- gsub(qty.column.key, "", colnames(y), fixed = TRUE)
  }
  
  # Precursor annotation in wide format
  Genes <- Report[, extra.columns, drop = FALSE]
  colnames(Genes) <- extra.columns
  row.names(Genes) <- Report$Peptide.Sequence
  
  # Output either unlogged EListRaw (with zeros) or logged Elist (with NAs)
  if (log) {
    y[y == 0L] <- NA
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