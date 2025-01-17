plotProtein <- function(y, protein, col = "black", cex = 2, lwd = 2, ...)
# Plot the sample-wise protein summary with error bars by dpcQuant. 
{
  #	Check input
	if(!is(y,"EList")) stop("y must be an EList output from dpcquant()")
    if (!identical(length(protein), 1L)) stop("Input 1 protein to plot at each time")
    if (is.character(protein)) {
        if (!(protein %in% rownames(y))) stop("Protein id not found in y, check input")
    }

  # Get protein values and standard errors
    y_protein <- y$E[protein, ]
    se_protein <- y$other$standard.error[protein, ]
    narrays <- ncol(y)

  # Plot
    x <- seq_len(narrays)
    y0 <- y_protein - se_protein
    y1 <- y_protein + se_protein
    ylim <- c(min(y0, y1), max(y0, y1))
    plot(y_protein, pch = 16, ylim = ylim, col = col, xaxt = "none", 
    cex = cex, xlab = "Sample", ylab = "Estimated log-intensity", ...)
    axis(1, at = x, labels = colnames(y))
    arrows(x, y0, x, y1, length = 0.1, angle = 90, code = 3, lwd = lwd, col = col)

}