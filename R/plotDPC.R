plotDPC <- function(dpcfit,
                    add.jitter = TRUE,
                    jitter.amount = NULL,
                    point.cex = 0.2,
                    lwd = 2,
                    ylim = c(0, 1), 
                    main = "Detection probability curve", 
                    ...)
# Plot detection probability curve (DPC) produced by dpc().
# Mengbo Li
# Created 16 May 2022 as part of protDP package.
# Migrated to limpa 11 Sept 2024. Last modified 28 Dec 2024.
{
  x <- (dpcfit$mu.obs + dpcfit$mu.mis)/2
  y <- dpcfit$prop.detected
  if (add.jitter) y <- jitter(y, amount = jitter.amount)
  
  plot(x = x, y = y,
       pch = 16, cex = point.cex, 
       ylim = ylim,
       xlab = "Intensity", ylab = "Detection probability", 
       main = main, ...)
 
   x <- x[order(x)]
   lines(x = x, y = plogis(dpcfit$dpc.start[1] + dpcfit$dpc.start[2]*x), lty = "dashed", lwd = lwd)
   lines(x = x, y = plogis(dpcfit$dpc[1] + dpcfit$dpc[2]*x), col = "blue", lwd = lwd)
   legend("bottomright", legend = c("Start", "Final"), col = c("black", "blue"), lty = c(2, 1), lwd = lwd, cex = 0.8)

  invisible(list(x=x,y=y))
}