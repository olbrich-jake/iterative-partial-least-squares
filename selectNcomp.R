selectNcomp<-function (object, method = c("randomization", "onesigma"), nperm = 999, 
          alpha = 0.01, ncomp = object$ncomp, plot = FALSE, ...) 
{
  if (!isTRUE(object$validation$method == "CV")) 
    stop("No cross-validation data available in model")
  if (dim(residuals(object))[2] > 1) 
    stop("Only univariate response supported")
  rmseps <- c(RMSEP(object, "CV")$val)
  maxIdx <- ncomp + 1
  absBest <- which.min(rmseps[seq_len(maxIdx)])
  if (absBest > 0) {
    method <- match.arg(method)
    if (is.null(origResponse <- object$y)) 
      origResponse <- c(predict(object, ncomp = 1) + residuals(object)[,1, 1])
    allresids <- cbind(origResponse - (sum(origResponse) - 
                                         origResponse)/(length(origResponse) - 1), object$validation$pred[,1,] - origResponse)
    if (method == "randomization") {
      pvals <- sapply(seq_len(absBest - 1), function(ii) randomiz.test(allresids[,ii], allresids[,absBest], nperm = nperm))
      idx <- which(pvals > alpha)
      selection <- min(c(idx, absBest)) - 1
    }
    else {
      residsds <- apply(allresids, 2, sd)/sqrt(nrow(allresids))
      uls <- rmseps - residsds
      selection <- min(which(uls < rmseps[absBest])) - 1
    }
    if (plot) {
      xvals <- seq_along(rmseps) - 1
      plot(xvals, rmseps, ylab = "RMSEP", xlab = "Number of components", type = "b", ...)
      if (method == "onesigma") {
        arrows(xvals, uls, xvals, rmseps + residsds, code = 3, col = "gray", angle = 90, length = 0.1)
      }
      else {
        points(xvals[idx], rmseps[idx], cex = 2, col = 4)
      }
      abline(h = rmseps[absBest], col = "gray", lty = 3)
      abline(v = absBest - 1, col = "gray", lty = 3)
      abline(v = selection, col = "blue", lty = 2)
      legend("topright", legend = c("Abs. minimum", "Selection"),lty = 3:2, col = c("gray", "blue"), bg = "white")
    }
    selection
  }
  else {
    warning("Lowest CV error found at 0 components, no testing performed")
    0
  }
}