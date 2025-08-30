get_intercepts <- function(pls){
  endogs <- unique(pls$smMatrix[, "target"])
  preds <- pls$smMatrix[pls$smMatrix[, "target"] == endogs, "source"]
  dat   <- as.data.frame(pls$fscores[, c(endogs, preds), drop = FALSE])
  fit   <- lm(reformulate(preds, response = endogs), data = dat)
  return(coef(fit))
}

