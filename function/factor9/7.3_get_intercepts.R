get_intercepts <- function(pls){
  endogs <- unique(pls$smMatrix[, "target"])
  sapply(endogs, function(y){
    preds <- pls$smMatrix[pls$smMatrix[, "target"] == y, "source"]
    dat   <- as.data.frame(pls$fscores[, c(y, preds), drop = FALSE])
    fit   <- lm(reformulate(preds, response = y), data = dat)
  })
  return(coef(fit))
}

