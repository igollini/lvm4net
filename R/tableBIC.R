tableBIC <- function(out){
  lout <- length(out) - 1
  bicG <- numeric(lout)
  names(bicG) <- names(out[- (lout + 1)])
  
  for(i in 1:lout) bicG[i] <- out[[i]]$BIC
  
  resBIC <- vector('list', 2)
  names(resBIC) <- c('Model Selection', 'Table of BIC Results')
  resBIC[[1]] <- names(bicG)[bicG == min(bicG)]
  names(resBIC[[1]]) <- 'Model with lower BIC:'
  resBIC[[2]] <- bicG
  resBIC
}
