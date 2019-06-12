ResTable <- function(bicG, restype)
{
  if (restype == 'll') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL (G-H Quadrature correction)'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'llva') {
    resBIC <- vector('list', 1)
    names(resBIC) <- 'Table of LL (variational approximation)'
    resBIC[[1]] <- bicG
  }
  
  if (restype == 'BIC') {
    resBIC <- vector('list', 2)
    names(resBIC) <- c('Table of BIC Results', 'Model Selection')
    resBIC[[1]] <- bicG
    resBIC[[2]] <-
      paste(colnames(bicG)[t(bicG == min(bicG)) * seq(1:ncol(bicG))], rownames(bicG)[(bicG ==
          min(bicG)) * seq(1:nrow(bicG))], sep = ', ')
    names(resBIC[[2]]) <- 'Model with lower BIC:'
  }
  
  resBIC
}
