AitkenAcc <- function(l, lA, iter){
  a <- (l[3] - l[2]) / (l[2] - l[1])
  lA[1] <- lA[2]
  
  if (a > 0 & a < 1) 
  {
    lA[2] <- l[2] + (l[3] - l[2]) / (1 - a)
    }else{
      lA[2] <- l[3]
    }
  
  diff <- abs(lA[2] - lA[1])

  if ((l[3] < l[2]) & (iter > 5))
  {
    stop("Decrease in log-likelihood")
  }	
  
  Out <- c(diff, lA)
  Out
}