buildYYh <- function(YY, mu, dimy, Ns){
  
  YYh <- matrix(NA, (dimy + 1) * (dimy + 1), Ns)
  
  for(ddi in 1:dimy){
    YYh[((ddi - 1) * dimy + ddi):(ddi - 1 + (ddi * dimy)),] <- YY[((ddi - 1) * dimy + 1):(ddi * dimy),]
    YYh[ddi * dimy + ddi,] <- mu[ddi,]
  }
  
  YYh[(dimy * (dimy + 1) + 1):((dimy + 1) * (dimy + 1) - 1),] <- mu[1:ddi,]
  
  YYh[(dimy + 1) * (dimy + 1),] <- 1
  
  YYh		
}