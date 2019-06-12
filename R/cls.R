cls<- function(lab) 
{
  gr <- sort(unique(lab))
  z <- matrix(0, length(lab), length(gr))
  for (i in 1:length(gr)) z[lab == gr[i], i] <- 1
  z
}