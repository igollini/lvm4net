f_mlta_nstarts <- function(S, counts, G, D, nstarts, tol, maxiter, pdGH)
{
  out <- f_mlta(S, counts, G, D, tol, maxiter, pdGH)
  
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_mlta(S, counts, G, D, tol, maxiter, pdGH)
      if(out1$LL > out$LL) out <- out1
    }
  }
  return(out)
}
