f_lta_nstarts <- function(S, counts, D, nstarts, tol, maxiter, pdGH)
{
  out <- f_lta(S, counts, D, tol, maxiter, pdGH)
  if(nstarts > 1){ 
    for(i in 2:nstarts){
      out1 <- f_lta(S, counts, D, tol, maxiter, pdGH)
      if(out1$LL > out$LL) out <- out1
    }
  }
  out
}
