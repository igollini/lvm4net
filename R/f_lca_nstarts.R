f_lca_nstarts <- function(S, counts, G, nstarts, tol, maxiter) {
  
  out <- f_lca(S, counts, G, tol, maxiter)
  
  for(i in 2:nstarts){
    out1 <- f_lca(S, counts, G, tol, maxiter)
    if(out1$LL > out$LL) out <- out1
  }
  out
}
