f_mlta_methods <-
  function(S, counts, G, D, nstarts, tol, maxiter, pdGH, wfix)
  {
    if (D == 0) {
      if (any(G == 1)) {
        out <- f_lca_nstarts(S, counts, G, nstarts, tol, maxiter)
      } else{
        if (length(G) == 1) {
          out <- f_lca_nstarts(S, counts, G, nstarts, tol, maxiter)
        } else{
          out <- vector("list", length(G) + 1)
          names(out) <- c(paste('G', G, sep = '='), 'BIC')
          i <- 0
          for (g in G) {
            i <- i + 1
            out[[i]] <- f_lca_nstarts(S, counts, g, nstarts, tol, maxiter)
          }
          out[[length(G) + 1]] <- tableBIC(out)
        }
      }
    } else {
      if (D > 0 && G == 1){
        if(length(D) == 1){ 
          out <- f_lta_nstarts(S, counts, D, nstarts, tol, maxiter, pdGH)
          out$eta <- 1
        }else{
          out<-vector("list", length(D) + 1)
          names(out) <- c(paste('Dim y', D, sep = '='), 'BIC')
          i<-0
          for(diy in D){
            i <- i + 1
            out[[i]] <- f_lta_nstarts(S, counts, diy, nstarts, tol, maxiter, pdGH)
            out[[i]]$eta <- 1
          }
          
          cat('BIC results',"\n \n")
          out[[length(D) + 1]]<-tableBIC(out)
        }
      }
      if (D > 0 && G > 1) {
        if (wfix == TRUE) {
          out <- f_mlta_nstarts_wfix(S, counts, G, D, nstarts, tol, maxiter, pdGH)
          class(out) <- c("mlta")
          
        } else{
          out <- f_mlta_nstarts(S, counts, G, D, nstarts, tol, maxiter, pdGH)
        }
      }
    }
    class(out) <- c("mlta")
    return(out)
  }
