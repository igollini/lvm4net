f_lca <-function(S, counts, G, tol, maxiter){
  
  Ns <- nrow(S)
  N <- sum(counts)
  M <- ncol(S) 
  rownames(S) <- NULL
  
  # Initialize EM Algorithm
  eta <- rep(1/G, G)
  lab <- sample(1:G, prob = eta, size = Ns, replace = TRUE)
  z <- counts * cls(lab)
  
  # Preliminary M-step
  
  eta <- colSums(z) / N
  
  # Set up
  
  p <- (t(z) %*% S) / (N * eta)
  
  # Initialize Aitken acceleration
  
  l <- numeric(3)
  lA <- numeric(2)
  
  # Iterative process
  
  v <- matrix(0, Ns, G)
  ll <- -Inf
  diff <- 1
  iter <- 0
  
  while(diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    
    # E-step
    
    v <- apply(cbind(eta, p), 1, fvg, S)
    v[is.nan(v)] <- 0
    v[v < 0] <- 0
    vsum <- rowSums(v)
    z <- counts * v / vsum
    ll <- counts %*% log(vsum)
    
    # M-step
    
    eta <- colSums(z)/N
    
    p <- (t(z) %*% S) / (N * eta)
    
    # Aitken acceleration
    
    l <- c(l[-1], ll)
    Out <- AitkenAcc(l, lA, iter)
    diff <- Out[1]
    lA <- Out[-1]
  }
  
  p <- (t(z) %*% S) / (N * eta)	# pi_gm proportion of positive response in group g
  
  BIC <- -2*ll + (G *(M + 1) - 1) * log(N)
  
  expected <- vsum * N
  
  colnames(p) <- NULL
  rownames(p) <- rownames(p, do.NULL = FALSE, prefix = "Group ")
  colnames(p) <- colnames(p, do.NULL = FALSE, prefix = "p_g")
  eta<-matrix(eta, byrow=T, G)
  rownames(eta) <- rownames(eta, do.NULL = FALSE, prefix = "Group ")
  colnames(eta) <- "eta"
  names(ll)<- "Log-Likelihood:"
  names(BIC)<-"BIC:"
  
  out <- list(p = p, eta = eta, LL = ll, BIC = BIC)
  list(p = p, eta = eta, LL = ll, BIC = BIC, z = z, expected = expected)
  
}
