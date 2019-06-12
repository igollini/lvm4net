ghLLlta <- function(S, counts, w, b, D, pdGH) {
  b <- as.vector(b)
  Ns <- nrow(S)
  N <- sum(counts)
  M <- ncol(S)
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints^D
  GaHer <- glmmML::ghq(npoints, modified = FALSE)
  ygh <- GaHer$zeros
  wgh <- GaHer$weights
  Ygh <- as.matrix(expand.grid(rep(list(ygh), D)))
  Wghm <- as.matrix(expand.grid(rep(list(wgh), D)))
  Wgh <- apply(Wghm, 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  Agh <- t(t(w) %*% t(Ygh) + b)
  pgh <- 1 / (1 + exp(-Agh))
  
  fxy <- exp(S %*% t(log(pgh)) + (1 - S) %*% t(log(1 - pgh)))
  
  llGH1 <- rowSums(rep(Beta, each = Ns) * fxy)
  llGH <- sum(counts * log(llGH1))
  expected <- llGH1 * N
  
  list(expected = expected, llGH = llGH)
}
