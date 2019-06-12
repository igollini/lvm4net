fvg <- function(etaP, S = S)
  {	
  PS <- t(t(S) * etaP[-1]) + t(t(1 - S)*(1 - etaP[-1]))
  etaP[1] * apply(PS, 1, prod)
}
