print.mmlta <- function(x){
  cat("Log-Likelihood:\n")
      print(x$LL$`Table of LL (G-H Quadrature correction)`)
      cat("BIC:\n")
      print(x$BIC$`Table of BIC Results`)
      print(x$BIC$`Model Selection`)
}
