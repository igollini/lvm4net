print.mlta <- function(x){
  stopifnot(inherits(x, 'mlta'))
  cat("b:\n")
    print(x$b) 
    cat("\nw:\n")  
    print(x$w)
    cat("\neta:\n")
    print(x$eta)
    cat("\n")
    print(x$LL)
    print(x$BIC) 
}
