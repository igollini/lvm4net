#' Lift
#'
#' The lift can be used to analyse the dependence within each groups found using the function \code{\link{mlta}}.The lift can be used to quantify the effect of the dependence on the probability of a sender nodes being liked to two receivers within each group compared to the probability of being liked to two receivers under an independence model. Two independent links to the receiver nodes have lift = 1: the more the links to receiver nodes are dependent, the further the value of the lift is from 1. 
#' @param x object of class \code{mlta}
#' @param pdGH number of quadrature points for the Gauss-Hermite quadrature. Default \code{pdGH = 21}
#'
#' @return The function returns an (\code{M} x \code{M} x \code{D}) array.
#' @seealso \code{\link{mlta}}
#' @references Gollini, I. (in press) 'A mixture model approach for clustering bipartite networks', Challenges in Social Network Research Volume in the Lecture Notes in Social Networks (LNSN - Series of Springer). Preprint: \url{https://arxiv.org/abs/1905.02659}.
#' @references Gollini, I., and Murphy, T. B. (2014), 'Mixture of Latent Trait Analyzers for Model-Based Clustering of Categorical Data', Statistics and Computing, 24(4), 569-588 \url{http://arxiv.org/abs/1301.2167}.
#' @export
#' @examples
#' ### Simulate Bipartite Network
#' set.seed(1)
#' X <- matrix(rbinom(4 * 12, size = 1, prob = 0.4), nrow = 12, ncol = 4)
#' res <- mlta(X, G = 2, D = 1)
#' res_lift <- lift(res)
lift <- function(x, pdGH = 21)
{
  stopifnot(inherits(x, "mlta"))
  
  z <- x$z
  w <- x$w
  b <- x$b
  
  Ns <- nrow(z)
  D <- nrow(w)
  M <- ncol(b)
  G <- nrow(b)
  
  if (length(dim(w)) == 2)
    w <- array(w, c(dim(w), G))
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <-
    apply(as.matrix(expand.grid(rep(
      list(GaHer$weights), D
    ))), 1, prod) * apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  Agh <- matrix(NA, ny, M)
  pgh <- matrix(NA, ny, M)
  fxy <- array(0, c(Ns, ny, G))
  
  px1 <- rep(NA, M)
  px1x1 <- matrix(0, M, M)
  lift <- array(NA, c(M, M, G))
  
  for (g in 1:G)
  {
    if (D == 1) {
      Agh <- t(tcrossprod(w[, , g], Ygh) + b[g, ])
    } else{
      Agh <- t(tcrossprod(t(w[, , g]), Ygh) + b[g, ])
    }
    pgh <- 1 / (1 + exp(-Agh))
    px1 <- colSums(Beta * pgh)
    px1x1 <- crossprod(sqrt(Beta) * pgh)
    
    Mat <- (px1x1 / tcrossprod(px1)) * upper.tri(matrix(1, M, M))
    Mat[lower.tri(Mat)] <- NA
    diag(Mat) <- NA
    
    lift[, , g] <-  Mat
  }
  
  rownames(lift) <- rownames(lift, do.NULL = FALSE, prefix = "")
  colnames(lift) <- colnames(lift, do.NULL = FALSE, prefix = "")
  dimnames(lift)[[3]] <- paste("g =", seq(1, G) , sep = ' ')
  
  lift
}
