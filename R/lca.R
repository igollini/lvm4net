#' Latent Class Analysis
#'
#' Latent class analysis (LCA) can be used to find groups in the sender nodes (with the condition of independence within the groups). For more details see Gollini, I. (in press) and Gollini, I., and Murphy, T. B. (2014).
#'
#' @param X (\code{N} x \code{M}) binary incidence matrix
#' @param G number of groups
#' @param nstarts integer number of different starts for the EM algorithm.  Default \code{nstarts = 3}.
#' @param tol desired tolerance for convergence. Default \code{tol = 0.1^2}
#' @param maxiter maximum number of iterations. Default \code{maxiter = 500}
#'
#' @return List containing the following information for each model fitted:
#' \itemize{
#' \item \code{p} (\code{G} x \code{M}) matrix containing the conditional probability of observing a link to sender nodes if the receiver nodes are from group g.
#' \item \code{eta} \eqn{\eta_g} is the mixing proportion for the group \eqn{g (g = 1,..., G)}, that corresponds to the prior probability that a randomly chosen sender node is in the g-th group. 
#' \item \code{z} (\code{N} x \code{G}) matrix containing posterior probability for each sender node to belong to each group
#' \item \code{LL} log likelihood
#' \item \code{BIC} Bayesian Information Criterion (BIC) (Schwarz (1978))
#' }
#' If multiple models are fitted the output contains also a table to compare the BIC for all models fitted.
#' 
#' @seealso \code{\link{mlta}}
#' @references Gollini, I. (in press) 'A mixture model approach for clustering bipartite networks', Challenges in Social Network Research Volume in the Lecture Notes in Social Networks (LNSN - Series of Springer). Preprint: \url{https://arxiv.org/abs/1905.02659}.
#' @references Gollini, I., and Murphy, T. B. (2014), 'Mixture of Latent Trait Analyzers for Model-Based Clustering of Categorical Data', Statistics and Computing, 24(4), 569-588 \url{http://arxiv.org/abs/1301.2167}.
#' @export
#' @examples
#' ### Simulate Bipartite Network
#' set.seed(1)
#' X <- matrix(rbinom(4 * 12, size = 1, prob = 0.4), nrow = 12, ncol = 4)
#'
#' resLCA <- lca(X, G = 2:3)
lca <- function(X, G, nstarts = 3, tol = 0.1^2, maxiter = 250) {
  
  if (any(G < 1)) {
    print("Specify G > 0!")
    return("Specify G > 0!")
  }
  
  XS <- XtoS(X)
  S <- XS$S
  counts <- XS$counts
  
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
  out
}
