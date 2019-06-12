#' Mixture of Latent Trait Analyzers
#'
#' Mixture of latent trait analyzers (MLTA) has been introduced by Gollini and Murphy (2014) and Gollini (in press) to identify groups assuming the existence of a latent trait describing the dependence structure between receiver nodes within groups of sender nodes and therefore capturing the heterogeneity of sender nodes' behaviour within groups. The function \code{mlta} makes use of a variational inferential approach. For more details see Gollini, I. (in press) and Gollini, I., and Murphy, T. B. (2014).
#'
#' @param X (\code{N} x \code{M}) binary incidence matrix
#' @param G number of groups
#' @param D dimension of the continuous latent variable
#' @param wfix Logical. Fit the parsiomonius model with the w parameters equal across groups. Default \code{wfix = FALSE}
#' @param nstarts number of starts.  Default \code{nstarts = 3}
#' @param tol desired tolerance for convergence. Default \code{tol = 0.1^2}
#' @param maxiter maximum number of iterations. Default \code{maxiter = 500}
#' @param pdGH number of quadrature points for the Gauss-Hermite quadrature. Default \code{pdGH = 21}
#'
#' @return List containing the following information for each model fitted:
#' \itemize{
#' \item \code{b} matrix containing intercepts for the logistic response function
#' \item \code{w} array containing slopes for the logistic response function
#' \item \code{eta} \eqn{\eta_g} is the mixing proportion for the group \eqn{g (g = 1,..., G)}, that corresponds to the prior probability that a randomly chosen sender node is in the g-th group. 
#' \item \code{mu} (\code{N} x \code{D} x \code{G}) array containing posterior means for the latent variable
#' \item \code{C} (\code{D} x \code{D} x \code{N} x \code{G}) array containing posterior variances for the latent variable
#' \item \code{z} (\code{N} x \code{G}) matrix containing posterior probability for each sender node to belong to each group
#' \item \code{LL} log likelihood
#' \item \code{BIC} Bayesian Information Criterion (BIC) (Schwarz (1978))
#' }
#' If multiple models are fitted the output contains also tables to compare the log likelihood and BIC for all models fitted.
#' 
#' @seealso \code{\link{lta}} \code{\link{lca}}
#' @references Gollini, I. (in press) 'A mixture model approach for clustering bipartite networks', Challenges in Social Network Research Volume in the Lecture Notes in Social Networks (LNSN - Series of Springer). Preprint: \url{https://arxiv.org/abs/1905.02659}.
#' @references Gollini, I., and Murphy, T. B. (2014), 'Mixture of Latent Trait Analyzers for Model-Based Clustering of Categorical Data', Statistics and Computing, 24(4), 569-588 \url{http://arxiv.org/abs/1301.2167}.
#' @export
#' @examples
#' ### Simulate Bipartite Network
#' set.seed(1)
#' X <- matrix(rbinom(4 * 12, size = 1, prob = 0.4), nrow = 12, ncol = 4)
#'
#' resMLTA <- mlta(X, G = 2, D = 1)
mlta <-
  function(X,
    G,
    D,
    wfix = FALSE,
    nstarts = 3,
    tol = 0.1 ^ 2,
    maxiter = 250,
    pdGH = 21)
  {
    if (any(G < 1)) {
      print("Specify G > 0!")
      return("Specify G > 0!")
    }
    
    if (any(D < 0)) {
      print("Specify D >= 0!")
      return("Specify D >= 0!")
    }
    
    ############################################
    
    #XS <- XtoS(X)
    #S <- XS$S
    #counts <- XS$counts
    
    S <- X
    counts <- rep(1, nrow(X))
    
    if (length(D) == 1 && length(G) == 1) {
      out <-
        f_mlta_methods(S, counts, G, D, nstarts, tol, maxiter, pdGH, wfix)
    } else{
      out <- vector("list", length(D) * length(G) + 3)
      names(out) <- c(t(outer(
        paste('G=', G, sep = ''),
        paste('dim y=', D, sep = ''),
        paste,
        sep = ','
      )),
        'BIC', 'LL', 'LLva')
      
      bictab <- matrix(0, length(D), length(G))
      lltab <- matrix(0, length(D), length(G))
      
      rownames(bictab) <- paste('dim y=', D, sep = '')
      colnames(bictab) <- paste('G=', G, sep = '')
      
      rownames(lltab) <- paste('dim y=', D, sep = '')
      colnames(lltab) <- paste('G=', G, sep = '')
      
      llvatab <- lltab
      
      i <- 0
      for (g in G){
        for (diy in D) {
          i <- i + 1
          
          out[[i]] <-
            f_mlta_methods(S, counts, g, diy, nstarts, tol, maxiter, pdGH, wfix)
          
          bictab[i] <- out[[i]]$BIC
          lltab[i] <- out[[i]]$LL
          
          if (diy == 0) {
            llvatab[i] <- out[[i]]$LL
          } else{
            llvatab[i] <- out[[i]]$LLva
          }
        }
      
      out[[length(G) * length(D) + 1]] <-
        ResTable(bictab, restype = 'BIC')
      out[[length(G) * length(D) + 2]] <-
        ResTable(lltab, restype = 'll')
      out[[length(G) * length(D) + 3]] <-
        ResTable(llvatab, restype = 'llva')
      }
      class(out) <- "mmlta"
    }
    
    out
  }
