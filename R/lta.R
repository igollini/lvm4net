#' Latent Trait Analysis
#'
#' Latent trait analysis (LTA) can be used to model the dependence in the receiver nodes by using a continuous D-dimensional latent variable. The function \code{lta} makes use of a variational inferential approach. For more details see Gollini, I. (in press) and Gollini, I., and Murphy, T. B. (2014).
#'
#' @param X (\code{N} x \code{M}) binary incidence matrix
#' @param D dimension of the continuous latent variable
#' @param nstarts number of starts.  Default \code{nstarts = 3}
#' @param tol desired tolerance for convergence. Default \code{tol = 0.1^2}
#' @param maxiter maximum number of iterations. Default \code{maxiter = 500}
#' @param pdGH number of quadrature points for the Gauss-Hermite quadrature. Default \code{pdGH = 21}
#'
#' @return List containing the following information for each model fitted:
#' \itemize{
#' \item \code{b} intercepts for the logistic response function
#' \item \code{w} slopes for the logistic response function
#' \item \code{mu} (\code{N} x \code{D}) matrix containing posterior means for the latent variable
#' \item \code{C} list of \code{N} (\code{D} x \code{D}) matrices containing posterior variances for the latent variable
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
#' resLTA <- lta(X, D = 1:2) 
lta <- function(X, D, nstarts = 3, tol = 0.1^2, maxiter = 250, pdGH = 21) {
  
  if(any(D == 0)) stop("D must be > 0")

  XS <- XtoS(X)
  S <- XS$S
  counts <- XS$counts
  
	if(length(D) == 1){ 
	  out <- f_lta_nstarts(S, counts, D, nstarts, tol, maxiter, pdGH)
		}else{
		  out<-vector("list", length(D) + 1)
		  names(out) <- c(paste('Dim y', D, sep = '='), 'BIC')
			i<-0
			for(diy in D){
			  i <- i + 1
			  out[[i]] <- f_lta_nstarts(S, counts, diy, nstarts, tol, maxiter, pdGH)
			}
			out[[length(D) + 1]]<-tableBIC(out)
		}
	out
}
