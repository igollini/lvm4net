#' Latent Variable Models for Networks
#'
#' \code{lvm4net} provides a range of tools for latent variable models for
#' network data. Most of the models are implemented using a fast
#' variational inference approach.
#' 
#' Latent space models for one-mode binary networks: the function \code{\link{lsm}} implements the latent space model (LSM) introduced by  Hoff et al. (2002) using variational inference and squared Euclidian distance; the function
#' \code{\link{lsjm}} implements latent space joint model (LSJM) for multiplex networks introduced by
#' Gollini and Murphy (2016).
#' These models assume that each node of a network has a latent position
#' in a latent space: the closer two nodes are in the latent space, the more likely
#' they are connected.
#' 
#' Latent variable models for binary bipartite networks: the function \code{\link{lca}} implements the latent class analysis (LCA) to find groups in the sender nodes (with the condition of independence within the groups); the function \code{\link{lta}} implements the latent trait analysis (LTA) to model the dependence in the receiver nodes by using a continuous latent variable; the function \code{\link{mlta}} implements the mixture of latent trait analyzers (MLTA) introduced by Gollini and Murphy (2014) and Gollini (in press) to identify groups assuming the existence of a latent trait describing the dependence structure between receiver nodes within groups of sender nodes and therefore capturing the heterogeneity of sender nodes' behaviour within groups. \code{\link{lta}} and \code{\link{mlta}} use variational inference.
#' @references Gollini, I. (in press) 'A mixture model approach for clustering bipartite networks', Challenges in Social Network Research Volume in the Lecture Notes in Social Networks (LNSN - Series of Springer). Preprint: \url{https://arxiv.org/abs/1905.02659}.
#' @references Gollini, I., and Murphy, T. B. (2014), 'Mixture of Latent Trait Analyzers for Model-Based Clustering of Categorical Data', Statistics and Computing, 24(4), 569-588 \url{http://arxiv.org/abs/1301.2167}.
#' @references Gollini, I., and Murphy, T. B. (2016), 'Joint Modelling of Multiple Network Views', Journal of Computational and Graphical Statistics, 25(1), 246-265 \url{http://arxiv.org/abs/1301.3759}.
#' @references Hoff, P., Raftery, A., and Handcock, M. (2002), "Latent Space Approaches to Social Network Analysis", Journal of the American Statistical Association, 97, 1090--1098.
#'
#' @name lvm4net-package
#' @aliases lvm4net
#' @import MASS
#' @import ergm
#' @import network
#' @importFrom stats as.dist cmdscale dist glm quantile rbinom rnorm sd weighted.mean
#' @importFrom utils combn
#' @importFrom graphics abline boxplot legend lines matlines matplot matpoints mtext par plot points polygon text
#' @importFrom grDevices rgb
#' @importFrom igraph layout.fruchterman.reingold graph.adjacency
#' @importFrom ellipse ellipse
#' @importFrom mvtnorm dmvnorm
#' @importFrom glmmML ghq
#' @importFrom corpcor wt.var 
#' @docType package
NULL