# `lvm4net`: Latent Variable Models for Networks <img src="man/figures/hex_lvm4net.png" align="left" style="padding-right:20px;padding-top:10px;"/>

`lvm4net` provides a range of tools for latent variable models for network data. Most of the models are implemented using a fast variational inference approach. 

**Latent space models for one-mode binary networks**: the function `lsm` implements the latent space model (LSM) introduced by  Hoff et al. (2002) using a variational inference and squared Euclidian distance; the function `lsjm` implements latent space joint model (LSJM) for multiplex networks introduced by Gollini and Murphy (2016). These models assume that each node of a network has a latent position in a latent space: the closer two nodes are in the latent space, the more likely they are connected.

**Latent variable models for binary bipartite networks**: the function `lca` implements the latent class analysis (LCA) to find groups in the sender nodes (with the condition of independence within the groups); the function `lta` implements the latent trait analysis (LTA) to model the dependence in the receiver nodes by using a continuous latent variable; the function `mlta` implements the mixture of latent trait analyzers (MLTA) introduced by Gollini and Murphy (2014) and Gollini (in press) to identify groups assuming the existence of a latent trait describing the dependence structure between receiver nodes within groups of sender nodes and therefore capturing the heterogeneity of sender nodes' behaviour within groups. `lta` and `mlta` use variational inference.

</br>

### References 

- Gollini, I. (in press) "A mixture model approach for clustering bipartite networks", Challenges in Social Network Research Volume in the *Lecture Notes in Social Networks* (LNSN - Series of Springer). Preprint: [arXiv:1905.02659](http://arxiv.org/abs/1905.02659).

- Gollini, I., and Murphy, T. B. (2014), ["Mixture of Latent Trait Analyzers for Model-Based Clustering of Categorical Data"](http://link.springer.com/article/10.1007/s11222-013-9389-1), *Statistics and Computing*, 24(4), 569-588, [arXiv:1301.2167](http://arxiv.org/abs/1301.2167).

- Gollini, I., and Murphy, T. B. (2016), ["Joint Modelling of Multiple Network Views"](http://www.tandfonline.com/doi/abs/10.1080/10618600.2014.978006#.VJGM9GSsW0c), *Journal of Computational and Graphical Statistics*, [arXiv:1301.3759](http://arxiv.org/abs/1301.3759).

- Hoff, P., Raftery, A., and Handcock, M. (2002), ["Latent Space Approaches to Social Network Analysis"](http://www.tandfonline.com/doi/abs/10.1198/016214502388618906#.VSglXRPF-e4), *Journal of the American Statistical Association*, 97, 1090--1098.
