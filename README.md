# `lvm4net`: Latent Variable Models for Networks <img src="man/figures/hex_lvm4net.png" align="left" style="padding-right:20px;padding-top:10px;"/>

`lvm4net` provides a range of tools for latent variable models for network data. Most of the models are implemented using a fast variational inference approach. Latent space models for binary networks: the function `lsm` implements the latent space model (LSM) introduced by  Hoff et al. (2002) using a variational inference and squared Euclidian distance; the function `lsjm` implements latent space joint model (LSJM) for multiplex networks introduced by Gollini and Murphy (2016). These models assume that each node of a network has a latent position in a latent space: the closer two nodes are in the latent space, the more likely they are connected.

*Functions for binary bipartite networks will be added soon.*

</br>

### References 
- Gollini, I., and Murphy, T. B. (2016), ["Joint Modelling of Multiple Network Views"](http://www.tandfonline.com/doi/abs/10.1080/10618600.2014.978006#.VJGM9GSsW0c), *Journal of Computational and Graphical Statistics*, [arXiv:1301.3759](http://arxiv.org/abs/1301.3759).

- Hoff, P., Raftery, A., and Handcock, M. (2002), ["Latent Space Approaches to Social Network Analysis"](http://www.tandfonline.com/doi/abs/10.1198/016214502388618906#.VSglXRPF-e4), *Journal of the American Statistical Association*, 97, 1090--1098.
