f_mlta_wfix <- function(S, counts, G, D, tol, maxiter, pdGH)
{
  Ns <- nrow(S)
  N <- sum(counts)
  M <- ncol(S) 
  rownames(S) <- NULL
  
  # Initialize EM Algorithm
  
  eta <- rep(1 / G, G)
  p <- matrix(1 / M, G, M)
  
  lab <- sample(1:G,
    prob = eta,
    size = Ns,
    replace = TRUE)
  z <- counts * cls(lab)
  
  # Preliminary M-step
  
  eta <- colSums(z) / N
  
  ###
  
  xi <- array(20, c(Ns, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- matrix(rnorm(M * D), D, M)
  b <- matrix(rnorm(M, G), G, M)
  wh <- matrix(0, D + G, M)
  
  C <- array(0, c(D * D, Ns, G))
  
  mu <- array(0, c(Ns, D, G))
  
  YY <- array(0, c(D * D, Ns, G))
  
  gam <- matrix(0, D + G, M)
  K <- array(diag(D + G), c(D + G, D + G, M))
  
  lxi <- matrix(0, G, Ns)
  
  # Iterative process
  
  v <- matrix(0, Ns, G)
  ll <- -Inf
  diff <- 1
  iter <- 0
  
  tol <- 0.1 ^ 2	#0.05
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w ^ 2, `*`)))
        mu[, , g] <-
          C[, , g] * rowSums(sweep(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2, w, `*`))
        
        YY[, , g] <- matrix(C[, , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <- YY[, , g] %*% w^2 + mu[, , g] %*% (2 * b[g, ] * w) + 
          matrix(b[g, ]^2, nrow = Ns, ncol = M, byrow = TRUE)
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[1, ] <- colSums(crossprod(z * mu[, 1, ], S - 0.5))
      gam[2:(1 + G), ] <- crossprod(z, S - 0.5)
      
      K[1, 1, ] <- apply(aperm(array(z * YY[1, , ], c(Ns, G, M)), c(1, 3, 2)) * lambda_xi, 2, sum)
      
      K[2:(G + 1), 1, ] <- t(apply(aperm(array(
          z * mu[, 1, ], c(Ns, G, M)
        ), c(1, 3, 2)) * lambda_xi, c(2, 3), sum))
      K[1, 2:(G + 1), ] <- K[2:(G + 1), 1, ]
      
      for (m in 1:M)
      {
        for (g in 1:G)
          K[1 + g, 1 + g, m] <- sum(z[, g] * lambda_xi[, m, g])
        
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      w <- as.matrix(t(w))
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        lxi[g, ] <-
          0.5 * log(C[1, , g]) + mu[, 1, g] ^ 2 / (2 * C[1, , g]) + rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^ 2 + 
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
              sweep(S - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
    }
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], 
            tcrossprod(w, S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`))), 
            2, function(x) matrix(x[1:4], nrow = D) %*% x[-(1:4)]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums(crossprod(w, matrix(x, ncol = D)) * t(w))) + 
              tcrossprod(2 * b[g, ] * t(w), mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE
              )
          )
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # wh<-rbind(w,b)
      gam[3:(2 + G), ] <- crossprod(z, S - 0.5)
      
      aa <- aperm(array(z, c(Ns, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(Ns, G, D * D)), c(3, 1, 2)) * YY
      
      K[3:(2 + G), 3:(2 + G), ] <-
        apply(apply(aperm(array(z, c(
          Ns, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[2 + g, 1:2, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:2, 2 + g, ] <- K[2 + g, (1:2), ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:2, 1:2, ] <- kk
      
      for (m in 1:M)
      {
        gam[1:D, m] <- apply(aa * (S[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- C[1, , g] * C[4, , g] - C[3, , g] * C[2, , g]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[4, , g] / detC, -C[2, , g] / detC, -C[3, , g] /
              detC, C[1, , g] / detC, t(mu[, , g])), 2, function(x)
                t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) + 
          rowSums(
                  log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g]^2 + 
                    sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                    sweep(S - 0.5, MARGIN = 2, b[g, ], `*`))
      } # end for (g in 1:G)
      
    }# end D=2
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , g] <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
        
        mu[, , g] <-
          t(apply(rbind(C[, , g], w %*% t(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% x[-(1:(D * D))]))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY[, , g] <- C[, , g] + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY[, , g], 2, function(x)
              rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) + 
              (2 * b[g, ] * t(w)) %*% t(mu[, , g]) + 
              matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE
              )
          )
        
      }
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      gam[(D + 1):(D + G), ] <- t(z) %*% (S - 0.5)
      
      aa <- aperm(array(z, c(Ns, G, D)), c(1, 3, 2)) * mu
      bb <- aperm(array(z, c(Ns, G, D * D)), c(3, 1, 2)) * YY
      
      K[(D + 1):(D + G), (D + 1):(D + G), ] <-
        apply(apply(aperm(array(z, c(
          Ns, G, M
        )), c(1, 3, 2)) * lambda_xi, c(2, 3), sum), 1, diag)
      
      kk <- 0
      for (g in 1:G) {
        K[D + g, 1:D, ] <- crossprod(aa[, , g], lambda_xi[, , g])
        K[1:D, D + g, ] <- K[D + g, 1:D, ]
        kk <- kk + bb[, , g] %*% lambda_xi[, , g]
      }
      
      K[1:D, 1:D, ] <- kk
      
      for (m in 1:M)	{
        gam[1:D, m] <- apply(aa * (S[, m] - 0.5), 2, sum)
        wh[, m] <- -solve(2 * K[, , m]) %*% gam[, m]
      }
      
      w <- wh[1:D, ]
      
      for (g in 1:G)
      {
        b[g, ] <- wh[D + g, ]
        
        # Approximation of log(p(x|z))
        
        detC <- apply(C[, , g], 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C[, , g], t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(
              1:(D * D))])) + rowSums(
              log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] *
                xi[, , g]^2 + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + 
                  sweep(S - 0.5, MARGIN = 2, b[g, ], `*`)
            )
      } # end for (g in 1:G)
      
    } # end if D > 2
    
    # E-step
    
    v <- t(eta * exp(lxi))
    v[is.nan(v)] <- 0
    vsum <- apply(v, 1, sum)
    z <- counts * v / vsum
    ll <- sum(counts * log(vsum))
    
    # M-step
    
    eta <- apply(z, 2, sum) / N
    
    # Stopping Criteria
    
    diff <- sum(abs(ll - ll_old))
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood 
  
  # Gauss-Hermite Quadrature
  
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
  
  fxy <- array(0, c(Ns, ny, G))
  
  for (g in 1:G)
  {
    Agh <- t(tcrossprod(t(w), Ygh) + b[g, ])
    pgh <- 1 / (1 + exp(-Agh))
    fxy[, , g] <-
      exp(tcrossprod(S, log(pgh)) + tcrossprod(1 - S, log(1 - pgh)))
  }
  
  eta <- as.vector(eta)
  
  LLva <- ll
  BICva <- -2 * LLva + {
    G * M + M * D - D * {
      D - 1
    } / 2 + G - 1
  } * log(N)
  
  llGH1 <- apply(rep(Beta, each = Ns) * fxy, c(1, 3), sum)
  LL <- sum(counts * log(colSums(eta * t(llGH1))))
  
  BIC <- -2 * LL + {
    G * M + M * D - D * {
      D - 1
    } / 2 + G - 1
  } * log(N)
  
  expected <- colSums(eta * t(llGH1)) * N
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  
  eta <- matrix(eta, byrow = T, G)
  rownames(eta) <- rownames(eta, do.NULL = FALSE, prefix = "Group ")
  colnames(eta) <- "eta"
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")

  out1 <-
    list(
      b = b,
      w = w,
      eta = eta,
      LL = LL,
      BIC = BIC,
      LLva = LLva,
      BICva = BICva,
      expected = expected,
      mu = mu,
      C = C,
      z = z
    )
  
  return(out1)
  
}
