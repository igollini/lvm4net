f_mlta <- function(S, counts, G, D, tol, maxiter, pdGH)
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
  
  # Initialize Variational Approximation
  
  xi <- array(20, c(Ns, M, G))
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  w <- array(rnorm(M * D * G), c(D, M, G))
  b <- matrix(rnorm(M, G), G, M)
  
  C <- array(0, c(D, D, Ns, G))
  mu <- array(0, c(Ns, D, G))
  
  lxi <- matrix(0, G, Ns)
  
  # Iterative process
  
  v <- matrix(0, Ns, G)
  ll <- -Inf
  diff <- 1
  iter <- 0
  
  xin <- xi
  
  while (diff > tol & iter < maxiter)
  {
    iter <- iter + 1
    ll_old <- ll
    
    if (D == 1) {
      for (g in 1:G)
      {
        w[, , g] <- as.matrix(t(w[, , g]))
        
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C[, , , g] <-
          1 / (1 - 2 * rowSums(sweep(lambda_xi[, , g], MARGIN = 2, w[, , g] ^ 2, `*`)))
        mu[, , g] <-
          C[, , , g] * rowSums(sweep(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`),
            MARGIN = 2,
            w[, , g],
            `*`
          ))
        
        mu[, , g] <- matrix(mu[, , g], Ns, D)
        YY <- matrix(C[, , , g] + mu[, , g] ^ 2, ncol = 1)
        
        xi[, , g] <-
          YY %*% (w[, , g] ^ 2) + mu[, , g] %*% t(2 * b[g, ] * w[, , g]) +
          matrix(b[g, ] ^ 2,
            nrow = Ns,
            ncol = M,
            byrow = TRUE)
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <- z[, g] * 2 * cbind(c(YY), mu[, , g], mu[, , g], 1)
        
        den <- t(YYh) %*% lambda_xi[, , g]
        num <- rbind(c(mu[, , g]), 1) %*% (z[, g] * (S - 0.5))
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(solve(matrix(x[1:4], 2, 2)), x[5:6]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(t(w[, , g]))
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(x¦z))
        
        lxi[g, ] <-
          0.5 * log(C[, , , g]) + c(mu[, , g]) ^ 2 / (2 * C[, , , g]) +
          rowSums(
            log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * xi[, , g] ^
              2 +
              sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) +
              sweep(S - 0.5, MARGIN = 2, b[g, ], `*`)
          )
        
      } # end for (g in 1:G)
      
      # E-step
      
      v <- t(eta * exp(lxi))
      v[is.nan(v)] <- 0
      vsum <- rowSums(v)
      z <- counts * v / vsum
      ll <- sum(counts * log(vsum))
      
      # M-step
      
      eta <- colSums(z) / N
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end if D==1
    
    if (D == 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C2 <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w[, , g]), t(w[, , g]))))
        C[, , , g] <- array(C2, c(D, D, Ns))
        
        mu[, , g] <-
          t(apply(rbind(C2, w[, , g] %*% t(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:4], nrow = D) %*% (x[-(1:4)])))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY <- C2 + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY, 2, function(x)
              rowSums((
                t(w[, , g]) %*% matrix(x, ncol = D)
              ) * t(w[, , g]))) + (2 * b[g, ] * t(w[, , g])) %*% t(mu[, , g]) + matrix(
                b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE
              )
          )
        
        # YY<-array(YY,c(D,D,Ns))
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <-
          z[, g] * 2 * cbind(YY[1, ], YY[2, ], mu[, 1, g], YY[3, ], YY[4, ], mu[, 2, g], mu[, 1, g], mu[, 2, g], 1)
        
        den <- t(YYh) %*% lambda_xi[, , g]
        num <- rbind(t(mu[, , g]), 1) %*% {
          z[, g] * {
            S - 0.5
          }
        }
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(solve(matrix(x[1:9], 3, 3)), x[10:12]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(w[, , g])
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(x¦z))
        
        detC <- C2[1, ] * C2[4, ] - C2[3, ] * C2[2, ]
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C2[4, ] / detC, -C2[2, ] / detC, -C2[3, ] /
              detC, C2[1, ] / detC, t(mu[, , g])), 2, function(x)
                t((x[-{
                  1:4
                }])) %*% matrix(x[1:4], nrow = D) %*% (x[-{
                  1:4
                }])) + rowSums(
                  log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * {
                    xi[, , g] ^ 2
                  } + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + sweep(S - 0.5, MARGIN =
                      2, b[g, ], `*`)
                )
        
      } # end for (g in 1:G)
      
      # E-step
      
      v <- t(eta * exp(lxi))
      v[is.nan(v)] <- 0
      vsum <- rowSums(v)
      z <- counts * v / vsum
      ll <- sum(counts * log(vsum))
      
      # M-step
      
      eta <- colSums(z) / N
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end if D=2
    
    if (D > 2) {
      for (g in 1:G)
      {
        # # STEP 1: Computing the Latent Posterior Statistics
        
        C2 <-
          apply(lambda_xi[, , g], 1, function(x)
            solve(diag(D) - 2 * crossprod(x * t(w[, , g]), t(w[, , g]))))
        C[, , , g] <- array(C2, c(D, D, Ns))
        
        mu[, , g] <-
          t(apply(rbind(C2, w[, , g] %*% t(
            S - 0.5 + 2 * sweep(lambda_xi[, , g], MARGIN = 2, b[g, ], `*`)
          )), 2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% (
              x[-(1:(D * D))])))
        
        # # STEP 2: Optimising the Variational Parameters (xi)
        
        YY <- C2 + apply(mu[, , g], 1, tcrossprod)
        
        xi[, , g] <-
          t(
            apply(YY, 2, function(x)
              rowSums((t(w[, , g]) %*% matrix(x, ncol = D)) * t(w[, , g]))) + (2 * b[g, ] * t(w[, , g])) %*% t(mu[, , g]) +
              matrix(b[g, ] ^ 2,
                nrow = M,
                ncol = Ns,
                byrow = FALSE)
          )
        
        xi[, , g] <- sqrt(xi[, , g])
        sigma_xi[, , g] <- 1 / (1 + exp(-xi[, , g]))
        lambda_xi[, , g] <- (0.5 - sigma_xi[, , g]) / (2 * xi[, , g])
        
        # # STEP 3: Optimising the Model Parameters (w and b)
        
        YYh <-
          sweep(buildYYh(YY, t(mu[, , g]), D, Ns), MARGIN = 2, 2 * z[, g], `*`)
        den <- (YYh) %*% lambda_xi[, , g]
        num <- rbind(t(mu[, , g]), 1) %*% (z[, g] * (S - 0.5))
        wh <-
          apply(rbind(den, num), 2, function(x)
            - crossprod(solve(matrix(x[1:((D + 1) * (D + 1))], D + 1, D + 1)), x[-(1:((D +
                1) * (D + 1)))]))
        
        w[, , g] <- wh[1:D, ]
        w[, , g] <- as.matrix(w[, , g])
        b[g, ] <- wh[D + 1, ]
        
        # Approximation of log(p(x¦z))
        
        detC <- apply(C2, 2, function(x)
          det(matrix(x, D, D)))
        
        lxi[g, ] <-
          0.5 * log(detC) + 0.5 * apply(rbind(C2, t(mu[, , g])), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*% (x[-(1:(D * D))])) + 
          rowSums(
              log(sigma_xi[, , g]) - 0.5 * xi[, , g] - lambda_xi[, , g] * (xi[, , g] ^
                  2) + sweep(lambda_xi[, , g], MARGIN = 2, b[g, ] ^ 2, `*`) + sweep(S - 0.5, MARGIN =
                      2, b[g, ], `*`)
            )
      } # end for (g in 1:G)
      
      # E-step
      
      v <- t(eta * exp(lxi))
      v[is.nan(v)] <- 0
      vsum <- rowSums(v)
      z <- counts * v / vsum
      ll <- sum(counts * log(vsum))
      
      # M-step
      
      eta <- colSums(z) / N
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end if D>2
    
  } # end while(diff>tol)
  
  # Correction to the log-likelihood
  
  # Gauss-Hermite Quadrature
  
  npoints <- round(pdGH ^ (1 / D))
  ny <- npoints ^ D
  GaHer <- glmmML::ghq(npoints, FALSE)
  Ygh <- expand.grid(rep(list(GaHer$zeros), D))
  Ygh <- as.matrix(Ygh)
  Wgh <- apply(as.matrix(expand.grid(rep(
    list(GaHer$weights), D
  ))), 1, prod) *
    apply(exp(Ygh ^ 2), 1, prod)
  
  Hy <- apply(Ygh, 1, mvtnorm::dmvnorm)
  Beta <- Hy * Wgh / sum(Hy * Wgh)
  
  fxy <- array(0, c(Ns, ny, G))
  
  for (g in 1:G)
  {
    if (D == 1) {
      Agh <- t(tcrossprod(w[, , g], Ygh) + b[g, ])
    } else{
      Agh <- t(tcrossprod(t(w[, , g]), Ygh) + b[g, ])
    }
    pgh <- 1 / (1 + exp(-Agh))
    
    fxy[, , g] <-
      exp(tcrossprod(S, log(pgh)) + tcrossprod(1 - S, log(1 - pgh)))
    fxy[is.nan(fxy[, , g])] <- 0
  }
  
  eta <- as.vector(eta)
  
  LLva <- ll
  BICva <-
    -2 * LLva + (G * (M * (D + 1) - D * (D - 1) / 2) + G - 1) * log(N)
  
  llGH1 <- apply(rep(Beta, each = Ns) * fxy, c(1, 3), sum)
  LL <- sum(counts * log(colSums(eta * t(llGH1))))
  
  BIC <-
    -2 * LL + (G * (M * (D + 1) - D * (D - 1) / 2) + G - 1) * log(N)
  
  expected <- colSums(eta * t(llGH1)) * N
  
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- rownames(b, do.NULL = FALSE, prefix = "Group ")
  
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  dimnames(w)[[3]] <- paste("Group ", seq(1, G) , sep = '')
  
  eta <- matrix(eta, byrow = T, G)
  rownames(eta) <- rownames(eta, do.NULL = FALSE, prefix = "Group ")
  colnames(eta) <- "eta"
  names(LLva) <- c("Log-Likelihood (variational approximation):")
  names(BICva) <- c("BIC (variational approximation):")
  names(LL) <- c("Log-Likelihood (G-H Quadrature correction):")
  names(BIC) <- c("BIC (G-H Quadrature correction):")
  
  out <- list(
    b = b,
    w = w,
    eta = eta,
    mu = mu,
    C = C,
    z = z,
    LL = LL,
    BIC = BIC,
    LLva = LLva,
    BICva = BICva,
    expected = expected
  )
  
  out
}
