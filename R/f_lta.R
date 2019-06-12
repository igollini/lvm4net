f_lta <- function(S, counts, D, tol, maxiter, pdGH) {
  Ns <- nrow(S)
  N <- sum(counts)
  M <- ncol(S) 
  rownames(S) <- NULL
  
  # Initialization
  
  xi <- matrix(20, Ns, M)
  w <- matrix(rnorm(M * D), D, M)
  b <- rnorm(M)
  
  sigma_xi <- 1 / (1 + exp(-xi))
  lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
  
  diff <- 1
  iter <- 0L
  ll <- -Inf
  
  while (diff > tol & iter < maxiter) {
    iter <- iter + 1
    ll_old <- ll
    
    if (D == 1) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <-
        1 / (1 - 2 * rowSums(sweep(lambda_xi, MARGIN = 2, w ^ 2, `*`)))
      mu <-
        C * rowSums(sweep(
          S - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`),
          MARGIN = 2,
          w,
          `*`
        ))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- matrix(C + mu ^ 2, ncol = 1)
      mu <- matrix(mu, Ns, D)
      
      xi <-
        YY %*% w ^ 2 + mu %*% (2 * b * w) + matrix(b ^ 2,
          nrow = Ns,
          ncol = M,
          byrow = TRUE)
      
      xi <- sqrt(xi)
      
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <- counts * 2 * cbind(c(YY), mu, mu, 1)
      
      den <- t(YYh) %*% lambda_xi
      num <- rbind(c(mu), 1) %*% (counts * (S - 0.5))
      
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(solve(matrix(x[1:4], 2, 2)), x[5:6]))
      
      w <- as.matrix(t(wh[1:D, ]))
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      lxi <- counts * (0.5 * log(C) + c(mu) ^ 2 / (2 * C) +
          rowSums(
            log(sigma_xi) - 0.5 * xi - lambda_xi * (xi ^ 2) +
              sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) +
              sweep(S - 0.5, MARGIN = 2, b, `*`)
          ))
      
      C <- array(C, c(D, D, Ns))
      
      ll <- sum(lxi)
      bw <- t(rbind(b, w))
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
    } # end if D=1
    
    if (D == 2) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <- apply(lambda_xi, 1, function(x)
        solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
      
      mu <- apply(rbind(C, w %*% t(
        S - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)
      )),
        2,
        function(x)
          matrix(x[1:4], nrow = D) %*% (x[-(1:4)]))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- C + apply(mu, 2, tcrossprod)
      
      xi <-
        apply(YY, 2, function(x)
          rowSums((t(w) %*% matrix(x, ncol = D)) * t(w))) +
        (2 * b * t(w)) %*% mu +
        matrix(b ^ 2,
          nrow = M,
          ncol = Ns,
          byrow = FALSE)
      
      xi <- t(sqrt(xi))
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <-
        counts * 2 * cbind(YY[1, ], YY[2, ], mu[1, ], YY[3, ], YY[4, ], mu[2, ], mu[1, ], mu[2, ], 1)
      
      den <- t(YYh) %*% lambda_xi
      num <- rbind(mu, 1) %*% (counts * (S - 0.5))
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(solve(matrix(x[1:9], 3, 3)), x[10:12]))
      
      w <- wh[1:D, ]
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      detC <- C[1, ] * C[4, ] - C[3, ] * C[2, ]
      
      lxi <- counts * (
        0.5 * log(detC) +
          0.5 * apply(rbind(C[4, ] / detC,-C[2, ] / detC,-C[3, ] / detC, C[1, ] / detC, mu), 2,
            function(x)
              t((x[-(1:4)])) %*% matrix(x[1:4], nrow = D) %*% (x[-(1:4)])) +
          rowSums(
            log(sigma_xi) - 0.5 * xi - lambda_xi * (xi ^ 2) +
              sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + sweep(S - 0.5, MARGIN = 2, b, `*`)
          )
      )
      
      C <- array(C, c(D, D, Ns))
      mu <- t(mu)
      
      ll <- sum(lxi)
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    }# end if D=2
    
    if (D > 2) {
      # # STEP 1: Computing the Latent Posterior Statistics
      
      C <-
        apply(lambda_xi, 1, function(x)
          solve(diag(D) - 2 * crossprod(x * t(w), t(w))))
      
      mu <-
        apply(rbind(C, w %*% t(
          S - 0.5 + 2 * sweep(lambda_xi, MARGIN = 2, b, `*`)
        )) ,
          2, function(x)
            matrix(x[1:(D * D)], nrow = D) %*% (x[-(1:(D * D))]))
      
      # # STEP 2: Optimising the Variational Parameters (xi)
      
      YY <- C + apply(mu, 2, tcrossprod)
      
      xi <-
        t(
          apply(YY, 2, function(x)
            rowSums((
              t(w) %*% matrix(x, ncol = D)
            ) * t(w))) +
            (2 * b * t(w)) %*% mu + matrix(
              b ^ 2,
              nrow = M,
              ncol = Ns,
              byrow = FALSE
            )
        )
      
      xi <- sqrt(xi)
      sigma_xi <- 1 / (1 + exp(-xi))
      lambda_xi <- (0.5 - sigma_xi) / (2 * xi)
      
      # # STEP 3: Optimising the Model Parameters (w and b)
      
      YYh <- sweep(buildYYh(YY, mu, D, Ns), MARGIN = 2, 2 * counts, `*`)
      
      den <- YYh %*% lambda_xi
      
      num <- rbind(mu, 1) %*% (counts * (S - 0.5))
      
      wh <-
        apply(rbind(den, num), 2, function(x)
          - crossprod(solve(matrix(x[1:((D + 1) * (D + 1))], D + 1, D + 1)), x[-(1:((D + 1) * (D + 1)))]))
      
      w <- wh[1:D, ]
      b <- wh[D + 1, ]
      
      # Approximation of Likelihood
      
      detC <- apply(C, 2, function(x)
        det(matrix(x, D, D)))
      
      lxi <- counts * (
        0.5 * log(detC) +
          0.5 * apply(rbind(C, mu), 2, function(x)
            t((x[-(1:(D * D))])) %*% solve(matrix(x[1:(D * D)], nrow = D)) %*%
              (x[-(1:(D * D))])) +
          rowSums(
            log(sigma_xi) - 0.5 * xi - lambda_xi * xi ^ 2 +
              sweep(lambda_xi, MARGIN = 2, b ^ 2, `*`) + sweep(S - 0.5, MARGIN = 2, b, `*`)
          )
      )
      
      C <- array(C, c(D, D, Ns))
      mu <- t(mu)
      
      ll <- sum(lxi)
      
      # Stopping Criteria
      
      diff <- sum(abs(ll - ll_old))
      
    } # end if D>2
    
  } # End of while statement
  
  # correction to the log-likelihood
  
  ghLL <- ghLLlta(S, counts, w, b, D, pdGH)
  expected <- ghLL$expected
  LL <- ghLL$llGH
  
  LLva <- ll
  
  BICva <- -2 * LLva + (M * (D + 1) - D * (D - 1) / 2) * log(N)
  BIC <- -2 * LL + (M * (D + 1) - D * (D - 1) / 2) * log(N)
  
  b <- t(as.matrix(b))
  rownames(b) <- NULL
  colnames(b) <- colnames(b, do.NULL = FALSE, prefix = "Item ")
  rownames(b) <- ""
  rownames(w) <- NULL
  colnames(w) <- colnames(w, do.NULL = FALSE, prefix = "Item ")
  rownames(w) <- rownames(w, do.NULL = FALSE, prefix = "Dim ")
  names(LLva) <- "Log-Likelihood (variational approximation):"
  names(BICva) <- "BIC (variational approximation):"
  names(LL) <- "Log-Likelihood (G-H Quadrature correction):"
  names(BIC) <- "BIC (G-H Quadrature correction):"
  
  list(
    b = b,
    w = w,
    LLva = LLva,
    BICva = BICva,
    LL = LL,
    BIC = BIC,
    mu = mu,
    C = C,
    expected = expected
  )
}
