
# verossimilhança da inflação de zeros
# w e uma matriz da forma (1, x1, ..., x_h)
# tem que entrar com gama = (gama0, gama1, ..., gama_h)

veroZero <- function(w, gamma) {
  
  W0 <- as.matrix(w[w$tempo == 0, 2:ncol(w)])
  W <- as.matrix(w[w$tempo > 0, 2:ncol(w)])
  aux1 <- sum((W0%*%gamma)) - sum(log( 1 + exp(W0%*%gamma))) - sum(log(1 + exp(W%*%gamma)))
  return(-aux1)
}

# verossimilhança do GTDL com fragilidade Gama
# matriz x = (censura, tempo, x1)
# par e um vetor dado por (alpha, log(lambda), log(theta), beta)

veroGTDL <- function(x, par) {
  
  alpha <- par[1]
  lambda <- exp(par[2])
  theta <- exp(par[3])
  beta <- par[4]
  
  cens <- x[x$tempo > 0,1]
  tempo <- x[x$tempo > 0,2]
  X <- as.matrix(x[x$tempo > 0,3])
  
  aux2 <- log(lambda)*sum(cens) +
    sum( cens*(alpha*tempo + X%*%beta)) -
    sum( cens*log( 1 + exp(alpha*tempo + X%*%beta) ) ) - 
    sum( (cens + (1/theta))*log( 1 + ((theta*lambda/alpha)*log( (1 + exp(alpha*tempo + X%*%beta))/(1 + exp(X%*%beta)) )) ) )
  return(-aux2)
}




##
# função de sobrevivência do modelo com regressão apenas na inflação de zero
sobrevGTDLv1 <- function(t, x, par) {
  
  gamma <- par[1:2]
  alpha <- par[3]
  lambda <- par[4]
  theta <- par[5]
  beta <- par[6]
  
  tempo <- t
  X <- as.matrix(x)
  X0 <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
  
  st <- (1/(1 + exp(X0%*%gamma) )) * (1 + ( (lambda*theta/alpha) * log( (1 + exp(alpha*tempo + X%*%beta) )/(1 + exp(X%*%beta))) ) )^(-1/theta)
  return(st)
}


# função que calcula a proporção de zeros inflacionados dados os níveis da covar
p0f <- function(w, gamma) {
  
  w <- as.matrix(w)
  W <- matrix(data = c(rep(1,nrow(w)),w), nrow = nrow(w), ncol = 2)
  
  p <- exp(W%*%gamma)/(1+exp(W%*%gamma))
  return(p)
}


# função que calcula a fração de cura dados os níveis da covar
p1f <- function(x, par) {
  gamma <- par[1:2]
  alpha <- par[3]
  lambda <- par[4]
  theta <- par[5]
  beta1 <- par[6]
  
  X <- as.matrix(x)
  W <- matrix(data = c(rep(1,nrow(X)),X), nrow = nrow(X), ncol = 2)
  
  p = (1/(1+exp(W%*%gamma))) * (1 + ((lambda*theta)/alpha)*log(1/(1+exp(X%*%beta1))) )^(-1/theta)
  return(p)
}


# função criada para usar no comando uniroot e encontrar a raiz
froot <- function(t,x,par,u) sobrevGTDLv1(t,x,par) - 1 + u



one_rep <- function(nk, params_fixed, interval_root = c(0, 1080)) {
  # unpack fixed params
  gama0  <- params_fixed$gama0
  gama1  <- params_fixed$gama1
  alpha  <- params_fixed$alpha
  lambda <- params_fixed$lambda
  theta  <- params_fixed$theta
  beta1  <- params_fixed$beta1
  
  # --- generate covariate
  X <- rbinom(nk, size = 1, prob = 0.5)
  
  # probabilities
  p0 <- p0f(X, gamma = c(gama0, gama1))
  p1 <- p1f(X, par = c(gama0, gama1, alpha, lambda, theta, beta1))
  
  # --- generate times
  t <- numeric(nk)
  delta <- numeric(nk)
  
  for (i in 1:nk) {
    u1 <- runif(1, 0, 1)
    u2 <- runif(1, min = p0[i], max = (1 - p1[i]))
    
    tf <- ifelse(u1 <= p0[i],
                 0,
                 ifelse((1 - p1[i]) <= u1,
                        Inf,
                        uniroot(
                          f = froot,
                          x = X[i],
                          u = u2,
                          par = c(gama0, gama1, alpha, lambda, theta, beta1),
                          interval = interval_root
                        )$root
                 )
    )
    
    tf <- as.numeric(tf)
    tc <- runif(1, 0, 65)
    
    t[i] <- min(tf, tc)
    delta[i] <- ifelse(t[i] == tf, 1, 0)
  }
  
  # --- estimation
  dadoSimulado <- data.frame(cens = delta, tempo = t, x1 = X)
  dado0 <- data.frame(tempo = t, w0 = rep(1, nk), w1 = X)
  
  # Fit GTDL
  emvg <- try(
    optim(par = c(-0.1, -0.5, -0.5, 1),
          fn = veroGTDL, x = dadoSimulado,
          hessian = TRUE, method = "BFGS"),
    TRUE
  )
  if (inherits(emvg, "try-error")) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 1, gtdl_solve = 0, gtdl_negvar = 0,
                        zero_fail = 0, zero_solve = 0, zero_negvar = 0)))
  }
  
  v_solver <- try(solve(emvg$hessian), TRUE)
  if (inherits(v_solver, "try-error")) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 0, gtdl_solve = 1, gtdl_negvar = 0,
                        zero_fail = 0, zero_solve = 0, zero_negvar = 0)))
  }
  if (min(diag(v_solver)) < 0) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 0, gtdl_solve = 0, gtdl_negvar = 1,
                        zero_fail = 0, zero_solve = 0, zero_negvar = 0)))
  }
  
  # Fit zero-inflation
  emv0 <- try(
    optim(par = c(1, 1),
          fn = veroZero, w = dado0,
          hessian = TRUE, method = "BFGS"),
    TRUE
  )
  if (inherits(emv0, "try-error")) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 0, gtdl_solve = 0, gtdl_negvar = 0,
                        zero_fail = 1, zero_solve = 0, zero_negvar = 0)))
  }
  
  v_solver0 <- try(solve(emv0$hessian), TRUE)
  if (inherits(v_solver0, "try-error")) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 0, gtdl_solve = 0, gtdl_negvar = 0,
                        zero_fail = 0, zero_solve = 1, zero_negvar = 0)))
  }
  if (min(diag(v_solver0)) < 0) {
    return(list(ok = FALSE,
                err = c(gtdl_fail = 0, gtdl_solve = 0, gtdl_negvar = 0,
                        zero_fail = 0, zero_solve = 0, zero_negvar = 1)))
  }
  
  # --- store estimates + SEs (same as your code)
  est_gtdl <- c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4])
  
  sdg <- sqrt(diag(v_solver))
  sdg[2] <- deltamethod(g = ~ exp(x1), mean = emvg$par[2], cov = v_solver[2,2, drop=FALSE])
  sdg[3] <- deltamethod(g = ~ exp(x1), mean = emvg$par[3], cov = v_solver[3,3, drop=FALSE])
  
  est_zero <- c(emv0$par[1], emv0$par[2])
  sdz <- sqrt(diag(v_solver0))
  
  dp_p0_1 <- deltamethod(
    g = ~ exp(x1 + x2)/(1 + exp(x1 + x2)),
    mean = emv0$par, cov = v_solver0
  )
  dp_p0_0 <- deltamethod(
    g = ~ exp(x1)/(1 + exp(x1)),
    mean = emv0$par, cov = v_solver0
  )
  
  # joint covariance block (same structure you had)
  matriz_aux <- diag(6)
  matriz_aux[1:4, 1:4] <- emvg$hessian
  matriz_aux[5:6, 5:6] <- emv0$hessian
  
  EMV <- c(emvg$par, emv0$par)
  cov_joint <- solve(matriz_aux)
  
  dp_p1_1 <- deltamethod(
    g = ~ (1 - exp(x5 + x6)/(1 + exp(x5 + x6))) *
      (1 - ((exp(x2) * exp(x3))/x1) * log(1 + exp(x4)))^(-1/exp(x3)),
    mean = EMV, cov = cov_joint
  )
  
  dp_p1_0 <- deltamethod(
    g = ~ (1 - exp(x5)/(1 + exp(x5))) *
      (1 - ((exp(x2) * exp(x3))/x1) * log(1 + exp(0)))^(-1/exp(x3)),
    mean = EMV, cov = cov_joint
  )
  
  list(
    ok = TRUE,
    est_gtdl = est_gtdl,
    se_gtdl  = sdg,
    est_zero = est_zero,
    se_zero  = sdz,
    dp_p0    = c(dp_p0_1, dp_p0_0),
    dp_p1    = c(dp_p1_1, dp_p1_0),
    err = c(gtdl_fail = 0, gtdl_solve = 0, gtdl_negvar = 0,
            zero_fail = 0, zero_solve = 0, zero_negvar = 0)
  )
}



#########


sobrevGTDL_alpha <- function(t, x, par) {
  gamma <- par[1:2]
  alpha0 <- par[3]; alpha1 <- par[4]
  lambda <- par[5]; theta <- par[6]; beta1 <- par[7]
  
  X <- as.matrix(x)
  W <- cbind(1, X)
  p0 <- as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
  
  alpha_i <- as.numeric(alpha0 + alpha1 * X)
  xb  <- as.numeric(X %*% beta1)
  eta <- as.numeric(alpha_i * t + xb)
  
  A    <- log((1 + exp(eta)) / (1 + exp(xb)))
  Spos <- (1 + (lambda*theta/alpha_i) * A)^(-1/theta)
  
  (1 - p0) * Spos
}



veroGTDL_alpha <- function(x, par) {
  # par = c(alpha0, alpha1, log(lambda), log(theta), beta1)
  
  alpha0 <- par[1]
  alpha1 <- par[2]
  lambda <- exp(par[3])
  theta  <- exp(par[4])
  beta1  <- par[5]
  
  cens  <- x[x$tempo > 0, 1]
  tempo <- x[x$tempo > 0, 2]
  X     <- as.matrix(x[x$tempo > 0, 3])
  
  alpha_i <- as.numeric(alpha0 + alpha1 * X)
  xb      <- as.numeric(X %*% beta1)
  eta     <- as.numeric(alpha_i * tempo + xb)
  
  term1 <- log(lambda) * sum(cens)
  term2 <- sum(cens * eta)
  term3 <- -sum(cens * log(1 + exp(eta)))
  
  A <- log((1 + exp(eta)) / (1 + exp(xb)))
  B <- 1 + (theta * lambda / alpha_i) * A
  
  term4 <- -sum((cens + (1/theta)) * log(B))
  
  return(-(term1 + term2 + term3 + term4))  # negative log-likelihood
}



veroGTDL_alpha_H0 <- function(x, par) {
  # par = c(alpha0, log(lambda), log(theta)), with alpha1=0, beta1=0
  
  alpha0 <- par[1]
  lambda <- exp(par[2])
  theta  <- exp(par[3])
  
  cens  <- x[x$tempo > 0, 1]
  tempo <- x[x$tempo > 0, 2]
  X     <- as.matrix(x[x$tempo > 0, 3])  # used only for dimension consistency
  
  alpha_i <- rep(alpha0, nrow(X))
  xb      <- rep(0,      nrow(X))
  eta     <- alpha_i * tempo + xb
  
  term1 <- log(lambda) * sum(cens)
  term2 <- sum(cens * eta)
  term3 <- -sum(cens * log(1 + exp(eta)))
  
  A <- log((1 + exp(eta)) / (1 + exp(xb)))  # xb=0
  B <- 1 + (theta * lambda / alpha_i) * A
  
  term4 <- -sum((cens + (1/theta)) * log(B))
  
  return(-(term1 + term2 + term3 + term4))
}

veroZero_H0 <- function(w, gamma0) {
  W0 <- as.matrix(w[w$tempo == 0, 2])  # intercept column only
  W  <- as.matrix(w[w$tempo > 0, 2])
  
  eta0 <- as.numeric(W0 * gamma0)
  eta  <- as.numeric(W  * gamma0)
  
  ll <- sum(eta0) - sum(log(1 + exp(eta0))) - sum(log(1 + exp(eta)))
  return(-ll)
}



one_rep_LR_type1 <- function(nk, params_fixed, interval_root = c(0, 1080), alpha_sig = 0.05) {
  
  # fixed params
  gamma0 <- params_fixed$gamma0
  gamma1 <- params_fixed$gamma1
  alpha0 <- params_fixed$alpha0
  alpha1 <- params_fixed$alpha1
  lambda <- params_fixed$lambda
  theta  <- params_fixed$theta
  beta1  <- params_fixed$beta1
  
  # covariate (group indicator)
  X <- rbinom(nk, size = 1, prob = 0.5)
  
  # --- helpers from your code
  p0f <- function(w, gamma) {
    w <- as.matrix(w)
    W <- cbind(1, w)
    as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
  }
  
  p1f_alpha <- function(x, par) {
    # par = c(gamma0,gamma1,alpha0,alpha1,lambda,theta,beta1)
    gamma <- par[1:2]
    alpha0 <- par[3]
    alpha1 <- par[4]
    lambda <- par[5]
    theta  <- par[6]
    beta1  <- par[7]
    
    X <- as.matrix(x)
    W <- cbind(1, X)
    p0 <- as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
    
    alpha_i <- as.numeric(alpha0 + alpha1 * X)
    xb      <- as.numeric(X %*% beta1)
    
    p_star <- (1 - (lambda*theta/alpha_i) * log(1 + exp(xb)))^(-1/theta)
    (1 - p0) * p_star
  }
  
  # survival for generator
  sobrevGTDL_alpha <- function(t, x, par) {
    gamma <- par[1:2]
    alpha0 <- par[3]; alpha1 <- par[4]
    lambda <- par[5]; theta <- par[6]; beta1 <- par[7]
    
    X <- as.matrix(x)
    W <- cbind(1, X)
    p0 <- as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
    
    alpha_i <- as.numeric(alpha0 + alpha1 * X)
    xb  <- as.numeric(X %*% beta1)
    eta <- as.numeric(alpha_i * t + xb)
    
    A    <- log((1 + exp(eta)) / (1 + exp(xb)))
    Spos <- (1 + (lambda*theta/alpha_i) * A)^(-1/theta)
    
    (1 - p0) * Spos
  }
  
  froot_alpha <- function(t, x, par, u) sobrevGTDL_alpha(t, x, par) - 1 + u
  
  # probabilities
  p0 <- p0f(X, gamma = c(gamma0, gamma1))
  p1 <- p1f_alpha(X, par = c(gamma0, gamma1, alpha0, alpha1, lambda, theta, beta1))
  
  # --- generate times
  t <- numeric(nk)
  delta <- numeric(nk)
  
  for (i in 1:nk) {
    u1 <- runif(1)
    u2 <- runif(1, min = p0[i], max = (1 - p1[i]))
    
    tf <- if (u1 <= p0[i]) {
      0
    } else if ((1 - p1[i]) <= u1) {
      Inf
    } else {
      uniroot(
        f = froot_alpha,
        x = X[i],
        u = u2,
        par = c(gamma0, gamma1, alpha0, alpha1, lambda, theta, beta1),
        interval = interval_root
      )$root
    }
    
    tc <- runif(1, 0, 65)
    
    t[i] <- min(tf, tc)
    delta[i] <- as.integer(t[i] == tf)
  }
  
  dadoSimulado <- data.frame(cens = delta, tempo = t, x1 = X)
  dado0        <- data.frame(tempo = t, w0 = rep(1, nk), w1 = X)
  
  # ---------- UNRESTRICTED FIT ----------
  fit_gtdl <- try(
    optim(
      par = c(-0.1, 0.0, log(0.5), log(0.4), 0.0),   # (alpha0, alpha1, loglambda, logtheta, beta1)
      fn = veroGTDL_alpha, x = dadoSimulado,
      method = "BFGS"
    ),
    TRUE
  )
  if (inherits(fit_gtdl, "try-error")) return(list(ok = FALSE))
  
  fit_zero <- try(
    optim(
      par = c(-1, 0),
      fn = veroZero, w = dado0,
      method = "BFGS"
    ),
    TRUE
  )
  if (inherits(fit_zero, "try-error")) return(list(ok = FALSE))
  
  ll_full <- -( fit_gtdl$value + fit_zero$value )
  
  # ---------- RESTRICTED FIT under H0: alpha1=beta1=gamma1=0 ----------
  fit_gtdl0 <- try(
    optim(
      par = c(-0.1, log(0.5), log(0.4)),             # (alpha0, loglambda, logtheta)
      fn = veroGTDL_alpha_H0, x = dadoSimulado,
      method = "BFGS"
    ),
    TRUE
  )
  if (inherits(fit_gtdl0, "try-error")) return(list(ok = FALSE))
  
  fit_zero0 <- try(
    optim(
      par = c(-1),                                   # gamma0 only
      fn = veroZero_H0, w = dado0,
      method = "BFGS"
    ),
    TRUE
  )
  if (inherits(fit_zero0, "try-error")) return(list(ok = FALSE))
  
  ll_red <- -( fit_gtdl0$value + fit_zero0$value )
  
  # LR test
  LR <- 2 * (ll_full - ll_red)
  df <- 3
  pval <- 1 - pchisq(LR, df = df)
  reject <- as.integer(pval < alpha_sig)
  
  list(ok = TRUE, LR = LR, pval = pval, reject = reject)
}






##standard GTDL model (no frailty term)

veroGTDL_nofrail <- function(x, par) {
  # par = c(alpha, log(lambda), beta)
  alpha  <- par[1]
  lambda <- exp(par[2])
  beta   <- par[3]
  
  cens  <- x[x$tempo > 0, 1]
  tempo <- x[x$tempo > 0, 2]
  X     <- as.matrix(x[x$tempo > 0, 3])
  
  eta <- as.numeric(alpha * tempo + X %*% beta)
  xb  <- as.numeric(X %*% beta)
  
  # log(1 + exp(z)) 
  log1pexp <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  
  # logA = log((1+exp(eta))/(1+exp(xb)))
  logA <- log1pexp(eta) - log1pexp(xb)
  
  # logS0 e logf0
  logS <- -(lambda / alpha) * logA
  logf <- log(lambda) + eta - log1pexp(eta) + logS
  
  ll <- sum(cens * logf + (1 - cens) * logS)
  return(-ll)  # negative log-likelihood
}


#alpha=alpha0+x*alpha1

veroGTDL_alpha_nofrail <- function(x, par) {
  # par = c(alpha0, alpha1, log(lambda), beta1)
  alpha0 <- par[1]
  alpha1 <- par[2]
  lambda <- exp(par[3])
  beta1  <- par[4]
  
  cens  <- x[x$tempo > 0, 1]
  tempo <- x[x$tempo > 0, 2]
  X     <- as.numeric(x[x$tempo > 0, 3])
  
  log1pexp <- function(z) ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  
  alpha_i <- alpha0 + alpha1 * X
  xb      <- beta1 * X
  eta     <- alpha_i * tempo + xb
  
  logA <- log1pexp(eta) - log1pexp(xb)
  logS <- -(lambda / alpha_i) * logA
  logf <- log(lambda) + eta - log1pexp(eta) + logS
  
  ll <- sum(cens * logf + (1 - cens) * logS)
  return(-ll)
}





one_rep_LR_theta0 <- function(nk, params_fixed, interval_root = c(0, 1080), alpha_sig = 0.05) {
  
  # fixed params (DATA GENERATION)
  gamma0 <- params_fixed$gamma0
  gamma1 <- params_fixed$gamma1
  alpha0 <- params_fixed$alpha0
  alpha1 <- params_fixed$alpha1
  lambda <- params_fixed$lambda
  theta  <- params_fixed$theta   # true theta used to generate
  beta1  <- params_fixed$beta1
  
  # --- helpers
  p0f <- function(w, gamma) {
    w <- as.matrix(w)
    W <- cbind(1, w)
    as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
  }
  
  # ---- survival GTDL NO-FRAIL (theta=0) for generator
  S0_GTDL_alpha <- function(t, x, alpha0, alpha1, lambda, beta1) {
    alpha_i <- alpha0 + alpha1 * x
    xb <- beta1 * x
    eta <- alpha_i * t + xb
    ((1 + exp(eta)) / (1 + exp(xb)))^(-lambda / alpha_i)
  }
  
  # cure prob under no-frail when alpha_i < 0: p(x) = lim_{t->inf} S0(t|x)
  p_cure_GTDL_alpha <- function(x, alpha0, alpha1, lambda, beta1) {
    alpha_i <- alpha0 + alpha1 * x
    xb <- beta1 * x
    if (alpha_i < 0) (1 + exp(xb))^(lambda / alpha_i) else 0
  }
  
  # ---- frailty generator pieces (your original)
  p1f_alpha <- function(x, par) {
    gamma <- par[1:2]
    alpha0 <- par[3]; alpha1 <- par[4]
    lambda <- par[5]; theta  <- par[6]; beta1 <- par[7]
    
    X <- as.matrix(x)
    W <- cbind(1, X)
    p0 <- as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
    
    alpha_i <- as.numeric(alpha0 + alpha1 * X)
    xb      <- as.numeric(X %*% beta1)
    
    p_star <- (1 - (lambda*theta/alpha_i) * log(1 + exp(xb)))^(-1/theta)
    (1 - p0) * p_star
  }
  
  sobrevGTDL_alpha <- function(t, x, par) {
    gamma <- par[1:2]
    alpha0 <- par[3]; alpha1 <- par[4]
    lambda <- par[5]; theta <- par[6]; beta1 <- par[7]
    
    X <- as.matrix(x)
    W <- cbind(1, X)
    p0 <- as.numeric(exp(W %*% gamma) / (1 + exp(W %*% gamma)))
    
    alpha_i <- as.numeric(alpha0 + alpha1 * X)
    xb  <- as.numeric(X %*% beta1)
    eta <- as.numeric(alpha_i * t + xb)
    
    A    <- log((1 + exp(eta)) / (1 + exp(xb)))
    Spos <- (1 + (lambda*theta/alpha_i) * A)^(-1/theta)
    
    (1 - p0) * Spos
  }
  
  # ---------- DATA GENERATION ----------
  # covariate (group indicator)
  X <- rbinom(nk, size = 1, prob = 0.5)
  
  t <- numeric(nk)
  delta <- integer(nk)
  
  if (theta <= 0) {
    # ===== Correct generator under H0: theta = 0 (NO-FRAIL) =====
    for (i in 1:nk) {
      x <- X[i]
      
      p0 <- p0f(x, gamma = c(gamma0, gamma1))
      pc <- p_cure_GTDL_alpha(x, alpha0, alpha1, lambda, beta1)  # p(x) in (0,1) if alpha_i<0
      
      u <- runif(1)
      
      tf <- if (u <= p0) {
        0
      } else if (u > 1 - (1 - p0) * pc) {
        Inf
      } else {
        # Inversion:
        # u = p0 + (1-p0)*(1 - S0(t))  =>  S0(t) = 1 - (u - p0)/(1-p0)
        s_target <- 1 - (u - p0) / (1 - p0)
        
        uniroot(
          f = function(tt) S0_GTDL_alpha(tt, x, alpha0, alpha1, lambda, beta1) - s_target,
          interval = interval_root
        )$root
      }
      
      tc <- runif(1, 0, 65)
      t[i] <- min(tf, tc)
      delta[i] <- as.integer(is.finite(tf) && tf <= tc)
    }
    
  } else {
    # ===== Your original generator for theta > 0 =====
    p0 <- p0f(X, gamma = c(gamma0, gamma1))
    p1 <- p1f_alpha(X, par = c(gamma0, gamma1, alpha0, alpha1, lambda, theta, beta1))
    
    froot_alpha <- function(t, x, par, u) sobrevGTDL_alpha(t, x, par) - 1 + u
    
    for (i in 1:nk) {
      u1 <- runif(1)
      u2 <- runif(1, min = p0[i], max = (1 - p1[i]))
      
      tf <- if (u1 <= p0[i]) {
        0
      } else if ((1 - p1[i]) <= u1) {
        Inf
      } else {
        uniroot(
          f = froot_alpha,
          x = X[i],
          u = u2,
          par = c(gamma0, gamma1, alpha0, alpha1, lambda, theta, beta1),
          interval = interval_root
        )$root
      }
      
      tc <- runif(1, 0, 65)
      
      t[i] <- min(tf, tc)
      delta[i] <- as.integer(t[i] == tf)
    }
  }
  
  dadoSimulado <- data.frame(cens = delta, tempo = t, x1 = X)
  dado0        <- data.frame(tempo = t, w0 = rep(1, nk), w1 = X)
  
  # ---------- UNRESTRICTED FIT (theta livre) ----------
  fit_gtdl <- try(
    optim(
      par = c(-0.1, 0.0, log(0.5), log(0.4), 0.0),   # (alpha0, alpha1, loglambda, logtheta, beta1)
      fn = veroGTDL_alpha, x = dadoSimulado,
      method = "BFGS",
      control = list(maxit = 500)
    ),
    TRUE
  )
  if (inherits(fit_gtdl, "try-error") || fit_gtdl$convergence != 0) return(list(ok = FALSE))
  
  fit_zero <- try(
    optim(
      par = c(-1, 0),
      fn = veroZero, w = dado0,
      method = "BFGS",
      control = list(maxit = 500)
    ),
    TRUE
  )
  if (inherits(fit_zero, "try-error") || fit_zero$convergence != 0) return(list(ok = FALSE))
  
  ll_full <- -( fit_gtdl$value + fit_zero$value )
  
  # ---------- RESTRICTED FIT under H0: theta = 0 (sem fragilidade) ----------
  fit_gtdl0 <- try(
    optim(
      par = c(-0.1, 0.0, log(0.5), 0.0),             # (alpha0, alpha1, loglambda, beta1)
      fn = veroGTDL_alpha_nofrail, x = dadoSimulado,
      method = "BFGS",
      control = list(maxit = 500)
    ),
    TRUE
  )
  if (inherits(fit_gtdl0, "try-error") || fit_gtdl0$convergence != 0) return(list(ok = FALSE))
  
  ll_red <- -( fit_gtdl0$value + fit_zero$value )
  
  # ---------- LR test (theta na fronteira) ----------
  LR <- 2 * (ll_full - ll_red)
  if (!is.finite(LR) || LR < 0) return(list(ok = FALSE))
  
  pval <- 0.5 * (1 - pchisq(LR, df = 1))
  reject <- as.integer(pval < alpha_sig)
  
  list(ok = TRUE, LR = LR, pval = pval, reject = reject)
}





one_rep_LR_theta0_retry <- function(nk, params_fixed, interval_root = c(0, 1080),
                                    alpha_sig = 0.05, max_tries = 200) {
  tries <- 0L
  repeat {
    tries <- tries + 1L
    
    out <- one_rep_LR_theta0(
      nk = nk,
      params_fixed = params_fixed,
      interval_root = interval_root,
      alpha_sig = alpha_sig
    )
    
    if (isTRUE(out$ok)) {
      out$tries <- tries
      return(out)
    }
    
    if (tries >= max_tries) {
      # retorna um "hard fail" (você pode decidir parar tudo se isso acontecer)
      return(list(ok = FALSE, tries = tries))
    }
  }
}

