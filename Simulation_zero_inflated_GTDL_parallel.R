
rm(list = ls())
library(msm)          # deltamethod
library(doParallel)
library(foreach)
library(xtable)


n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# (optional) reproducible parallel RNG
library(doRNG)
registerDoRNG(2024)

source("functions.R")



set.seed(2024)


# fixed parameters
params_fixed <- list(
  gama0 = -1.0, gama1 = 0,
  alpha = 0, lambda = 0.5, theta = 1.5, beta1 = 0
)


n <- c(100, 200, 300, 500, 800, 1000, 2000)
repeticao <- 1000

EMV_gtdl <- array(NA, dim = c(repeticao, 4, 2, length(n)))
EMV_zero <- array(NA, dim = c(repeticao, 2, 2, length(n)))
dp_p0    <- array(NA, dim = c(repeticao, 2, length(n)))
dp_p1    <- array(NA, dim = c(repeticao, 2, length(n)))

erro_optim1 <- erro_optim2 <- erro_optim3 <- 0
erro_optim11 <- erro_optim12 <- erro_optim13 <- 0

for (k in seq_along(n)) {
  cat("n = ", n[k], "\n")
  
  res <- foreach(j = 1:repeticao,
                 .packages = c("msm"),
                 .export = c("veroZero","veroGTDL","sobrevGTDLv1","p0f","p1f","froot","one_rep")
  ) %dorng% {
    one_rep(nk = n[k], params_fixed = params_fixed)
  }
  
  # keep only successful fits (like your while-loop logic)
  ok_idx <- which(vapply(res, `[[`, logical(1), "ok"))
  
  if (length(ok_idx) < repeticao) {
    # If you want EXACTLY repeticao successful runs, you can rerun failures;
    # but simplest is: fill what you have, keep NAs for failures.
    cat("  successful =", length(ok_idx), "of", repeticao, "\n")
  }
  
  # errors
  err_mat <- do.call(rbind, lapply(res, `[[`, "err"))
  erro_optim1  <- erro_optim1  + sum(err_mat[, "gtdl_fail"])
  erro_optim2  <- erro_optim2  + sum(err_mat[, "gtdl_solve"])
  erro_optim3  <- erro_optim3  + sum(err_mat[, "gtdl_negvar"])
  erro_optim11 <- erro_optim11 + sum(err_mat[, "zero_fail"])
  erro_optim12 <- erro_optim12 + sum(err_mat[, "zero_solve"])
  erro_optim13 <- erro_optim13 + sum(err_mat[, "zero_negvar"])
  
  # store results in j=1..repeticao positions (successful ones)
  for (j in ok_idx) {
    EMV_gtdl[j, , 1, k] <- res[[j]]$est_gtdl
    EMV_gtdl[j, , 2, k] <- res[[j]]$se_gtdl
    
    EMV_zero[j, , 1, k] <- res[[j]]$est_zero
    EMV_zero[j, , 2, k] <- res[[j]]$se_zero
    
    dp_p0[j, , k] <- res[[j]]$dp_p0
    dp_p1[j, , k] <- res[[j]]$dp_p1
  }
}

stopCluster(cl)









# parâmetros verdadeiros (use os mesmos que você fixou na simulação)
parametros <- c(params_fixed$alpha, params_fixed$lambda, params_fixed$theta, params_fixed$beta1, params_fixed$gama0, params_fixed$gama1)


alpha=params_fixed$alpha
lambda=params_fixed$lambda
theta=params_fixed$theta
beta1=params_fixed$beta1
gama0=params_fixed$gama0
gama1=params_fixed$gama1

# n_valid por tamanho amostral (quantas repetições realmente têm estimativas)
n_valid <- sapply(seq_along(n), function(i) sum(!is.na(EMV_gtdl[,1,1,i]) & !is.na(EMV_zero[,1,1,i])))
n_valid



vies <- matrix(NA, nrow = length(n), ncol = 10)

# 1) parâmetros do GTDL (alpha, lambda, theta, beta1)
for (j in 1:4) {
  for (i in 1:length(n)) {
    vies[i, j] <- mean(EMV_gtdl[, j, 1, i], na.rm = TRUE) - parametros[j]
  }
}

# 2) parâmetros da inflação (gama0, gama1)
for (j in 5:6) {
  for (i in 1:length(n)) {
    vies[i, j] <- mean(EMV_zero[, (j-4), 1, i], na.rm = TRUE) - parametros[j]
  }
}

# 3) vies de p0 (proporção de zeros) para x=1 e x=0
p0_verd <- c(p0f(1, c(gama0, gama1)), p0f(0, c(gama0, gama1)))

for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_zero[,1,1,i]) & !is.na(EMV_zero[,2,1,i]))
  if (length(idx) > 0) {
    soma1 <- sum(sapply(idx, function(k) p0f(1, EMV_zero[k,1:2,1,i])))
    soma0 <- sum(sapply(idx, function(k) p0f(0, EMV_zero[k,1:2,1,i])))
    vies[i, 7] <- soma1/length(idx) - p0_verd[1]
    vies[i, 8] <- soma0/length(idx) - p0_verd[2]
  }
}

# 4) vies de p1 (fração de cura) para x=1 e x=0
p1_verd <- c(
  p1f(1, c(gama0, gama1, alpha, lambda, theta, beta1)),
  p1f(0, c(gama0, gama1, alpha, lambda, theta, beta1))
)

for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_gtdl[,1,1,i]) & !is.na(EMV_zero[,1,1,i]))
  if (length(idx) > 0) {
    soma1 <- sum(sapply(idx, function(k)
      ifelse(EMV_gtdl[k,1,1,i] < 0,
             p1f(1, c(EMV_zero[k,1:2,1,i], EMV_gtdl[k,1:4,1,i])),
             NA_real_)
    ), na.rm = TRUE)
    
    soma0 <- sum(sapply(idx, function(k)
      ifelse(EMV_gtdl[k,1,1,i] < 0,
             p1f(0, c(EMV_zero[k,1:2,1,i], EMV_gtdl[k,1:4,1,i])),
             NA_real_)
    ), na.rm = TRUE)
    
    # denom: quantos realmente entraram no if(alpha<0)
    denom1 <- sum(sapply(idx, function(k) !is.na(ifelse(EMV_gtdl[k,1,1,i] < 0, 1, NA))))
    if (denom1 > 0) {
      vies[i, 9]  <- soma1/denom1 - p1_verd[1]
      vies[i, 10] <- soma0/denom1 - p1_verd[2]
    }
  }
}


xtable(t(vies), digits = 3)





eqm <- matrix(NA, nrow = length(n), ncol = 10)

# 1) GTDL
for (j in 1:4) {
  for (i in 1:length(n)) {
    eqm[i, j] <- sqrt(mean((EMV_gtdl[, j, 1, i] - parametros[j])^2, na.rm = TRUE))
  }
}

# 2) inflação
for (j in 5:6) {
  for (i in 1:length(n)) {
    eqm[i, j] <- sqrt(mean((EMV_zero[, (j-4), 1, i] - parametros[j])^2, na.rm = TRUE))
  }
}

# 3) EQM de p0
for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_zero[,1,1,i]) & !is.na(EMV_zero[,2,1,i]))
  if (length(idx) > 0) {
    soma1 <- mean(sapply(idx, function(k) (p0f(1, EMV_zero[k,1:2,1,i]) - p0_verd[1])^2))
    soma0 <- mean(sapply(idx, function(k) (p0f(0, EMV_zero[k,1:2,1,i]) - p0_verd[2])^2))
    eqm[i, 7] <- sqrt(soma1)
    eqm[i, 8] <- sqrt(soma0)
  }
}

# 4) EQM de p1
for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_gtdl[,1,1,i]) & !is.na(EMV_zero[,1,1,i]))
  if (length(idx) > 0) {
    vals1 <- sapply(idx, function(k)
      ifelse(EMV_gtdl[k,1,1,i] < 0,
             (p1f(1, c(EMV_zero[k,1:2,1,i], EMV_gtdl[k,1:4,1,i])) - p1_verd[1])^2,
             NA_real_)
    )
    vals0 <- sapply(idx, function(k)
      ifelse(EMV_gtdl[k,1,1,i] < 0,
             (p1f(0, c(EMV_zero[k,1:2,1,i], EMV_gtdl[k,1:4,1,i])) - p1_verd[2])^2,
             NA_real_)
    )
    eqm[i, 9]  <- sqrt(mean(vals1, na.rm = TRUE))
    eqm[i, 10] <- sqrt(mean(vals0, na.rm = TRUE))
  }
}

xtable(t(eqm), digits = 3)




PC <- matrix(NA, nrow = length(n), ncol = 10)

# GTDL
for (j in 1:4) {
  for (i in 1:length(n)) {
    est <- EMV_gtdl[, j, 1, i]
    se  <- EMV_gtdl[, j, 2, i]
    ok  <- !is.na(est) & !is.na(se)
    PC[i, j] <- mean((est[ok] - 1.96*se[ok] < parametros[j]) &
                       (est[ok] + 1.96*se[ok] > parametros[j]))
  }
}

# zero
for (j in 5:6) {
  for (i in 1:length(n)) {
    est <- EMV_zero[, (j-4), 1, i]
    se  <- EMV_zero[, (j-4), 2, i]
    ok  <- !is.na(est) & !is.na(se)
    PC[i, j] <- mean((est[ok] - 1.96*se[ok] < parametros[j]) &
                       (est[ok] + 1.96*se[ok] > parametros[j]))
  }
}



for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_zero[,1,1,i]) & !is.na(EMV_zero[,2,1,i]) &
                 !is.na(dp_p0[,1,i]) & !is.na(dp_p0[,2,i]))
  if (length(idx) > 0) {
    cover1 <- sapply(idx, function(k) {
      g0 <- EMV_zero[k,1,1,i]; g1 <- EMV_zero[k,2,1,i]
      LI <- exp(g0+g1)/(1+exp(g0+g1)) - 1.96*dp_p0[k,1,i]
      LS <- exp(g0+g1)/(1+exp(g0+g1)) + 1.96*dp_p0[k,1,i]
      (LI < p0_verd[1]) & (LS > p0_verd[1])
    })
    cover0 <- sapply(idx, function(k) {
      g0 <- EMV_zero[k,1,1,i]
      LI <- exp(g0)/(1+exp(g0)) - 1.96*dp_p0[k,2,i]
      LS <- exp(g0)/(1+exp(g0)) + 1.96*dp_p0[k,2,i]
      (LI < p0_verd[2]) & (LS > p0_verd[2])
    })
    PC[i, 7] <- mean(cover1)
    PC[i, 8] <- mean(cover0)
  }
}


#PC de p1 (x=1 e x=0) usando dp_p1

for (i in 1:length(n)) {
  idx <- which(!is.na(EMV_gtdl[,1,1,i]) & !is.na(EMV_zero[,1,1,i]) &
                 !is.na(dp_p1[,1,i]) & !is.na(dp_p1[,2,i]))
  if (length(idx) > 0) {
    
    cover1 <- sapply(idx, function(k) {
      a <- EMV_gtdl[k,1,1,i]
      if (is.na(a) || a >= 0) return(NA)
      lam <- EMV_gtdl[k,2,1,i]; th <- EMV_gtdl[k,3,1,i]; b1 <- EMV_gtdl[k,4,1,i]
      g0 <- EMV_zero[k,1,1,i]; g1 <- EMV_zero[k,2,1,i]
      p1hat <- (1 - exp(g0+g1)/(1+exp(g0+g1))) * (1 - ((lam*th)/a)*log(1+exp(b1)))^(-1/th)
      LI <- p1hat - 1.96*dp_p1[k,1,i]
      LS <- p1hat + 1.96*dp_p1[k,1,i]
      (LI < p1_verd[1]) & (LS > p1_verd[1])
    })
    
    cover0 <- sapply(idx, function(k) {
      a <- EMV_gtdl[k,1,1,i]
      if (is.na(a) || a >= 0) return(NA)
      lam <- EMV_gtdl[k,2,1,i]; th <- EMV_gtdl[k,3,1,i]
      g0 <- EMV_zero[k,1,1,i]
      p1hat <- (1 - exp(g0)/(1+exp(g0))) * (1 - ((lam*th)/a)*log(1+exp(0)))^(-1/th)
      LI <- p1hat - 1.96*dp_p1[k,2,i]
      LS <- p1hat + 1.96*dp_p1[k,2,i]
      (LI < p1_verd[2]) & (LS > p1_verd[2])
    })
    
    PC[i, 9]  <- mean(cover1, na.rm = TRUE)
    PC[i, 10] <- mean(cover0, na.rm = TRUE)
  }
}

xtable(t(PC), digits = 3)
