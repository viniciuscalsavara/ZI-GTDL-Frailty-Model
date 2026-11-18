################################
# Hypothesis testing: H0: theta=0 (LR with boundary correction)
################################

library(foreach)
library(doParallel)
library(doRNG)

registerDoRNG(2024)
source("functions.R")

## ---- 1) Fixed parameters (base) ----
params_base <- list(
  gamma0=-1, gamma1=-0.6, alpha0=-0.1, alpha1=0,
  lambda=0.5, beta1=-0.9
)

thetas <- c(0)
n <- c(200, 1000)

repeticao <- 1000
alpha_sig <- 0.05
max_tries <- 200   # ajuste se precisar

## ---- 2) Storage: arrays [rep, n, theta] ----
LR_stat <- array(NA_real_, dim = c(repeticao, length(n), length(thetas)))
pval_LR <- array(NA_real_, dim = c(repeticao, length(n), length(thetas)))
rej_LR  <- array(NA_integer_, dim = c(repeticao, length(n), length(thetas)))
n_try   <- array(NA_integer_, dim = c(repeticao, length(n), length(thetas)))

## ---- 3) Parallel backend ----
ncores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

## ---- 4) Loops: theta scenario -> sample size -> parallel reps ----
for (th in seq_along(thetas)) {
  cat("\n===== theta TRUE =", thetas[th], "=====\n")
  
  params_fixed <- params_base
  params_fixed$theta <- thetas[th]
  
  for (k in seq_along(n)) {
    cat("n =", n[k], "\n")
    
    res <- foreach(j = 1:repeticao,
                   .packages = c("stats"),
                   .export = c(
                     ## wrapper (new)
                     "one_rep_LR_theta0_retry",
                     ## main function
                     "one_rep_LR_theta0",
                     ## likelihood pieces it calls
                     "veroGTDL_alpha", "veroGTDL_alpha_nofrail",
                     "veroZero"
                   )) %dorng% {
                     one_rep_LR_theta0_retry(
                       nk = n[k],
                       params_fixed = params_fixed,
                       alpha_sig = alpha_sig,
                       max_tries = max_tries
                     )
                   }
    
    # checa se houve "hard fails" (atingiu max_tries)
    ok_vec <- vapply(res, `[[`, logical(1), "ok")
    if (any(!ok_vec)) {
      stop("Algumas replicações atingiram max_tries sem convergir. Aumente max_tries ou revise chutes/otimização.")
    }
    
    for (j in 1:repeticao) {
      LR_stat[j, k, th] <- res[[j]]$LR
      pval_LR[j, k, th] <- res[[j]]$pval
      rej_LR[j, k, th]  <- res[[j]]$reject
      n_try[j, k, th]   <- res[[j]]$tries
    }
    
    cat("  mean tries =", round(mean(n_try[, k, th]), 2),
        " | max tries =", max(n_try[, k, th]), "\n")
  }
}

parallel::stopCluster(cl)

## ---- 5) Summary tables ----
summ_tbl <- do.call(rbind, lapply(seq_along(thetas), function(th) {
  data.frame(
    theta = thetas[th],
    n     = n,
    B     = repeticao,
    power = apply(rej_LR[, , th], 2, mean),
    mean_tries = apply(n_try[, , th], 2, mean),
    max_tries  = apply(n_try[, , th], 2, max),
    stringsAsFactors = FALSE
  )
}))

type1_tbl <- subset(summ_tbl, theta == 1e-6)
power_tbl <- subset(summ_tbl, theta != 1e-6)

type1_tbl$type1_pct <- round(type1_tbl$power, 3) * 100
power_tbl$power_pct <- round(power_tbl$power, 3) * 100

type1_tbl
power_tbl
