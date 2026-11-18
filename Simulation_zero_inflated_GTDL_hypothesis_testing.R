## =================================
## LR test Type I error 
## H0: alpha1 = beta1 = gamma1 = 0  
## =================================

library(foreach)
library(doParallel)
library(doRNG)
registerDoRNG(2024)
source("functions.R")

## ---- 1) Fixed parameters under H0 ----
params_fixed <- list(
  gamma0 = -1.0, gamma1 = -0.05,     # gamma1=0 under H0
  alpha0 = -0.1, alpha1 = -0.05,     # alpha1=0 under H0
  lambda =  0.5,
  theta  =  0.4,
  beta1  =  -0.05                    # beta1=0 under H0
)

n <- c(100, 200, 300, 500, 800, 1000, 2000,5000)
repeticao <- 1000
alpha_sig <- 0.05

## ---- 2) Storage ----
LR_stat <- array(NA_real_, dim = c(repeticao, length(n)))
pval_LR <- array(NA_real_, dim = c(repeticao, length(n)))
rej_LR  <- array(NA_integer_, dim = c(repeticao, length(n)))
ok_mat  <- array(FALSE, dim = c(repeticao, length(n)))

## Optional: count failures
erro_fit <- 0

## ---- 3) Parallel backend ----
ncores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)


## ---- 4) Loop over sample sizes (parallel over replicates) ----
for (k in seq_along(n)) {
  cat("n = ", n[k], "\n")
  
  res <- foreach(j = 1:repeticao,
                 .packages = c("stats"),
                 .export = c(
                   ## main function
                   "one_rep_LR_type1",
                   ## likelihood pieces it calls
                   "veroGTDL_alpha", "veroGTDL_alpha_H0",
                   "veroZero", "veroZero_H0"
                 )) %dorng% {
                   one_rep_LR_type1(
                     nk = n[k],
                     params_fixed = params_fixed,
                     alpha_sig = alpha_sig
                   )
                 }
  
  ## which runs succeeded
  ok_idx <- which(vapply(res, `[[`, logical(1), "ok"))
  ok_mat[ok_idx, k] <- TRUE
  
  if (length(ok_idx) < repeticao) {
    cat("  successful =", length(ok_idx), "of", repeticao, "\n")
    erro_fit <- erro_fit + (repeticao - length(ok_idx))
  }
  
  ## store results for successful runs
  for (j in ok_idx) {
    LR_stat[j, k] <- res[[j]]$LR
    pval_LR[j, k] <- res[[j]]$pval
    rej_LR[j, k]  <- res[[j]]$reject
  }
}

parallel::stopCluster(cl)

## ---- 5) Type I error summary table ----
type1_tbl <- data.frame(
  n     = n,
  B     = repeticao,
  n_ok  = colSums(ok_mat),
  type1 = colMeans(rej_LR, na.rm = TRUE),
  stringsAsFactors = FALSE
)

type1_tbl
erro_fit


round(type1_tbl$type1,3)*100


###Plot (type I error rate)


dat <- data.frame(
  effect = c(0, 0.01, 0.02, 0.03, 0.04, 0.05),
  n100   = c(7.2, 8.0,  9.9, 10.3, 12.9, 16.3),
  n200   = c(6.0, 7.0,  8.4, 10.6, 15.1, 19.9),
  n300   = c(4.8, 6.5, 10.7, 14.5, 20.2, 28.6),
  n500   = c(6.1, 8.5, 12.1, 19.1, 29.9, 43.9),
  n800   = c(5.3, 7.6, 16.3, 30.1, 46.8, 66.4),
  n1000  = c(5.3, 7.1, 15.0, 32.0, 54.6, 75.4),
  n2000  = c(5.3,10.7, 29.4, 60.0, 86.8, 97.4),
  n5000  = c(5.8,19.2, 66.5, 96.4,100.0,100.0)
)

## ----- Reshape to long format -----
long <- reshape(
  dat,
  varying = list(names(dat)[-1]),
  v.names = "rej",
  timevar = "n",
  times = as.numeric(gsub("n","", names(dat)[-1])),
  direction = "long"
)
row.names(long) <- NULL

## ----- Plot (base R; no extra packages) -----
long$effect_f <- factor(long$effect)

plot(NA, xlim = range(long$n), ylim = c(0, 100),
     xlab = "Sample size (n)", ylab = "Rejection rate (%)",
     log = "x")  # log-x helps spacing across 100..5000

for (ef in levels(long$effect_f)) {
  d <- long[long$effect_f == ef, ]
  lines(d$n, d$rej, type = "b")
}

legend("bottomright", legend = paste0("Effect = ", levels(long$effect_f)),
       bty = "n")

abline(h = 5, lty = 2)  # nominal 5% line (type I target)

