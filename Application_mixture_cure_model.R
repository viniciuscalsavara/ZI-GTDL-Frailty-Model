############################
## Mixture cure rate model
############################


#Functions#

S_weib <- function(t, shape, scale) {
  exp(- (t / scale)^shape)
}

f_weib <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1) * exp(- (t / scale)^shape)
}


loglik_mix_ZI_weibull <- function(par, tempo, status, X0, X1, Xb, eps = 1e-12) {
  
  p0dim <- ncol(X0)
  p1dim <- ncol(X1)
  bdim  <- ncol(Xb)
  
  # desempacota
  i1 <- 1
  beta0 <- par[i1:(i1+p0dim-1)]; i1 <- i1+p0dim
  beta1 <- par[i1:(i1+p1dim-1)]; i1 <- i1+p1dim
  log_shape <- par[i1];          i1 <- i1+1
  beta_scale <- par[i1:(i1+bdim-1)]
  
  shape <- exp(log_shape)  # k>0
  scale <- exp(drop(Xb %*% beta_scale))  # lambda(x)>0
  
  # multinomial logit para (p0, p1, ps)
  z0 <- drop(X0 %*% beta0)
  z1 <- drop(X1 %*% beta1)
  den <- 1 + exp(z0) + exp(z1)
  p0  <- exp(z0) / den
  p1  <- exp(z1) / den
  ps  <- 1 / den  # suscetível
  
  t <- tempo
  d <- status
  
  ll <- 0
  
  # t=0 (evento em zero)
  idx0 <- which(t == 0)
  if (length(idx0) > 0) {
    ll <- ll + sum(log(pmax(p0[idx0], eps)))
  }
  
  # t>0
  idx <- which(t > 0)
  if (length(idx) > 0) {
    ti <- t[idx]; di <- d[idx]
    p0i <- p0[idx]; p1i <- p1[idx]; psi <- ps[idx]
    sc  <- scale[idx]
    
    S0 <- S_weib(ti, shape = shape, scale = sc)
    f0 <- f_weib(ti, shape = shape, scale = sc)
    
    # evento
    if (any(di == 1)) {
      ll <- ll + sum(log(pmax(psi[di==1], eps)) + log(pmax(f0[di==1], eps)))
    }
    
    # censura
    if (any(di == 0)) {
      surv_mix <- p1i[di==0] + psi[di==0] * S0[di==0]
      ll <- ll + sum(log(pmax(surv_mix, eps)))
    }
  }
  
  return(-ll)
}




`%||%` <- function(a, b) if (!is.null(a)) a else b

summary_optim_mix <- function(fit, X0, X1, Xb){
  
  est <- fit$par
  V   <- solve(fit$hessian)
  se  <- sqrt(diag(V))
  
  z  <- est / se
  p  <- 2 * (1 - pnorm(abs(z)))
  
  ci_lo <- est - 1.96 * se
  ci_hi <- est + 1.96 * se
  
  
  p0dim <- ncol(X0)
  p1dim <- ncol(X1)
  bdim  <- ncol(Xb)
  
  nm <- c(
    paste0("beta0_p0:", colnames(X0) %||% paste0("X0_", seq_len(p0dim))),
    paste0("beta1_p1:", colnames(X1) %||% paste0("X1_", seq_len(p1dim))),
    "log_shape",
    paste0("beta_scale:", colnames(Xb) %||% paste0("Xb_", seq_len(bdim)))
  )
  
  tab <- data.frame(
    Parameter = nm,
    Estimate  = est,
    SE        = se,
    CI_lower  = ci_lo,
    CI_upper  = ci_hi,
    z_value   = z,
    p_value   = p
  )
  
  # shape na escala natural
  idx_shape <- p0dim + p1dim + 1
  log_shape_hat <- est[idx_shape]
  se_log_shape  <- se[idx_shape]
  
  shape_hat <- exp(log_shape_hat)
  se_shape  <- shape_hat * se_log_shape
  z_shape   <- log_shape_hat / se_log_shape
  p_shape   <- 2 * (1 - pnorm(abs(z_shape)))
  
  tab_shape_nat <- data.frame(
    Parameter = "shape (exp(log_shape))",
    Estimate  = shape_hat,
    SE        = se_shape,
    CI_lower  = exp(log_shape_hat - 1.96 * se_log_shape),
    CI_upper  = exp(log_shape_hat + 1.96 * se_log_shape),
    z_value   = z_shape,
    p_value   = p_shape
  )
  
  tab_final <- rbind(tab, tab_shape_nat)
  
  # arredondar para 3 casas
  tab_final[, -1] <- round(tab_final[, -1], 3)
  
  tab_final
}




S_weib <- function(t, shape, scale) exp(-(t/scale)^shape)

S_mix_weib <- function(t, X0row, X1row, Xbrow, par){
  # par = c(beta0, beta1, log_shape, beta_scale)
  p0dim <- length(X0row)
  p1dim <- length(X1row)
  bdim  <- length(Xbrow)
  
  beta0 <- par[1:p0dim]
  beta1 <- par[(p0dim+1):(p0dim+p1dim)]
  log_shape <- par[p0dim+p1dim+1]
  beta_scale <- par[(p0dim+p1dim+2):(p0dim+p1dim+1+bdim)]
  
  shape <- exp(log_shape)
  scale <- exp(sum(Xbrow * beta_scale))
  
  z0 <- sum(X0row * beta0)
  z1 <- sum(X1row * beta1)
  den <- 1 + exp(z0) + exp(z1)
  
  p0 <- exp(z0)/den
  p1 <- exp(z1)/den
  ps <- 1/den
  
  # Para t>0: S(t)=p1 + ps*S0(t)
  p1 + ps * S_weib(t, shape = shape, scale = scale)
}




########
# Start
########


library(ggplot2)
library(survival)
library(msm) 
library(latex2exp)
library(tidyr)
library(survminer)
library(gridExtra)
library(dplyr)
library(gt)
library(patchwork)

dados <- read.table(file = "dados_censura_versao4.txt", header = TRUE) %>%
  mutate(
    x1=as.factor(x1),
    x2=as.factor(x2),
    x3=as.factor(x3),
    x4=as.factor(x4),
    status2=ifelse(status==1,0,1)
  )


tempo  <- dados$tempo
status <- dados$status

##X1 


Xall <- as.matrix(cbind(1, dados$x1))

X0 <- Xall
X1 <- Xall
Xb <- Xall  


p0dim <- ncol(X0); p1dim <- ncol(X1); bdim <- ncol(Xb)

par0 <- c(
  rep(-2, p0dim),  
  rep(-2, p1dim),  
  log(1.0),        
  rep(0.0, bdim)   
)

fit_mix_weib <- optim(
  par = par0,
  fn  = loglik_mix_ZI_weibull,
  tempo = tempo,
  status = status,
  X0 = X0, X1 = X1, Xb = Xb,
  method = "BFGS",
  hessian = TRUE
)

fit_mix_weib$convergence
fit_mix_weib$value   # -logLik
k <- length(fit_mix_weib$par)
AIC <- 2*fit_mix_weib$value + 2*k
AIC


tab_res <- summary_optim_mix(fit_mix_weib, X0 = X0, X1 = X1, Xb = Xb)
print(tab_res, row.names = FALSE)





# KM por X1

km <- survfit(Surv(tempo, status) ~ x1, data = dados)

km_df <- data.frame(
  tempo = km$time,
  surv  = km$surv,
  x1    = rep(as.numeric(gsub("x1=", "", names(km$strata))), km$strata),
  metodo = "Kaplan-Meier"
)

tgrid <- seq(0, max(dados$tempo), length.out = 800)

Xrow_0 <- c(1, 0)
Xrow_1 <- c(1, 1)

fit0 <- sapply(tgrid, function(tt) S_mix_weib(tt, X0row = Xrow_0, X1row = Xrow_0, Xbrow = Xrow_0, par = fit_mix_weib$par))
fit1 <- sapply(tgrid, function(tt) S_mix_weib(tt, X0row = Xrow_1, X1row = Xrow_1, Xbrow = Xrow_1, par = fit_mix_weib$par))

fit_df <- rbind(
  data.frame(tempo = tgrid, surv = fit0, x1 = 0, metodo = "Mixture Weibull"),
  data.frame(tempo = tgrid, surv = fit1, x1 = 1, metodo = "Mixture Weibull")
)

df_all <- bind_rows(km_df, fit_df)

# Plot
p1 <- ggplot(df_all, aes(x = tempo, y = surv)) +
  geom_step(data = subset(df_all, metodo == "Kaplan-Meier"),
            aes(color = factor(x1), linetype = metodo),
            linewidth = 1) +
  geom_line(data = subset(df_all, metodo == "Mixture Weibull"),
            aes(color = factor(x1), linetype = metodo),
            linewidth = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (month)", y = "S(t)",
       color = "X1",
       linetype = "Estimation method") +
  scale_linetype_manual(values = c("Kaplan-Meier" = "solid",
                                   "Mixture Weibull" = "dashed")) +
  theme_classic(base_size = 16)

print(p1)







##X2 

dados$x2 <- factor(dados$x2)   

Xall <- model.matrix(~ x2, data = dados)  

X0 <- Xall
X1 <- Xall
Xb <- Xall  

p0dim <- ncol(X0); p1dim <- ncol(X1); bdim <- ncol(Xb)

par0 <- c(
  rep(-2, p0dim),  
  rep(-2, p1dim),   
  log(1.0),       
  rep(0.0, bdim)    
)

fit_mix_weib <- optim(
  par = par0,
  fn  = loglik_mix_ZI_weibull,
  tempo = tempo,
  status = status,
  X0 = X0, X1 = X1, Xb = Xb,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 2000)
)

fit_mix_weib$convergence
fit_mix_weib$message
fit_mix_weib$value
AIC <- 2*fit_mix_weib$value + 2*length(fit_mix_weib$par)
AIC

tab_res <- summary_optim_mix(fit_mix_weib, X0 = X0, X1 = X1, Xb = Xb)
print(tab_res, row.names = FALSE)




# KM por X2
dados$x2 <- factor(dados$x2)

# KM por X2
km <- survfit(Surv(tempo, status) ~ x2, data = dados)

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  x2     = rep(gsub("^x2=", "", names(km$strata)), km$strata),
  metodo = "Kaplan-Meier"
)


tgrid <- seq(0, max(dados$tempo), length.out = 800)

levs <- levels(dados$x2)


form_x2 <- as.formula("~ x2")

make_Xrow <- function(level, X_template_cols) {
  nd <- data.frame(x2 = factor(level, levels = levs))
  Xr <- model.matrix(form_x2, data = nd)
  Xr <- Xr[, X_template_cols, drop = FALSE]
  as.numeric(Xr)
}

X0_cols <- colnames(X0)
X1_cols <- colnames(X1)
Xb_cols <- colnames(Xb)

fit_list <- lapply(levs, function(lv){
  
  X0row <- make_Xrow(lv, X0_cols)
  X1row <- make_Xrow(lv, X1_cols)
  Xbrow <- make_Xrow(lv, Xb_cols)
  
  surv_fit <- sapply(tgrid, function(tt)
    S_mix_weib(tt, X0row = X0row, X1row = X1row, Xbrow = Xbrow, par = fit_mix_weib$par)
  )
  
  data.frame(tempo = tgrid, surv = surv_fit, x2 = lv, metodo = "Mixture Weibull")
})

fit_df <- bind_rows(fit_list)

df_all <- bind_rows(km_df, fit_df)

# Plot
p2 <- ggplot(df_all, aes(x = tempo, y = surv)) +
  geom_step(data = subset(df_all, metodo == "Kaplan-Meier"),
            aes(color = x2, linetype = metodo),
            linewidth = 1) +
  geom_line(data = subset(df_all, metodo == "Mixture Weibull"),
            aes(color = x2, linetype = metodo),
            linewidth = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Time (month)", y = "S(t)",
       color = "X2",
       linetype = "Estimation method") +
  scale_linetype_manual(values = c("Kaplan-Meier" = "solid",
                                   "Mixture Weibull" = "dashed")) +
  theme_classic(base_size = 16)

print(p2)








# -------------------------
# X3
# -------------------------


dados$x3 <- as.numeric(as.character(dados$x3))

Xall <- as.matrix(cbind(1, dados$x3))
colnames(Xall) <- c("(Intercept)", "x3")

X0 <- Xall
X1 <- Xall
Xb <- Xall  

p0dim <- ncol(X0); p1dim <- ncol(X1); bdim <- ncol(Xb)

par0 <- c(
  rep(-2, p0dim),  
  rep(-2, p1dim),  
  log(1.0),        
  rep(0.0, bdim)   
)

fit_mix_weib <- optim(
  par = par0,
  fn  = loglik_mix_ZI_weibull,
  tempo = tempo,
  status = status,
  X0 = X0, X1 = X1, Xb = Xb,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 2000)
)

fit_mix_weib$convergence
fit_mix_weib$message
fit_mix_weib$value   # -logLik
AIC <- 2*fit_mix_weib$value + 2*length(fit_mix_weib$par)
AIC

tab_res <- summary_optim_mix(fit_mix_weib, X0 = X0, X1 = X1, Xb = Xb)
print(tab_res, row.names = FALSE)

# -------------------------
# KM vs Mixture Weibull: X3
# -------------------------

km <- survfit(Surv(tempo, status) ~ x3, data = dados)

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  x3     = rep(as.numeric(gsub("^x3=", "", names(km$strata))), km$strata),
  metodo = "Kaplan-Meier"
)

tgrid <- seq(0, max(dados$tempo), length.out = 800)

Xrow_0 <- c(1, 0)
Xrow_1 <- c(1, 1)

fit0 <- sapply(tgrid, function(tt)
  S_mix_weib(tt, X0row = Xrow_0, X1row = Xrow_0, Xbrow = Xrow_0, par = fit_mix_weib$par)
)
fit1 <- sapply(tgrid, function(tt)
  S_mix_weib(tt, X0row = Xrow_1, X1row = Xrow_1, Xbrow = Xrow_1, par = fit_mix_weib$par)
)

fit_df <- rbind(
  data.frame(tempo = tgrid, surv = fit0, x3 = 0, metodo = "Mixture Weibull"),
  data.frame(tempo = tgrid, surv = fit1, x3 = 1, metodo = "Mixture Weibull")
)

df_all <- bind_rows(km_df, fit_df)

p3 <- ggplot(df_all, aes(x = tempo, y = surv)) +
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = factor(x3), linetype = metodo),
    linewidth = 1
  ) +
  geom_line(
    data = subset(df_all, metodo == "Mixture Weibull"),
    aes(color = factor(x3), linetype = metodo),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Time (month)", y = "S(t)",
    color = "X3",
    linetype = "Estimation method"
  ) +
  scale_linetype_manual(values = c("Kaplan-Meier" = "solid",
                                   "Mixture Weibull" = "dashed")) +
  theme_classic(base_size = 16)

print(p3)






##X4


dados$x4 <- as.numeric(as.character(dados$x4))

# -------------------------
# X4
# -------------------------

Xall <- as.matrix(cbind(1, dados$x4))
colnames(Xall) <- c("(Intercept)", "x4")

X0 <- Xall
X1 <- Xall
Xb <- Xall 

p0dim <- ncol(X0); p1dim <- ncol(X1); bdim <- ncol(Xb)

par0 <- c(
  rep(-2, p0dim),  
  rep(-2, p1dim),  
  log(1.0),        
  rep(0.0, bdim)  
)

fit_mix_weib <- optim(
  par = par0,
  fn  = loglik_mix_ZI_weibull,
  tempo = tempo,
  status = status,
  X0 = X0, X1 = X1, Xb = Xb,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 2000)
)

fit_mix_weib$convergence
fit_mix_weib$message
fit_mix_weib$value   # -logLik
AIC <- 2*fit_mix_weib$value + 2*length(fit_mix_weib$par)
AIC

tab_res <- summary_optim_mix(fit_mix_weib, X0 = X0, X1 = X1, Xb = Xb)
print(tab_res, row.names = FALSE)

# -------------------------
# KM vs Mixture Weibull: X4
# -------------------------

km <- survfit(Surv(tempo, status) ~ x4, data = dados)

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  x4     = rep(as.numeric(gsub("^x4=", "", names(km$strata))), km$strata),
  metodo = "Kaplan-Meier"
)

tgrid <- seq(0, max(dados$tempo), length.out = 800)

Xrow_0 <- c(1, 0)
Xrow_1 <- c(1, 1)

fit0 <- sapply(tgrid, function(tt)
  S_mix_weib(tt, X0row = Xrow_0, X1row = Xrow_0, Xbrow = Xrow_0, par = fit_mix_weib$par)
)
fit1 <- sapply(tgrid, function(tt)
  S_mix_weib(tt, X0row = Xrow_1, X1row = Xrow_1, Xbrow = Xrow_1, par = fit_mix_weib$par)
)

fit_df <- rbind(
  data.frame(tempo = tgrid, surv = fit0, x4 = 0, metodo = "Mixture Weibull"),
  data.frame(tempo = tgrid, surv = fit1, x4 = 1, metodo = "Mixture Weibull")
)

df_all <- bind_rows(km_df, fit_df)

p4 <- ggplot(df_all, aes(x = tempo, y = surv)) +
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = factor(x4), linetype = metodo),
    linewidth = 1
  ) +
  geom_line(
    data = subset(df_all, metodo == "Mixture Weibull"),
    aes(color = factor(x4), linetype = metodo),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Time (month)", y = "S(t)",
    color = "X4",
    linetype = "Estimation method"
  ) +
  scale_linetype_manual(values = c("Kaplan-Meier" = "solid",
                                   "Mixture Weibull" = "dashed")) +
  theme_classic(base_size = 16)

print(p4)







#Multivariable mixture model#

Xall <- as.matrix(cbind(1, dados$x1, dados$x2_1, dados$x2_2, dados$x3, dados$x4))

X0 <- Xall
X1 <- Xall
Xb <- Xall  


p0dim <- ncol(X0); p1dim <- ncol(X1); bdim <- ncol(Xb)

par0 <- c(
  rep(-2, p0dim),  
  rep(-2, p1dim),  
  log(1.0),       
  rep(0.0, bdim)   
)

fit_mix_weib <- optim(
  par = par0,
  fn  = loglik_mix_ZI_weibull,
  tempo = tempo,
  status = status,
  X0 = X0, X1 = X1, Xb = Xb,
  method = "BFGS",
  hessian = TRUE
)

fit_mix_weib$convergence
fit_mix_weib$value   # -logLik
k <- length(fit_mix_weib$par)
AIC <- 2*fit_mix_weib$value + 2*k
AIC
