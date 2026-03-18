##################################
## Credit risk data - Application
##################################



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


df <- read.table(file = "dados_censura_versao4.txt", header = TRUE) %>%
  mutate(
    x1=as.factor(x1),
    x2=as.factor(x2),
    x3=as.factor(x3),
    x4=as.factor(x4),
    status2=ifelse(status==1,0,1)
  )




labs_x1 <- c(
  "Without consultation",
  "With consultation"
)

labs_x2 <- c(
  "Up to 3 years",
  "From 3 to 5 years",
  "Above 5 years"
)

labs_x3 <- c(
  "Does not have",
  "It has"
)

labs_x4 <- c(
  "Banks",
  "Other segments"
)


#median follow-up time

km1 <- survfit(Surv(time = tempo, event = status2) ~ 1, data = df)
km1
# K-M curves



ekm <- survfit(Surv(time = tempo, event = status) ~ 1, data = df)
ekm_df <- data.frame(tempo = ekm$time,
                     km = ekm$surv)

p<-ggplot(data = ekm_df, mapping = aes(tempo, km)) + 
  geom_step() + 
  labs(x = "Time (months)", y = "Estimated survival function") +
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(0,60,5)) +
  geom_hline(yintercept = 0.7623639, linetype = 2, linewidth = 0.5) +
  geom_hline(yintercept = 0.2161, linetype = 2, linewidth = 0.5) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    strip.text   = element_text(size = 18),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 18)
  )+
  annotate("text", x =1,  y = 0.8023,  label = expression(1-p[0]), size = 5) +
  annotate("text", x = 61, y = 0.2561, label = expression(p[1]), size = 5)

p

ggsave("KM_curve_version2.pdf", p, width = 9, height = 7, dpi = 300, bg = "white")






summarize_counts_by_cov <- function(df, cov) {
  df %>%
    mutate(
      cov_group = .data[[cov]],
      is_event  = (status == 1),
      is_t0     = (tempo == 0),
      is_post0  = (tempo > 0)
    ) %>%
    group_by(cov_group) %>%
    summarise(
      n_total = n(),
      n_event_t0 = sum(is_event & is_t0, na.rm = TRUE),
      n_event_post0 = sum(is_event & is_post0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(variable = cov) %>%
    select(variable, group = cov_group, n_total, n_event_t0, n_event_post0)
}

# Aplicar para as 4 covariáveis
tab_x1 <- summarize_counts_by_cov(df, "x1")
tab_x2 <- summarize_counts_by_cov(df, "x2")
tab_x3 <- summarize_counts_by_cov(df, "x3")
tab_x4 <- summarize_counts_by_cov(df, "x4")

tab_x1
tab_x2
tab_x3
tab_x4

# (Opcional) juntar tudo num único data.frame
tab_all <- bind_rows(tab_x1, tab_x2, tab_x3, tab_x4)
tab_all







# --- helper: KM + log-rank no estilo do seu código ---

km_by_cov <- function(df, cov, legend_title, legend_labels, x_max = 60) {
  
  fml <- as.formula(paste0("Surv(time = tempo, event = status) ~ ", cov))
  fit <- survfit(fml, data = df)
  
  strata_vec <- rep(names(fit$strata), fit$strata)
  km_df <- data.frame(
    tempo = fit$time,
    km    = fit$surv,
    strata = factor(strata_vec)
  )
  
  km_df$strata <- gsub(paste0("^", cov, "="), "", km_df$strata)
  
  lr <- survdiff(fml, data = df)
  pval <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  ptxt <- ifelse(pval < 0.001, "Log-rank p < 0.001",
                 paste0("Log-rank p = ", format.pval(pval, digits = 3, eps = 1e-3)))
  
  p <- ggplot(km_df, aes(x = tempo, y = km, group = strata, linetype = strata)) +
    geom_step(linewidth = 1) +
    labs(
      x = "Time (months)",
      y = "Estimated survival function",
      linetype = legend_title
    ) +
    scale_linetype_discrete(labels = legend_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    theme_bw(base_size = 18) +
    theme(
      # legenda dentro do painel:
      legend.position = c(0.98, 0.98),      # (x,y) de 0 a 1
      legend.justification = c(1, 1),       # ancora no canto superior direito
      legend.background = element_rect(),
      legend.key = element_rect(),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18)
    ) +
    annotate("text", x = 1, y = 0.05, label = ptxt, hjust = 0, size = 5)
  
  return(p)
}



# --- Rodar para as 4 covariáveis ---
p1 <- km_by_cov(df, "x1", "Consultation Information", labs_x1)
p2 <- km_by_cov(df, "x2", "Employment Time", labs_x2)
p3 <- km_by_cov(df, "x3", "Previous Settled Debts", labs_x3)
p4 <- km_by_cov(df, "x4", "Type of Debt", labs_x4)
# mostrar um por vez
p1
p2
p3
p4



# --- "Shared-looking" axes: y only on left column, x only on bottom row ---

# remove eixo-y (título + ticks + labels) dos painéis da direita
p2 <- p2 + theme(axis.title.y = element_blank(),
                 axis.text.y  = element_blank(),
                 axis.ticks.y = element_blank())

p4 <- p4 + theme(axis.title.y = element_blank(),
                 axis.text.y  = element_blank(),
                 axis.ticks.y = element_blank())

# remove eixo-x (título + ticks + labels) dos painéis de cima
p1 <- p1 + theme(axis.title.x = element_blank(),
                 axis.text.x  = element_blank(),
                 axis.ticks.x = element_blank())

p2 <- p2 + theme(axis.title.x = element_blank(),
                 axis.text.x  = element_blank(),
                 axis.ticks.x = element_blank())

# recombina com legenda coletada
p_all <- ((p1 | p2) / (p3 | p4)) 

p_all

ggsave("KM_4covariates.pdf",
       p_all, width = 16, height = 12, dpi = 300, bg = "white")

# Cox model

cox <- coxph(Surv(tempo, status) ~ x1+x2+x3+x4, data = df, x=TRUE)
summary(cox)
# Schoenfeld residuals
zph <- cox.zph(cox, transform = "identity")


p1 <- ggcoxzph(
  zph,
  font.main   = 12,
  point.size  = 2,
  point.shape = 19,
  point.alpha = 1,
  point.col   = "black",
  xlab = "Time (months)",
  #ylab = "Standardized Schoenfeld residuals",
  ggtheme = theme_bw(base_size = 18)+
    theme(
      legend.position = "right",
      strip.text   = element_text(size = 18),
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 18)
    )
)



for (i in seq_along(p1)) {
  for (j in seq_along(p1[[i]]$layers)) {
    if (inherits(p1[[i]]$layers[[j]]$geom, "GeomLine")) {
      p1[[i]]$layers[[j]]$aes_params$linewidth <- 1  # ggplot2 novo
      p1[[i]]$layers[[j]]$aes_params$size      <- 1  # fallback
    }
  }
}

p1


g <- arrangeGrob(grobs = p1, ncol = 2)  # ajuste ncol
ggsave("Schoenfeld_residuals.pdf", g, width = 12, height = 12, dpi = 300, bg = "white")



p1[[4]]

ggsave("Schoenfeld_residuals_x4.pdf", p1[[4]], width = 6, height = 6, dpi = 300, bg = "white")




#Estimated survival rates


surv_table_by_cov_wide <- function(df, cov, times = c(12, 24, 36, 48, 60), digits = 1) {
  
  fml <- as.formula(paste0("Surv(tempo, status) ~ ", cov))
  fit <- survfit(fml, data = df)
  
  sfit <- summary(fit, times = times)
  
  long <- data.frame(
    Covariate = cov,
    Level     = gsub(paste0("^", cov, "="), "", sfit$strata),
    Time      = sfit$time,
    surv      = sfit$surv,
    low       = sfit$lower,
    high      = sfit$upper
  ) %>%
    mutate(
      surv = round(100 * surv, digits),
      low  = round(100 * low,  digits),
      high = round(100 * high, digits),
      `S(t) (95% CI)` = paste0(surv, "% (", low, "–", high, "%)"),
      Time = paste0("t=", Time)
    ) %>%
    select(Covariate, Level, Time, `S(t) (95% CI)`)
  
  wide <- long %>%
    pivot_wider(names_from = Time, values_from = `S(t) (95% CI)`) %>%
    arrange(Covariate, Level)
  
  return(wide)
}


times_eval <- c(0,12, 24, 36, 48, 60)

tab_x1 <- surv_table_by_cov_wide(df, "x1", times = times_eval)
tab_x2 <- surv_table_by_cov_wide(df, "x2", times = times_eval)
tab_x3 <- surv_table_by_cov_wide(df, "x3", times = times_eval)
tab_x4 <- surv_table_by_cov_wide(df, "x4", times = times_eval)

tab_all <- bind_rows(tab_x1, tab_x2, tab_x3, tab_x4)

tab_all




gt_tab <- tab_all %>%
  gt() %>%
  tab_header(
    title = "Kaplan–Meier survival probabilities at selected time points",
    subtitle = "Values are S(t) with 95% confidence intervals"
  ) %>%
  cols_label(
    Covariate = "Covariate",
    Level     = "Level",
    `t=0` = "0 month",
    `t=12` = "12 months",
    `t=24` = "24 months",
    `t=36` = "36 months",
    `t=48` = "48 months",
    `t=60` = "60 months"
  )

gt_tab
# Para exportar LaTeX:
gt::as_latex(gt_tab)

##Hazard rate ratio###

#------------------------------------------------------------
# 0) Você já tem:
# emvg <- optim(..., hessian=TRUE)
# emv_par <- c(emvg$par[1:2], exp(emvg$par[3:4]), emvg$par[5])
# e a sua h_t(t, x, par) que usa par = (alpha1, alpha2, lambda, theta, beta)
#------------------------------------------------------------

# 1) HR(t) (escala natural) usando sua h_t
HR_t <- function(t, par_nat){
  h1 <- h_t(t = t, x = 1, par = par_nat)
  h0 <- h_t(t = t, x = 0, par = par_nat)
  h1 / h0
}

# 2) log(HR(t)) como função dos parâmetros NA ESCALA DO optim:
# par_raw = (alpha1, alpha2, loglambda, logtheta, beta)
logHR_raw <- function(t, par_raw){
  par_nat <- c(par_raw[1:2], exp(par_raw[3:4]), par_raw[5])
  log(HR_t(t, par_nat))
}

# 3) gradiente numérico (diferenças centrais) de logHR_raw em relação a par_raw
grad_logHR_raw <- function(t, par_raw, eps = 1e-6){
  p <- length(par_raw)
  gr <- numeric(p)
  
  for(j in 1:p){
    step <- eps * (abs(par_raw[j]) + 1)  # passo escalado
    up <- par_raw; up[j] <- up[j] + step
    dn <- par_raw; dn[j] <- dn[j] - step
    gr[j] <- (logHR_raw(t, up) - logHR_raw(t, dn)) / (2 * step)
  }
  gr
}

# 4) IC95% pontual para HR(t) via delta method (Wald) na escala log
HR_CI_pointwise <- function(t_grid, par_raw_hat, vcov_raw_hat){
  out <- data.frame(time = t_grid, HR = NA_real_, low = NA_real_, high = NA_real_)
  
  for(i in seq_along(t_grid)){
    t <- t_grid[i]
    g  <- logHR_raw(t, par_raw_hat)
    gr <- grad_logHR_raw(t, par_raw_hat)
    se <- sqrt(as.numeric(t(gr) %*% vcov_raw_hat %*% gr))
    
    out$HR[i]   <- exp(g)
    out$low[i]  <- exp(g - 1.96 * se)
    out$high[i] <- exp(g + 1.96 * se)
  }
  out
}


t_grid <- seq(0, 60, by = 0.01)


#####
# X1
#####


data <- dados[dados$tempo > 0, c(1,2,3)]
data <- data.frame(status = data$status,
                   tempo = data$tempo,
                   X1 = data$x1)

data_alpha <- dados[dados$tempo > 0, 3]

emvg <- optim(par=c(-0.1, -0.1, -0.5, -0.5, 0.1), veroGTDLv2, x_alpha = data_alpha, x_beta = data, hessian = TRUE, method = "BFGS")
emvg



par_raw_hat  <- emvg$par
vcov_raw_hat <- solve(emvg$hessian)

ci_hr <- HR_CI_pointwise(t_grid, par_raw_hat, vcov_raw_hat)

head(ci_hr)


hr_1 <- ggplot(ci_hr, aes(x = time)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = "Pointwise 95%CI"), alpha = 0.2) +
  geom_line(aes(y = HR, color = "Pointwise estimate"), linewidth = 1.1) +
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,0.5)) +
  scale_color_manual(values = "black", name = "") +
  scale_fill_manual(values = "grey20", name = "") +
  labs(title="Consultation information",
    x = "Time (months)",
    y = "Hazard rate ratio"
  ) +
  annotate("text", x = 0, y = 0.13, label = "Reference: Without consultation",
            hjust = 0, vjust = 1, size = 6) +
  theme_bw(base_size = 26) +
  theme(legend.position = "none")

hr_1



#####
# X2
#####



data <- dados[dados$tempo > 0, c(1,2,4)]
data <- data.frame(status = data$status,
                   tempo = data$tempo,
                   X21 = ifelse(data$x2 ==1, 1, 0), 
                   X22 = ifelse(data$x2 ==2, 1, 0)) 


data_alpha <- data[, c(3,4)]

emvg <- optim(par=c(-0.1, -0.1, -0.1, -0.5, -0.5, 0.1, 0.1), veroGTDLv2_, x_alpha = data_alpha, x_beta = data, hessian = TRUE, method = "BFGS")
emvg


par_raw_hat    <- emvg$par
vcov_raw_hat   <- solve(emvg$hessian)

idx2 <- c(1, 2, 4, 5, 6)  # alpha1 alpha2 loglambda logtheta beta1
idx3 <- c(1, 3, 4, 5, 7)  # alpha1 alpha3 loglambda logtheta beta2

vcov_raw_hat_level2 <- vcov_raw_hat[idx2, idx2]
vcov_raw_hat_level3 <- vcov_raw_hat[idx3, idx3]

# IMPORTANT: mantenha TUDO na escala raw (sem exp aqui)
par_level2_raw <- c(par_raw_hat[1], par_raw_hat[2], par_raw_hat[4], par_raw_hat[5], par_raw_hat[6])
par_level3_raw <- c(par_raw_hat[1], par_raw_hat[3], par_raw_hat[4], par_raw_hat[5], par_raw_hat[7])


ci_hr_level2 <- HR_CI_pointwise(t_grid, par_level2_raw, vcov_raw_hat_level2)
ci_hr_level2$contrast <- "3 to 5 years"

ci_hr_level3 <- HR_CI_pointwise(t_grid, par_level3_raw, vcov_raw_hat_level3)
ci_hr_level3$contrast <- "Above 5 years"

ci_hr_both <- rbind(ci_hr_level2, ci_hr_level3)



hr_2 <- ggplot(ci_hr_both,
               aes(x = time, y = HR,
                  # linetype = contrast,
                   color = contrast,
                   fill = contrast)) +
  
  geom_ribbon(aes(ymin = low, ymax = high),
              alpha = 0.15, color = NA) +
  
  geom_line(linewidth = 1.1) +
  
  scale_color_manual(values = c(
    "3 to 5 years" = "black",
    "Above 5 years" = "#1f77b4"
  )) +
  
  scale_fill_manual(values = c(
    "3 to 5 years" = "black",
    "Above 5 years" = "#1f77b4"
  )) +
  
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,0.5)) +
  
  labs(
    title = "Time of employment",
    x = "Time (months)",
    y = "",
    color = "",
   # linetype = "",
    fill = ""
  ) +
  
  annotate("text", x = 0, y = 0.14,
           label = "Reference: Up to 3 years",
           hjust = 0, vjust = 1, size = 6) +
  
  theme_bw(base_size = 26) +
  theme(
    legend.position = c(0.4, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_text(size = 26),
    legend.text  = element_text(size = 26)
  )

hr_2


#####
# X3
#####


data <- dados[dados$tempo > 0, c(1,2,5)]
data <- data.frame(status = data$status,
                   tempo = data$tempo,
                   X1 = data$x3)

data_alpha <- dados[dados$tempo > 0, 5]

emvg <- optim(par=c(-0.1, -0.1, -0.5, -0.5, 0.1), veroGTDLv2, x_alpha = data_alpha, x_beta = data, hessian = TRUE, method = "BFGS")
emvg

par_raw_hat  <- emvg$par
vcov_raw_hat <- solve(emvg$hessian)


ci_hr <- HR_CI_pointwise(t_grid, par_raw_hat, vcov_raw_hat)

head(ci_hr)


hr_3 <- ggplot(ci_hr, aes(x = time)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = "Pointwise 95%CI"), alpha = 0.2) +
  geom_line(aes(y = HR, color = "Pointwise estimate"), linewidth = 1.1) +
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_y_continuous(limits = c(0,3.5), breaks = seq(0,3.5,0.5)) +
  scale_color_manual(values = "black", name = "") +
  scale_fill_manual(values = "grey20", name = "") +
  labs(title="Previously Settled Debts",
       x = "Time (months)",
       y = "Hazard rate ratio"
  ) +
  annotate("text", x = 0, y = 0.1, label = "Reference: Does not have",
           hjust = 0, vjust = 1, size = 6) +
  theme_bw(base_size = 26) +
  theme(legend.position = "bottom")

hr_3



####
# X4
####

data <- dados[dados$tempo > 0, c(1,2,6)]
data <- data.frame(status = data$status,
                   tempo = data$tempo,
                   X1 = data$x4)

data_alpha <- dados[dados$tempo > 0, 6]

emvg <- optim(par=c(-0.1, -0.1, -0.5, -0.5, 0.1), veroGTDLv2, x_alpha = data_alpha, x_beta = data, hessian = TRUE, method = "BFGS")
emvg



par_raw_hat  <- emvg$par
vcov_raw_hat <- solve(emvg$hessian)

ci_hr <- HR_CI_pointwise(t_grid, par_raw_hat, vcov_raw_hat)

head(ci_hr)


hr_4 <- ggplot(ci_hr, aes(x = time)) +
  geom_ribbon(aes(ymin = low, ymax = high, fill = "Pointwise 95%CI"), alpha = 0.2) +
  geom_line(aes(y = HR, color = "Pointwise estimate"), linewidth = 1.1) +
  scale_x_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_y_continuous(limits = c(0,0.8), breaks = seq(0,1,0.2)) +
  scale_color_manual(values = "black", name = "") +
  scale_fill_manual(values = "grey20", name = "") +
  labs(title="Type of Debt",
    x = "Time (months)",
    y = ""
  ) +
  theme_bw(base_size = 26) +
  annotate("text", x = 0, y = 0.025, label = "Reference: Bank",
           hjust = 0, vjust = 1, size =6) +
  theme(legend.position = "bottom")

hr_4

#ggsave("HR_x4.pdf",hr_4, width = 8, height = 6, dpi = 300, bg = "white")


p1 <- hr_1 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

p2 <- hr_2 + theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

p4 <- hr_4 + theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

p4
p1
p2



p_all <- ((hr_1 | hr_2) / (hr_3 | hr_4))
p_all


p_all

ggsave("HRs_4covariates.pdf",plot=p_all,width = 20, height = 17, dpi = 300,bg = "white")



#####################
#Multivariable model#
#####################

veroGTDLv2 <- function(par, x_alpha, x_beta) {
  
  # parameters
  alph <- par[1:2]
  lambd <- exp(par[3])
  thet <- exp(par[4])
  bet <- par[5:length(par)]
  
  # covariates
  X_beta <- as.matrix(x_beta[x_beta$tempo > 0, 3:ncol(x_beta)])
  xalpha <- cbind(rep(1, nrow(as.matrix(data_alpha))), x_alpha)
  cens <- x_beta[x_beta$tempo > 0, 1]
  tempo <- x_beta[x_beta$tempo > 0, 2]
  
  # model
  aux2 <- log(lambd)*sum(cens) +
    sum( cens*((xalpha%*%alph)*tempo + X_beta%*%bet)) -
    sum( cens*log( 1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet) ) ) - 
    sum( (cens + (1/thet))*log( 1 + ((thet*lambd/(xalpha%*%alph))*log( (1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet))/(1 + exp(X_beta%*%bet)) )) ) )
  
  return(-aux2)
}


##############################################
# GTDL with frailty - likeligood function
##############################################

veroGTDLv2_ <- function(par, x_alpha, x_beta) {
  
  # parameters
  alph <- par[1:3]
  lambd <- exp(par[4])
  thet <- exp(par[5])
  bet <- par[6:length(par)]
  
  # covariates
  X_beta <- as.matrix(x_beta[x_beta$tempo > 0, 3:ncol(x_beta)])
  xalpha <- as.matrix(cbind(rep(1, nrow(as.matrix(data_alpha))), x_alpha))
  cens <- x_beta[x_beta$tempo > 0, 1]
  tempo <- x_beta[x_beta$tempo > 0, 2]
  
  # model
  aux2 <- log(lambd)*sum(cens) +
    sum( cens*((xalpha%*%alph)*tempo + X_beta%*%bet)) -
    sum( cens*log( 1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet) ) ) - 
    sum( (cens + (1/thet))*log( 1 + ((thet*lambd/(xalpha%*%alph))*log( (1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet))/(1 + exp(X_beta%*%bet)) )) ) )
  
  return(-aux2)
}



# Zero-inflated likelihood


veroZero <- function(x, beta) {
  
  # covariates
  X0 <- as.matrix(x[x$tempo == 0, 3:ncol(x)])
  X <- as.matrix(x[x$tempo > 0, 3:ncol(x)])
  
  # model
  aux1 <- sum((X0%*%beta)) - sum(log( 1 + exp(X0%*%beta))) - sum(log(1 + exp(X%*%beta)))
  return(-aux1)
}

veroGTDLg <- function(x, par) {
  
  # parameters
  alph <- par[1]
  lambd <- exp(par[2])
  thet <- exp(par[3])
  bet <- par[4:length(par)]
  
  # covariates
  X <- as.matrix(x[, 3:ncol(x)])
  cens <- x[,1]
  tempo <- x[,2]
  
  # model
  aux2 <- log(lambd)*sum(cens) +
    sum( cens*(X%*%bet + alph*tempo)) -
    sum( cens*log( 1 + exp(X%*%bet + alph*tempo)) ) - 
    sum( (cens + (1/thet))*log( 1 + ((thet*lambd/alph)*log( (1 + exp(X%*%bet + alph*tempo))/(1 + exp(X%*%bet)) )) ) )
  
  return(-aux2)
}


# GTDL Gamma likelihood

veroGTDLv2 <- function(par, x_alpha, x_beta) {
  
  # parameters
  alph <- par[1:6]
  lambd <- exp(par[7])
  thet <- exp(par[8])
  bet <- par[9:length(par)]
  
  # covariates
  X_beta <- as.matrix(x_beta[x_beta$tempo > 0, 3:ncol(x_beta)])
  xalpha <- cbind(rep(1, nrow(as.matrix(data_alpha))), as.matrix(x_alpha))
  cens <- x_beta[x_beta$tempo > 0, 1]
  tempo <- x_beta[x_beta$tempo > 0, 2]
  
  # model
  aux2 <- log(lambd)*sum(cens) +
    sum( cens*((xalpha%*%alph)*tempo + X_beta%*%bet)) -
    sum( cens*log( 1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet) ) ) - 
    sum( (cens + (1/thet))*log( 1 + ((thet*lambd/(xalpha%*%alph))*log( (1 + exp((xalpha%*%alph)*tempo + X_beta%*%bet))/(1 + exp(X_beta%*%bet)) )) ) )
  
  return(-aux2)
}






dados <- read.table(file = "dados_censura_versao4.txt", header = TRUE)

# etapa 1 
# todas as variáveis 

head(dados)
summary(dados$x1)
summary(dados$x2)

dados$x2_1 <- ifelse(dados$x2==1, 1, 0) 
dados$x2_2 <- ifelse(dados$x2==2, 1, 0) 

head(dados)

#data <- dados[dados$tempo > 0, c(1,2,3)]
#data <- dados[dados$tempo > 0, c(1,2,3,5)]
#data <- dados[dados$tempo > 0, c(1,2,3,5,6)]
data <- dados[dados$tempo > 0, c(1,2,3,5,6,7,8)]

head(data)

#p0 <- c(-0.1, -0.5, -0.5, 0.1)
#p0 <- c(-0.1380946, -0.8513506, -0.4538198, -0.1316308, 0.1)
#p0 <- c(-0.13817337, -0.85890109, -0.44937695, -0.14602869,  0.08015024, 0.1)
p0 <- c(-0.14021270, -0.69345937, -0.93122346, -0.05761031,  0.08064615, -0.92491018, 0.1, 0.1)

emvg <- optim(par=p0, veroGTDLg, x=data, hessian = TRUE)
emvg

c(emvg$par[1], exp(emvg$par[2]), exp(emvg$par[3]), emvg$par[4:5])
dp <- sqrt(diag(solve(emvg$hessian)))
dp[2] <- deltamethod(g = ~ exp(x2), mean = emvg$par, cov = solve(emvg$hessian))
dp[3] <- deltamethod(g = ~ exp(x3), mean = emvg$par, cov = solve(emvg$hessian))

emvg$par - 1.96*dp
emvg$par + 1.96*dp

2 * emvg$value + 2 * length(emvg$par)
2 * emvg$value + log(nrow(data)) * length(emvg$par)

### Etapa 2 

data <- dados[dados$tempo > 0, c(1,2,3,7,8,5,6)]

head(data)

p0 <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, -0.70758009, -0.91478664, -0.0596565, 0.08554532, 0.08922454, 0.08411799, -0.93193246)


data <- dados[dados$tempo > 0, c(1,2,3,7,8,5,6)]
data_alpha <- dados[dados$tempo > 0, c(3,7,8,5,6)]

emvg <- optim(par=p0, veroGTDLv2, x_alpha = data_alpha, x_beta = data, hessian = TRUE, method = "BFGS")
emvg

2 * emvg$value + 2 * length(emvg$par)
2 * emvg$value + log(nrow(data)) * length(emvg$par)


emv <- c(emvg$par[1:6], exp(emvg$par[7]), exp(emvg$par[8]), emvg$par[9:13])
round(emv, 4)

dp <- sqrt(diag(solve(emvg$hessian)))
round(dp, 4)

dp[7] <- deltamethod(g = ~ exp(x7), mean = emvg$par, cov = solve(emvg$hessian))
dp[8] <- deltamethod(g = ~ exp(x8), mean = emvg$par, cov = solve(emvg$hessian))

round(emv - 1.96*dp, 4)
round(emv + 1.96*dp, 4)

2 * emvg$value + 2 * length(emvg$par)
2 * emvg$value + log(nrow(data)) * length(emvg$par)

### Etapa 3  

data <- dados[,c(1,2,3,4,5,6,7,8)]
#data0 <- data.frame(status = data$status,
#                    tempo = data$tempo,
#                    interc = rep(1,nrow(data)),
#                    X1 = data$x1)

#data0 <- data.frame(status = data$status,
#                    tempo = data$tempo,
#                    interc = rep(1,nrow(data)),
#                    X1 = data$x1,
#                    X3 = data$x3)

#data0 <- data.frame(status = data$status,
#                    tempo = data$tempo,
#                    interc = rep(1,nrow(data)),
#                    X1 = data$x1,
#                    X3 = data$x3,
#                    X4 = data$x4)

data0 <- data.frame(status = data$status,
                    tempo = data$tempo,
                    interc = rep(1,nrow(data)),
                    X1 = data$x1,
                    X21 = data$x2_1,
                    X22 = data$x2_2,
                    X3 = data$x3,
                    X4 = data$x4)

head(data0)

#p0 <- c(0.5, 0.5)
#p0 <- c(-0.8545611, -0.321457, 0.1)
#p0 <- c(-0.86090645, -0.33076808,  0.04252194, 0.1)
p0 <- c(-0.72829667, -0.32738892, 0.1, 0.1,  0.04340229, -0.30193754)
emv0 <- optim(par=p0, veroZero, x=data0, hessian = TRUE)
round(emv0$par, 4)

dp0 <- sqrt(diag(solve(emv0$hessian)))
round(dp0, 4)

round(emv0$par - 1.96*dp0, 4)
round(emv0$par + 1.96*dp0, 4)

2 * emv0$value + 2 * length(emv0$par) # AIC 
2 * emv0$value + log(nrow(data0)) * length(emv0$par) # BIC


# -log-verossimilhança, AIC e BIC do modelo completo 
emvg$value + emv0$value # -log-verossimilhança
2 * emvg$value + 2 * length(emvg$par) + 2 * emv0$value + 2 * length(emv0$par) # AIC
2 * emvg$value + log(nrow(data)) * length(emvg$par) + 2 * emv0$value + log(nrow(data0)) * length(emv0$par) # BIC










sobrevGTDLg <- function(t,x, par, par0) {
  
  alpha <- par[1:2]
  lambda <- par[3]
  theta <- par[4]
  beta <- par[5]
  
  gama0 <- par0[1]
  gama1 <- par0[2]
  p0 <- exp(gama0+gama1*x)/(1+exp(gama0+gama1*x))
  
  st <- (1-p0)*(1 + ( (lambda*theta/(alpha[1]+alpha[2]))*log( (1+ exp(x*beta + (alpha[1] + alpha[2]*x)*t))/(1 + exp(x*beta))) ) )^(-1/theta)
  return(st)
}


sobrevGTDLg_ <- function(t,x, par, par0) {
  
  alpha <- par[1:3]
  lambda <- par[4]
  theta <- par[5]
  beta <- par[6:7]
  
  
  gama0 <- par0[1]
  gama1 <- par0[2]
  gama2 <- par0[3]
  p0 <- exp(gama0+gama1*x[1] +gama2*x[2])/(1+exp(gama0+gama1*x[1]+ gama2*x[2]))
  
  st <- (1-p0)*(1 + ( (lambda*theta/(alpha[1]+alpha[2]*x[1]+alpha[3]*x[2]))*log( (1+ exp(x[1]*beta[1] + x[2]*beta[2] + (alpha[1] + alpha[2]*x[1] + alpha[3]*x[2])*t))/(1 + exp(x[1]*beta[1]+x[2]*beta[2]))) ) )^(-1/theta)
  return(st)
}







#Estimated p0 and p1 according to the covariate combination


gamma_hat <- emv0$par
V_gamma   <- solve(emv0$hessian)
names(gamma_hat) <- c("interc","X1","X21","X22","X3","X4")
dimnames(V_gamma) <- list(names(gamma_hat), names(gamma_hat))


par_surv_hat <- emvg$par
V_surv_raw   <- solve(emvg$hessian)

alpha_hat <- par_surv_hat[1:6]
names(alpha_hat) <- c("interc","X1","X21","X22","X3","X4")

loglambda_hat <- par_surv_hat[7]
logtheta_hat  <- par_surv_hat[8]

beta_hat <- par_surv_hat[9:13]
names(beta_hat) <- c("X1","X21","X22","X3","X4")

surv_hat <- c(alpha_hat, loglambda_hat, logtheta_hat, beta_hat)
names(surv_hat) <- c(paste0("a_", names(alpha_hat)), "loglambda", "logtheta", paste0("b_", names(beta_hat)))


idx_surv <- c(1:6, 7, 8, 9:13)
V_surv <- V_surv_raw[idx_surv, idx_surv]
dimnames(V_surv) <- list(names(surv_hat), names(surv_hat))

# ===== logit function =====
expit <- function(z) 1/(1+exp(-z))

# ===== p0(x) =====
p0_of_x <- function(x_row, gamma) {
  # x_row: vetor nomeado com interc, X1, X21, X22, X3, X4
  expit(sum(gamma * x_row))
}

# ===== p1(x) fração de cura =====
p1_of_x <- function(x_row, gamma, survpar) {
  # x_row: interc, X1, X21, X22, X3, X4
  # survpar: a_* (inclui intercepto), loglambda, logtheta, b_* (sem intercepto)
  
  # p0
  p0 <- p0_of_x(x_row, gamma)
  
  # alpha(x)
  a_vec <- survpar[paste0("a_", names(x_row))]
  alpha_x <- sum(a_vec * x_row)
  
  
  b_vec <- survpar[paste0("b_", names(beta_hat))]
  
  x_beta <- x_row[names(beta_hat)]
  eta <- sum(b_vec * x_beta)
  
  lambda <- exp(survpar["loglambda"])
  theta  <- exp(survpar["logtheta"])
  
  
  A <- 1 - (lambda*theta/alpha_x) * log(1 + exp(eta))
  
  
  if(!is.finite(A) || A <= 0 || !is.finite(alpha_x)) return(NA_real_)
  
  p_star <- A^(-1/theta)
  p1 <- (1 - p0) * p_star
  return(as.numeric(p1))
}

# ===== delta method: p0 e p1 =====
delta_ci <- function(g_fun, psi_hat, V_psi, level = 0.95) {
  g_hat <- g_fun(psi_hat)
  if(is.na(g_hat)) return(c(est=NA, lo=NA, hi=NA, se=NA))
  
  grad_g <- numDeriv::grad(func = g_fun, x = psi_hat)
  var_g  <- as.numeric(t(grad_g) %*% V_psi %*% grad_g)
  se_g   <- sqrt(max(var_g, 0))
  
  z <- qnorm(1 - (1-level)/2)
  c(est = g_hat, lo = g_hat - z*se_g, hi = g_hat + z*se_g, se = se_g)
}


grid <- expand.grid(
  X1 = c(0,1),
  X2 = c(0,1,2),
  X3 = c(0,1),
  X4 = c(0,1)
)

grid$X21 <- as.integer(grid$X2 == 1)
grid$X22 <- as.integer(grid$X2 == 2)
grid$interc <- 1


xmat <- as.matrix(grid[, c("interc","X1","X21","X22","X3","X4")])


psi_hat <- c(gamma_hat, surv_hat)
V_psi   <- as.matrix(Matrix::bdiag(V_gamma, V_surv))
dimnames(V_psi) <- list(names(psi_hat), names(psi_hat))


res <- lapply(1:nrow(xmat), function(i) {
  x_row <- xmat[i,]
  names(x_row) <- colnames(xmat)
  
  
  g_p0 <- function(psi) {
    gamma <- psi[names(gamma_hat)]
    p0_of_x(x_row, gamma)
  }
  
  # g(psi) para p1: depende de gamma e survpar
  g_p1 <- function(psi) {
    gamma <- psi[names(gamma_hat)]
    surv  <- psi[names(surv_hat)]
    p1_of_x(x_row, gamma, surv)
  }
  
  ci0 <- delta_ci(g_p0, psi_hat, V_psi)
  ci1 <- delta_ci(g_p1, psi_hat, V_psi)
  
  data.frame(
    X1 = grid$X1[i], X2 = grid$X2[i], X3 = grid$X3[i], X4 = grid$X4[i],
    p0_est = ci0["est"], p0_lo = ci0["lo"], p0_hi = ci0["hi"],
    p1_est = ci1["est"], p1_lo = ci1["lo"], p1_hi = ci1["hi"]
  )
})

out <- do.call(rbind, res)
out[, grep("p0_|p1_", names(out))] <- round(out[, grep("p0_|p1_", names(out))], 3)

estimates_multi<-out[order(out$X2, out$X3, out$X4), ]
estimates_multi

estimates_multi_sorted <- estimates_multi %>%
  arrange(p0_est)

estimates_multi_sorted


#plot

lab_X1 <- c(`0`="Consultation information: 0: No  ($X_1$)",
            `1`="1: With consultation")
lab_X2 <- c(`0`="Time of employment: 0: Up to 3 years ($X_2$)",
            `1`="1: 3 to 5 years",
            `2`="2: More than 5 years")
lab_X3 <- c(`0`="Previously settled debts: 0: None ($X_3$)",
            `1`="1: Yes")
lab_X4 <- c(`0`="Type of debt: 0: Bank ($X_4$)",
            `1`="1: Other segments")


plot_df <- out %>%
  mutate(
    X1_lab = lab_X1[as.character(X1)],
    X2_lab = lab_X2[as.character(X2)],
    X3_lab = lab_X3[as.character(X3)],
    X4_lab = lab_X4[as.character(X4)],
    combo  = paste0("X1=",X1," | X2=",X2," | X3=",X3," | X4=",X4)
  ) %>%
  # ordenar as combinações para ficar “bonito”
  arrange(X2, X3, X4, X1) %>%
  mutate(combo = factor(combo, levels = unique(combo))) %>%
  select(combo, X1, X2, X3, X4, X1_lab, X2_lab, X3_lab, X4_lab,
         p0_est, p0_lo, p0_hi, p1_est, p1_lo, p1_hi) %>%
  pivot_longer(
    cols = c(p0_est, p0_lo, p0_hi, p1_est, p1_lo, p1_hi),
    names_to = c("param","stat"),
    names_pattern = "(p0|p1)_(est|lo|hi)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    param = factor(param, levels = c("p0","p1"),
                   labels = c("Probability of immediate settlement (p0)","Probability of debt non-settlement (p1)"))
  )

caption_txt <- paste0(
  "X1: ", gsub("\\s*\\(\\$X_1\\$\\)", "", lab_X1["0"]), " / ", gsub("\\s*\\(\\$X_1\\$\\)", "", lab_X1["1"]), "; ",
  "X2: ", gsub("\\s*\\(\\$X_2\\$\\)", "", lab_X2["0"]), " / ", lab_X2["1"], " / ", lab_X2["2"], "; ",
  "X3: ", gsub("\\s*\\(\\$X_3\\$\\)", "", lab_X3["0"]), " / ", lab_X3["1"], "; ",
  "X4: ", gsub("\\s*\\(\\$X_4\\$\\)", "", lab_X4["0"]), " / ", lab_X4["1"], "."
)


# ---- 3) Figura (dois painéis: p0 e p1) ----
p <- ggplot(plot_df, aes(x = est, y = combo)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = 0.6) +
  geom_point(size = 2.2) +
  facet_wrap(~param, ncol = 2, scales = "free_x") +
  coord_cartesian(xlim = c(0, 0.4)) +
  labs(
    x = "Estimate (95% CI)",
    y = "Covariate combination",
    caption = caption_txt
  ) +
  theme_bw(base_size = 18) +
  theme(
    strip.text = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    plot.caption = element_text(size = 10, hjust = 2)
  )


print(p)


p1 <- ggplot(plot_df, aes(x = est, y = combo)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), width = 0.5, linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~param, ncol = 2, scales = "free_x") +
  ggh4x::facetted_pos_scales(
    x = list(
      param == "Probability of immediate settlement (p0)" ~ scale_x_continuous(limits = c(0.15, 0.45), breaks = seq(0.15, 0.45, 0.05)),
      param == "Probability of debt non-settlement (p1)" ~ scale_x_continuous(limits = c(0.05, 0.35),breaks =seq(0.05,0.35,0.05))
    )
  ) +
  labs(x = "Estimate (95% CI)", y = "Covariate combination", caption = caption_txt) +
  theme_bw(base_size = 26) +
  theme(
    strip.text = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    plot.caption = element_text(size = 15, hjust =1)
  )

print(p1)



ggsave("Estimated_p0_p1_ZI_GTDL_model.pdf", p1, width = 22, height = 15, dpi = 300)




combo_levels_p0 <- plot_df %>%
  filter(param == "Probability of immediate settlement (p0)") %>%
  arrange(est) %>%                      # est = p0_est (já está no formato long)
  pull(combo) %>%
  unique()

plot_df2 <- plot_df %>%
  mutate(combo = factor(combo, levels = combo_levels_p0))

# 2) Figura (mesma que você mandou), agora com a ordem do eixo-y baseada em p0
p2 <- ggplot(plot_df2, aes(x = est, y = combo)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.5, linewidth = 0.9) +
  geom_point(size = 3) +
  facet_wrap(~param, ncol = 2, scales = "free_x") +
  ggh4x::facetted_pos_scales(
    x = list(
      param == "Probability of immediate settlement (p0)" ~
        scale_x_continuous(limits = c(0.15, 0.45), breaks = seq(0.15, 0.45, 0.05)),
      param == "Probability of debt non-settlement (p1)" ~
        scale_x_continuous(limits = c(0.05, 0.35), breaks = seq(0.05, 0.35, 0.05))
    )
  ) +
  labs(x = "Estimate (95% CI)", y = "Covariate combination", caption = caption_txt) +
  theme_bw(base_size = 26) +
  theme(
    strip.text = element_text(size = 26),
    axis.text.y = element_text(size = 26),
    plot.caption = element_text(size = 15, hjust = 1)
  )

print(p2)

ggsave("Estimated_p0_p1_ZI_GTDL_model_v2.pdf", p2, width = 22, height = 15, dpi = 300)






###Comparison with other model###

## Defective Gompertz model

source("functions_defective_Gompertz_model.R")

#data structure is different in defective Gompertz approach

df <- df %>%
  mutate(delta0 = ifelse(tempo > 0, 1, 0),
         delta=status,
         s=tempo)


par(mar=c(5,4,1,1)+.1)
km<-survfit(Surv(tempo,status)~x1,data=df)
plot(km,lty=c(1,2),lwd=c(2),col=c(1),xlab=" Time (weeks) ", ylab="Survival function",conf.int=FALSE,mark.time=F,cex.axis=1.4,cex.lab=1.4, bty="n",ylim=c(0,1))



attach(df)
n<-dim(df)




###
#X1
###


x = as.numeric(x1) - 1
v1 = try(mle(log.like_cov2,
             start  = list(a=-1,b0=-1,b1=0.1,b2=1,b3=1),
             method = "BFGS", T))

est1 = coef(v1)

aux1  = c(0,1)
p0x   = px(coef(v1)[2], coef(v1)[3], aux1)
b0x   = bx(coef(v1)[4], coef(v1)[5], aux1)
curax = p1x(coef(v1)[1], p0x, b0x, aux1)

Estimativa = round(c(est1, p0x, curax), 3)
EP         = ep_g(coef(v1), vcov(v1))
ic         = ic_g(Estimativa, EP)

print(data.frame(Estimativa, EP, ic))
print(c(logLik(v1), AIC(v1), AIC(v1, k = log(n))))

## =========================
## 2) Kaplan–Meier (df longo p/ ggplot)
## =========================
km <- survfit(Surv(s, delta) ~ x1, data = df)

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  covar  = c(rep(0, km$strata[1]), rep(1, km$strata[2])),
  metodo = "Kaplan-Meier"
)


t_grid <- seq(0.001, max(df$s, na.rm = TRUE), by = 0.01)

sp <- sg_cov2(est1[1], p0x, b0x, t_grid)  

df_mod <- data.frame(
  tempo   = t_grid,
  sobre_0 = sp[,1],
  sobre_1 = sp[,2]
)

df_mod_long <- pivot_longer(
  df_mod,
  cols      = c(sobre_0, sobre_1),
  names_to  = "covar",
  values_to = "surv"
) %>%
  mutate(
    covar  = ifelse(covar == "sobre_0", 0, 1),
    metodo = "Parametric"
  )


df_all <- bind_rows(km_df, df_mod_long)

p1<-ggplot(df_all, aes(x = tempo, y = surv)) +
  
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  geom_line(
    data = subset(df_all, metodo == "Parametric"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  labs(
    x = "Time (month)",
    y = "S(t)",
    color = "Consultation information",
    linetype = "Estimation method"
  ) +
  
  scale_color_manual(
    labels = c("without information", "with information"),
    values = c("red", "blue")
  ) +
  
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "solid",
               "Parametric"  = "dashed")
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, ceiling(max(df_all$tempo, na.rm=TRUE)), 5)) +
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_rect(color = "white"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )


p1





####
##X2
####


df$x2 <- factor(df$x2)
df$x2 <- relevel(df$x2, ref = "0")

x2_dum  <- model.matrix(~ x2 - 1, data = df)  
x2_dum2 <- x2_dum[, -1, drop = FALSE]         


cov1 = x2_dum[, -2][, 1]
cov2 = x2_dum[, -2][, 2]


v1 = try(mle(log.like_cov3,
             start  = list(a=-1,b0=-1,b1=0.1,b2=1,b3=1,c1=0,c2=0),
             method = "BFGS", T))

est1 = coef(v1)

p0x_0 = px_2(coef(v1)[2], coef(v1)[3], coef(v1)[6], 0, 0)
p0x_1 = px_2(coef(v1)[2], coef(v1)[3], coef(v1)[6], 1, 0)
p0x_2 = px_2(coef(v1)[2], coef(v1)[3], coef(v1)[6], 0, 1)

b0x_0 = bx_2(coef(v1)[4], coef(v1)[5], coef(v1)[7], 0, 0)
b0x_1 = bx_2(coef(v1)[4], coef(v1)[5], coef(v1)[7], 1, 0)
b0x_2 = bx_2(coef(v1)[4], coef(v1)[5], coef(v1)[7], 0, 1)

curax_0 = p1x_2(coef(v1)[1], p0x_0, b0x_0)
curax_1 = p1x_2(coef(v1)[1], p0x_1, b0x_1)
curax_2 = p1x_2(coef(v1)[1], p0x_2, b0x_2)

Estimativa = round(c(est1, p0x_0, p0x_1, p0x_2, curax_0, curax_1, curax_2), 3)
EP         = ep_g_2(coef(v1), vcov(v1))
ic         = ic_g(Estimativa, EP)

print(data.frame(Estimativa, EP, ic))
print(c(logLik(v1), AIC(v1), AIC(v1, k = log(n))))

km <- survfit(Surv(s, delta) ~ x2, data = df)


levs <- levels(df$x2)
n_str <- as.integer(km$strata)  # tamanhos por estrato, na ordem do survfit

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  covar  = rep(levs, times = n_str),
  metodo = "Kaplan-Meier"
)

t_grid <- seq(0.001, max(df$s, na.rm = TRUE), by = 0.01)

sp <- sg_cov2_2(est1[1],
                p0x_0, p0x_1, p0x_2,
                b0x_0, b0x_1, b0x_2,
                t_grid)  

df_mod <- data.frame(
  tempo   = t_grid,
  sobre_0 = sp[,1],
  sobre_1 = sp[,2],
  sobre_2 = sp[,3]
)

df_mod_long <- pivot_longer(
  df_mod,
  cols      = c(sobre_0, sobre_1, sobre_2),
  names_to  = "covar",
  values_to = "surv"
) %>%
  mutate(
    covar  = recode(covar,
                    "sobre_0" = levs[1],
                    "sobre_1" = levs[2],
                    "sobre_2" = levs[3]),
    metodo = "Parametric"
  )


df_all <- bind_rows(km_df, df_mod_long)

p2<-ggplot(df_all, aes(x = tempo, y = surv)) +
  
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  geom_line(
    data = subset(df_all, metodo == "Parametric"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  labs(
    x = "Time (month)",
    y = "S(t)",
    color = "Time of employment",
    linetype = "Estimation method"
  ) +
  
  scale_color_manual(
    breaks = levs,
    labels = c("Up to 3 years", "From 3 to 5 years", "Above 5 years"),
    values = c("red", "blue", "green")
  ) +
  
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "solid",
               "Parametric"  = "dashed")
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, ceiling(max(df_all$tempo, na.rm = TRUE)), 5)) +
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_rect(color = "white"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

p2


####
# X3 
####

x = as.numeric(x3) - 1
v1 = try(mle(log.like_cov2,
             start  = list(a=-1,b0=-1,b1=0.1,b2=1,b3=1),
             method = "BFGS", T))

est1 = coef(v1)

aux1  = c(0,1)
p0x   = px(coef(v1)[2], coef(v1)[3], aux1)
b0x   = bx(coef(v1)[4], coef(v1)[5], aux1)
curax = p1x(coef(v1)[1], p0x, b0x, aux1)

Estimativa = round(c(est1, p0x, curax), 3)
EP         = ep_g(coef(v1), vcov(v1))
ic         = ic_g(Estimativa, EP)

print(data.frame(Estimativa, EP, ic))
print(c(logLik(v1), AIC(v1), AIC(v1, k = log(n))))


t_grid <- seq(0.001, max(df$s, na.rm = TRUE), by = 0.01)

sp <- sg_cov2(est1[1], p0x, b0x, t_grid)

df_mod <- data.frame(
  tempo   = t_grid,
  sobre_0 = sp[,1],
  sobre_1 = sp[,2]
)

df_mod_long <- pivot_longer(
  df_mod,
  cols      = c(sobre_0, sobre_1),
  names_to  = "covar",
  values_to = "surv"
) %>%
  mutate(
    covar  = ifelse(covar == "sobre_0", 0, 1),
    metodo = "Parametric"
  )


km <- survfit(Surv(s, delta) ~ x3, data = df)

km_df <- data.frame(
  tempo = km$time,
  surv  = km$surv,
  covar = c(rep(0, km$strata[1]), rep(1, km$strata[2])),
  metodo = "Kaplan-Meier"
)




df_all <- bind_rows(km_df, df_mod_long)

p3<-ggplot(df_all, aes(x = tempo, y = surv)) +
  
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  geom_line(
    data = subset(df_all, metodo == "Parametric"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  labs(
    x = "Time (month)",
    y = "S(t)",
    color = "Debt",
    linetype = "Estimation method"
  ) +
  
  scale_color_manual(
    labels = c("Does not have", "It has"),
    values = c("red", "blue")
  ) +
  
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "solid",
               "Parametric"  = "dashed")
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, ceiling(max(df_all$tempo)), 5)) +
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_rect(color = "white"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )


p3


###
#X4
###

x = as.numeric(x4) - 1
v1 = try(mle(log.like_cov2,
             start  = list(a=-1,b0=-1,b1=0.1,b2=1,b3=1),
             method = "BFGS", T))

est1 = coef(v1)

aux1  = c(0,1)
p0x   = px(coef(v1)[2], coef(v1)[3], aux1)
b0x   = bx(coef(v1)[4], coef(v1)[5], aux1)
curax = p1x(coef(v1)[1], p0x, b0x, aux1)

Estimativa = round(c(est1, p0x, curax), 3)
EP         = ep_g(coef(v1), vcov(v1))
ic         = ic_g(Estimativa, EP)

print(data.frame(Estimativa, EP, ic))
print(c(logLik(v1), AIC(v1), AIC(v1, k = log(n))))


km <- survfit(Surv(s, delta) ~ x4, data = df)

km_df <- data.frame(
  tempo  = km$time,
  surv   = km$surv,
  covar  = c(rep(0, km$strata[1]), rep(1, km$strata[2])),
  metodo = "Kaplan-Meier"
)


t_grid <- seq(0.001, max(df$s, na.rm = TRUE), by = 0.01)

sp <- sg_cov2(est1[1], p0x, b0x, t_grid)  # matriz 2 colunas: x=0 e x=1

df_mod <- data.frame(
  tempo   = t_grid,
  sobre_0 = sp[,1],
  sobre_1 = sp[,2]
)

df_mod_long <- pivot_longer(
  df_mod,
  cols      = c(sobre_0, sobre_1),
  names_to  = "covar",
  values_to = "surv"
) %>%
  mutate(
    covar  = ifelse(covar == "sobre_0", 0, 1),
    metodo = "Parametric"
  )


df_all <- bind_rows(km_df, df_mod_long)

p4<-ggplot(df_all, aes(x = tempo, y = surv)) +
  
  geom_step(
    data = subset(df_all, metodo == "Kaplan-Meier"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  geom_line(
    data = subset(df_all, metodo == "Parametric"),
    aes(color = as.factor(covar), linetype = metodo),
    linewidth = 1
  ) +
  
  labs(
    x = "Time (month)",
    y = "S(t)",
    color = "Type of debt",
    linetype = "Estimation method"
  ) +
  
  scale_color_manual(
    labels = c("Bank", "Other segments"),
    values = c("red", "blue")
  ) +
  
  scale_linetype_manual(
    values = c("Kaplan-Meier" = "solid",
               "Parametric"  = "dashed")
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, ceiling(max(df_all$tempo, na.rm=TRUE)), 5)) +
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_rect(color = "white"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )


p4



size_font=20
plot_final <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title = "Zero-adjusted defective Gompertz model") &
  theme(
    plot.title   = element_text(size = size_font, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = size_font),
    axis.title.y = element_text(size = size_font),
    axis.text.x  = element_text(size = size_font),
    axis.text.y  = element_text(size = size_font),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 18)
   # strip.text   = element_text(size = size_font)
  )

plot_final

ggsave("Fitted_ZI_defective_Gompertz.pdf", plot_final, width = 16, height = 12, dpi = 300, bg = "white")


##########################################
# Multivariable Defective Gompertz model #
##########################################

cova <- cbind(
  1,
  as.numeric(x1) - 1,
  x2_dum[, -2, drop = FALSE],   # remove a coluna de referência (ajuste conforme seu ref)
  as.numeric(x3) - 1,
  as.numeric(x4) - 1
)

beta <- rep(1, ncol(cova))
gama <- rep(1, ncol(cova))
par  <- c(-1, beta, gama)

fit_mult <- optim(par = par, fn = log.like_gomp, method = "BFGS", hessian = TRUE)

summary.optim(fit_mult)

AIC <- 2 * fit_mult$value + 2 * length(fit_mult$par)  # ou length(par); prefiro usar o estimado
AIC


##With covariate in a(x)=alpha*x


cova <- cbind(
  1,
  as.numeric(x1) - 1,
  x2_dum[, -2, drop = FALSE],   # remove a coluna de referência (ajuste conforme seu ref)
  as.numeric(x3) - 1,
  as.numeric(x4) - 1
)


alpha <- rep(-0.1, ncol(cova))   # sugestão: começar em 0 (=> a(x)=1)
beta  <- rep(0, ncol(cova))   # sugestão: começar em 0 (=> p0=0.5)
gama  <- rep(0, ncol(cova))   # sugestão: começar em 0 (=> b(x)=1)

par <- c(alpha, beta, gama)

fit_mult <- optim(
  par     = par,
  fn      = log.like_gomp_complete,  # sua log-vero atualizada p/ incluir a(x)
  method  = "BFGS",
  hessian = TRUE
)


summary.optim_complete(fit_mult)

AIC <- 2 * fit_mult$value + 2 * length(fit_mult$par)  
AIC

