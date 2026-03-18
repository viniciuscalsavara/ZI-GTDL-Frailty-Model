

library(tidyverse)
library(ggh4x)
library(ggplot2)
library(dplyr)
source("functions.R")

#Fixed cure rate#

estimates=c(-1,-0.6,-0.1,0.5,0.4,-0.9)
x=c(0,1)
round(p1f(x,estimates),3)

#zero-inflated

round(p0f(x,estimates[1:2]),3)

# -----------------------------------------------
# Plot - Simulation study (asymptotic properties)
# -----------------------------------------------

sizes <- c(100, 200, 300, 500, 800, 1000, 2000)

bias <- tribble(
  ~param, ~`100`, ~`200`, ~`300`, ~`500`, ~`800`, ~`1000`, ~`2000`,
  "alpha",    0.0144, -0.0003,  0.0030, -0.0010, -0.0003, -0.0014, -0.0002,
  "lambda",   0.0412,  0.0105,  0.0134,  0.0031,  0.0061,  0.0026,  0.0002,
  "theta",    0.0540, -0.0139, -0.0014, -0.0132, -0.0003, -0.0092,  0.0009,
  "beta",    -0.0018,  0.0220, -0.0019,  0.0044, -0.0063,  0.0005,  0.0049,
  "gamma0",  -0.0269, -0.0048,  0.0024, -0.0006, -0.0026,  0.0001,  0.0004,
  "gamma1",  -0.0393, -0.0150, -0.0238, -0.0098, -0.0044, -0.0005, -0.0002,
  "p0(x=1)", -0.0017,  0.0008, -0.0003, -0.0001, -0.0002,  0.0006,  0.0004,
  "p0(x=0)", -0.0003,  0.0012,  0.0020,  0.0008,  0.0001,  0.0005,  0.0003,
  "p1(x=1)", -0.0088, -0.0021, -0.0011, -0.0008, -0.0000,  0.0005, -0.0004,
  "p1(x=0)", -0.0027,  0.0004, -0.0012, -0.0010, -0.0008, -0.0002,  0.0003
)

mse <- tribble(
  ~param, ~`100`, ~`200`, ~`300`, ~`500`, ~`800`, ~`1000`, ~`2000`,
  "alpha",    0.0206, 0.0055, 0.0131, 0.0007, 0.0005, 0.0002, 0.0001,
  "lambda",   0.0374, 0.0139, 0.0098, 0.0051, 0.0037, 0.0028, 0.0013,
  "theta",    0.2836, 0.1222, 0.0955, 0.0577, 0.0444, 0.0330, 0.0156,
  "beta",     0.3128, 0.1127, 0.0679, 0.0431, 0.0266, 0.0219, 0.0105,
  "gamma0",   0.1117, 0.0485, 0.0347, 0.0212, 0.0128, 0.0102, 0.0055,
  "gamma1",   0.2811, 0.1279, 0.0954, 0.0505, 0.0299, 0.0253, 0.0122,
  "p0(x=1)",  0.0028, 0.0014, 0.0011, 0.0006, 0.0003, 0.0003, 0.0001,
  "p0(x=0)",  0.0040, 0.0018, 0.0013, 0.0008, 0.0005, 0.0004, 0.0002,
  "p1(x=1)",  0.0059, 0.0022, 0.0017, 0.0009, 0.0005, 0.0004, 0.0002,
  "p1(x=0)",  0.0018, 0.0008, 0.0006, 0.0003, 0.0002, 0.0002, 0.0001
)

cp <- tribble(
  ~param, ~`100`, ~`200`, ~`300`, ~`500`, ~`800`, ~`1000`, ~`2000`,
  "alpha",    0.9260, 0.9140, 0.9140, 0.9350, 0.9490, 0.9440, 0.9520,
  "lambda",   0.9340, 0.9190, 0.9280, 0.9290, 0.9330, 0.9440, 0.9510,
  "theta",    0.6860, 0.8010, 0.8320, 0.9130, 0.9420, 0.9330, 0.9470,
  "beta",     0.9470, 0.9390, 0.9450, 0.9430, 0.9490, 0.9460, 0.9490,
  "gamma0",   0.9610, 0.9620, 0.9420, 0.9400, 0.9500, 0.9500, 0.9420,
  "gamma1",   0.9550, 0.9570, 0.9450, 0.9520, 0.9540, 0.9480, 0.9600,
  "p0(x=1)",  0.9260, 0.9380, 0.9250, 0.9410, 0.9460, 0.9470, 0.9430,
  "p0(x=0)",  0.9470, 0.9550, 0.9390, 0.9390, 0.9490, 0.9480, 0.9450,
  "p1(x=1)",  0.9270, 0.9550, 0.9520, 0.9420, 0.9670, 0.9470, 0.9450,
  "p1(x=0)",  0.8740, 0.9400, 0.9300, 0.9330, 0.9410, 0.9360, 0.9570
)

# -----------------------------
# 2) Convert to long format
#    - RMSE = sqrt(MSE)
# -----------------------------
to_long <- function(df, metric_name) {
  df %>%
    pivot_longer(-param, names_to = "n", values_to = "value") %>%
    mutate(
      n = as.numeric(n),
      metric = metric_name
    )
}

df_long <- bind_rows(
  to_long(bias, "Bias"),
  to_long(mse,  "RMSE") %>% mutate(value = sqrt(value)),
  to_long(cp,   "Coverage probability")
) %>%
  mutate(
    metric = factor(metric, levels = c("Bias", "RMSE", "Coverage probability")),
    # keep legend/order exactly as in the table
    param  = factor(param, levels = bias$param),
    # x index for plotting but label with actual sample sizes
    x      = match(n, sizes)
  )

# Greek-style labels (preserving order)
param_labs <- c(
  alpha     = expression(alpha),
  lambda    = expression(lambda),
  theta     = expression(theta),
  beta      = expression(beta),
  gamma0    = expression(gamma[0]),
  gamma1    = expression(gamma[1]),
  "p0(x=1)" = expression(p[0](x==1)),
  "p0(x=0)" = expression(p[0](x==0)),
  "p1(x=1)" = expression(p[1](x==1)),
  "p1(x=0)" = expression(p[1](x==0))
)

param_order <- levels(df_long$param)

# -----------------------------
# 3) Dashed CP bands ONLY in CP panel
# -----------------------------
cp_bands <- tibble(
  metric = factor("Coverage probability", levels = levels(df_long$metric)),
  y = c(0.936, 0.964)   # (you wrote 0.961 in the comment, but used 0.964 in code)
)

# Line types / point shapes
ltys <- c(1,2,3,4,5,6,7,8,9,10)[seq_along(param_order)]
pchs <- c(15,16,17,18,19,0,1,2,3,4)[seq_along(param_order)]

# -----------------------------
# 4) Plot
# -----------------------------

p <- ggplot(
  df_long,
  aes(x = x, y = value, color = param, linetype = param, shape = param, group = param)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Coverage probability" ~
        scale_y_continuous(
          trans  = "logit",
          limits = c(0.68, 0.975),
          breaks = c(0.7, 0.8,  0.90, 0.936,0.964,0.975)
        ),
      metric == "RMSE" ~
        scale_y_continuous(
          breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6),
          limits = c(0, 0.6)
        )
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(sizes),
    labels = sizes,
    name = "Sample size (n)"
  ) +
  scale_color_discrete(breaks = param_order, labels = param_labs) +
  scale_linetype_manual(breaks = param_order, values = ltys, labels = param_labs) +
  scale_shape_manual(breaks = param_order, values = pchs, labels = param_labs) +
  geom_hline(
    data = cp_bands,
    aes(yintercept = y),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(
    y = "Value",
    color = "Parameter",
    linetype = "Parameter",
    shape = "Parameter"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    axis.text.x = element_text(angle = 25, vjust = 0.8, hjust = 1)
  )

p


ggsave(file = "simulation_Study1_asymptotic_properties.pdf", plot=p,width = 16, height = 7,dpi=300,bg="white")








## Plot RMSE associated to the cure rates


library(readxl)
library(tidyr)
library(dplyr)
library(tidyr)
library(ggplot2)


df <- read_excel("results_RMSE_as .xlsx")


# assumes df exists with columns: theta, Parameter, `100`,`200`,`300`,`500`,`800`,`1000`,`2000`
sizes <- c("100","200","300","500","800","1000","2000")

df_long <- df %>%
  filter(Parameter %in% c("p01","p00","p11","p10")) %>%
  pivot_longer(cols = all_of(sizes), names_to = "n", values_to = "RMSE") %>%
  mutate(
    n = as.numeric(n),
    theta = factor(theta),
    Parameter = factor(Parameter, levels = c("p01","p00","p11","p10")),
    x = match(as.character(n), sizes)  # equally spaced positions
  ) %>%
  mutate(
    Parameter = recode(
      Parameter,
      p01 = "p0(x=1)",
      p00 = "p0(x=0)",
      p11 = "p1(x=1)",
      p10 = "p1(x=0)"
    )
  )

thetas <- levels(df_long$theta)
pchs <- c(15, 16, 17, 18, 19, 0, 1, 2, 3, 4)[seq_along(thetas)]
ltys <- c(1, 2, 3, 4, 5, 6, 7, 8)[seq_along(thetas)]

plot1<-ggplot(df_long, aes(x = x, y = RMSE, group = theta,
                           color = theta, linetype = theta, shape = theta)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  facet_wrap(~Parameter, ncol = 2, scales = "free_y") +
  scale_x_continuous(
    breaks = seq_along(sizes),
    labels = sizes,
    name = "Sample size (n)"
  ) +
  scale_shape_manual(values = pchs) +
  scale_linetype_manual(values = ltys) +
  labs(y = "RMSE", color = expression(theta), linetype = expression(theta), shape = expression(theta)) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 16),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14)
  )


plot1



#Plot 2


df_long <- df %>%
  filter(Parameter %in% c("p10","p11")) %>%
  pivot_longer(cols = all_of(sizes), names_to = "n", values_to = "RMSE") %>%
  mutate(
    n = as.numeric(n),
    theta = factor(theta),
    x = match(as.character(n), sizes),
    Parameter = recode(
      Parameter,
      p11 = "p[1](x==1)",
      p10 = "p[1](x==0)"
    )
  )

thetas <- levels(df_long$theta)
pchs <- c(15, 16, 17, 18, 19, 0, 1, 2, 3, 4)[seq_along(thetas)]
ltys <- c(1, 2, 3, 4, 5, 6, 7, 8)[seq_along(thetas)]

p <- ggplot(df_long, aes(x = x, y = RMSE, group = theta,
                             color = theta, linetype = theta, shape = theta)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  facet_wrap(~Parameter, ncol = 2, labeller = label_parsed) +
  scale_x_continuous(
    breaks = seq_along(sizes),
    labels = sizes,
    name = "Sample size (n)"
  ) +
  scale_shape_manual(values = pchs) +
  scale_linetype_manual(values = ltys) +
  labs(y = "RMSE",
       color = expression(theta),
       linetype = expression(theta),
       shape = expression(theta)) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    axis.text.x = element_text(angle = 25, vjust = 0.8, hjust = 1)
  )

p

ggsave(file = "simulation_Study2_RMSE_cure_rate.pdf", plot=p,width = 13, height = 6,dpi=300,bg="white")



#### Simulation study - Hypothesis testing (rejection rate plot)


df <- tribble(
  ~effect, ~`100`, ~`200`, ~`300`, ~`500`, ~`800`, ~`1000`, ~`2000`, ~`5000`,
  -0.05,    13.9,   16.7,   25.0,   40.5,   60.8,    73.5,    96.3,   100.0,
  -0.04,    11.5,   12.5,   18.7,   27.0,   43.2,    53.2,    83.2,    99.8,
  -0.03,    10.1,    9.8,   11.8,   16.8,   26.1,    31.5,    61.7,    95.2,
  -0.02,     8.9,    8.3,    8.6,    9.5,   14.5,    15.2,    27.9,    66.0,
  -0.01,     7.7,    6.3,    6.3,    7.1,    7.2,     7.9,    11.1,    18.6,
  0.00,     7.2,    6.0,    4.8,    6.1,    5.3,     5.3,     5.3,     5.8,
  0.01,     8.0,    7.0,    6.5,    8.5,    7.6,     7.1,    10.7,    19.2,
  0.02,     9.9,    8.4,   10.7,   12.1,   16.3,    15.0,    29.4,    66.5,
  0.03,    10.3,   10.6,   14.5,   19.1,   30.1,    32.0,    60.0,    96.4,
  0.04,    12.9,   15.1,   20.2,   29.9,   46.8,    54.6,    86.8,   100.0,
  0.05,    16.3,   19.9,   28.6,   43.9,   66.4,    75.4,    97.4,   100.0
)

# --- Long format ---
sizes <- c(100, 200, 300, 500, 800, 1000, 2000, 5000)


df_long <- df %>%
  pivot_longer(cols = -effect, names_to = "n", values_to = "rej_rate") %>%
  mutate(
    effect = as.numeric(effect),
    n = factor(as.numeric(n), levels = sizes),
    rej_rate = rej_rate / 100
  )

# Aesthetics for each sample size line
n_levels <- levels(df_long$n)
pchs <- c(15, 16, 17, 18, 19, 0, 1, 2)[seq_along(n_levels)]
ltys <- c(1, 2, 3, 4, 5, 6, 7, 8)[seq_along(n_levels)]

p <- ggplot(df_long, aes(
  x = effect, y = rej_rate,
  group = n, color = n, linetype = n, shape = n
)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_continuous(
    breaks = sort(unique(df_long$effect)),
    name = expression(paste("Effect magnitude (", beta, " = ", alpha[1], " = ", gamma[1], ")"))
  ) +
  scale_y_sqrt(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0.0365, 1),
    breaks = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1),
    name = "Empirical rejection rate"
  ) +
  scale_shape_manual(values = pchs, name = "Sample size (n)") +
  scale_linetype_manual(values = ltys, name = "Sample size (n)") +
  scale_color_discrete(name = "Sample size (n)") +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18),
    axis.text.x = element_text(angle = 25, vjust = 0.8, hjust = 1)
  )

p

ggsave(file = "simulation_Study_rejection_rate.pdf", plot=p,width = 12, height = 8,dpi=300,bg="white")









## True survival function used in the simulation study (hypothesis testing - two-groups)



t  <- seq(0, 60, by = 0.01)
x0 <- rep(0, length(t))
x1 <- rep(1, length(t))

effect_magnitude <- seq(-0.05, 0.05, by = 0.01)  # 11 values

# build data for ALL effect sizes
df_plot <- do.call(rbind, lapply(effect_magnitude, function(em) {
  par <- c(-1, em, -0.1, em, 0.5, 0.4, em)
  
  s0 <- sobrevGTDL_alpha(t, x0, par)
  s1 <- sobrevGTDL_alpha(t, x1, par)
  
  data.frame(
    t = rep(t, 2),
    S = c(s0, s1),
    group = factor(rep(c("x=0", "x=1"), each = length(t)),
                   levels = c("x=0", "x=1")),
    effect = em
  )
}))


p <- ggplot(df_plot, aes(x = t, y = S, color = group, linetype = group)) +
  geom_line(linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(
    ~ effect,
    ncol = 3,
    labeller = labeller(effect = function(x) paste0("Effect: ", x))
  ) +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  scale_y_continuous(limits = c(0, 1),  breaks = seq(0, 1, 0.1)) +
  scale_color_discrete(
    labels = c("Control", "Intervention")
  ) +
  scale_linetype_discrete(
    labels = c("Control", "Intervention")
  ) +
  labs(
    x = "Time (months)",
    y = "Survival function",
    color = "Group",
    linetype = "Group"
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18)
  )

p


ggsave(file = "simulation_Study_rejection_rate_survival_plots.pdf", plot=p,width = 12, height = 16,dpi=300,bg="white")


#Median survival time

median_surv <- df_plot %>%
  arrange(effect, group, t) %>%
  group_by(effect, group) %>%
  summarise(
    median_t = {
      idx <- which(S <= 0.5)
      if (length(idx) == 0) NA_real_ else t[min(idx)]
    },
    S_at_end = S[which.max(t)],   # útil p/ checar se não cruzou 0.5
    .groups = "drop"
  )

median_surv


median_surv2 <- df_plot %>%
  arrange(effect, group, t) %>%
  group_by(effect, group) %>%
  summarise(
    median_t = {
      idx <- which(S <= 0.5)
      if (length(idx) == 0) Inf else t[min(idx)]
    },
    .groups = "drop"
  ) %>%
  mutate(median_t = ifelse(is.infinite(median_t), ">60", as.character(median_t)))

median_surv2


median_wide <- median_surv %>%
  mutate(median_t = round(median_t, 3)) %>%
  select(effect, group, median_t) %>%
  tidyr::pivot_wider(names_from = group, values_from = median_t)

median_wide



