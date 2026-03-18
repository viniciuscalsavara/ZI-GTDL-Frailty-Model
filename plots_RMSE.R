library(readxl)
library(tidyr)
library(dplyr)
library(dplyr)
library(tidyr)
library(ggplot2)


df <- read_excel("results_RMSE_as .xlsx")

#View(df)

library(dplyr)
library(tidyr)
library(ggplot2)

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



#tLot 2


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

plot2 <- ggplot(df_long, aes(x = x, y = RMSE, group = theta,
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
    strip.text   = element_text(size = 18),
    axis.title   = element_text(size = 18),
    axis.text    = element_text(size = 18)
  )

plot2
ggsave(file = "simulation_Study2_RMSE_cure_rate.pdf", plot=plot2,width = 12, height = 6,dpi=300,bg="white")

