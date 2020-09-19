library(ggplot2)
library(dplyr)

# gpr_data <- read.csv("gpr_dlmt_best_points.csv") %>%
#     select(slowdown, method, measurement_order, time_per_pixel)

gpr_data <- read.csv("best_points.csv") %>%
    select(slowdown, method, measurement_order, time_per_pixel)

names(gpr_data) <- c("slowdown", "method", "point_number", "time_per_pixel")

gpr_lin_data <- read.csv("gpr_dlmt_best_points.csv") %>%
    select(slowdown, method, measurement_order, time_per_pixel) %>%
    mutate(method = "GPRM")

names(gpr_lin_data) <- c("slowdown", "method", "point_number", "time_per_pixel")

df_all_methods <- read.csv("dopt_anova_experiments/data/complete_1000.csv",
                           strip.white = TRUE, header = TRUE)

df_all_methods <- df_all_methods %>%
    select(slowdown, method, point_number, time_per_pixel) %>%
    bind_rows(gpr_data, gpr_lin_data) %>%
    mutate(method = factor(method,
                           levels = c("RS", "LHS", "GS", "GSR",
                                      "GA","LM", "LMB", "LMBT",
                                      "RQ", "DOPT", "DLM", "DLMT",
                                      "GPR", "GPRM"))) %>%
    filter(method %in% c("RS", "LHS", "GS", "GSR",
                         "GA", "LM", "DLMT", "GPR", "GPRM")) %>%
    group_by(method) %>%
    mutate(mean = mean(slowdown),
           median = median(slowdown),
           ci95 = 1.96 * sd(slowdown) / sqrt(n()),
           max = max(slowdown)) %>%
    ungroup()

ggplot(df_all_methods) +
    facet_grid(method ~ .) +
    theme_bw(base_size = 18) +
    coord_cartesian(xlim = c(.9, 4),
                    ylim = c(0, 1000)) +
    geom_histogram(aes(slowdown),
                   binwidth = .05,
                   fill = "gray48") +
    scale_y_continuous(breaks = c(0, 1000),
                       labels = c("0", "1000")) +
    geom_curve(aes(x = max + .1,
                   y = 500,
                   xend = max,
                   yend = 5),
               arrow = arrow(length = unit(0.05, "npc")),
               curvature = 0.3,
               stat = "unique") +
    geom_text(aes(x = max + .2,
                  y = 550,
                  label = "max"),
              stat = "unique") +
    geom_rect(aes(xmin = mean - ci95,
                  xmax = mean + ci95,
                  ymin = 0,
                  ymax = 1000,
                  fill = "red"),
              alpha = 0.3,
              stat = "unique") +
    geom_vline(aes(xintercept = median),
               color = "darkgreen",
               linetype = 3,
               stat = "unique") +
    geom_vline(aes(xintercept = mean),
               color = "red",
               linetype = 2,
               stat = "unique") +
    labs(y = "Count",
         x = "Slowdown compared to the optimal solution") +
    scale_fill_discrete(name = "",
                        breaks = c("red"),
                        labels = c("Mean error")) +
    ggtitle("") +
    theme(legend.position = "none",
          text = element_text(family="serif"),
          strip.background = element_rect(fill = "white"))
