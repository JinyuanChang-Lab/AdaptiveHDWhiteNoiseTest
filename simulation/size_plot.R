
library(ggplot2)
library(dplyr)
library(ggnewscale)
library(latex2exp)


df_res <- read.csv("./result/plot_for_size/figure_simulation_results_new.csv") 

plot_data <- df_res %>%
  mutate(empirical_size = empirical_size / 100) %>%
  filter(df_t %in% c(3, 5, 10, Inf))


plot_data$df_t <- factor(plot_data$df_t, 
                         levels = c(3, 5, 10, Inf))

create_and_save_plot <- function(target_p_type, file_name) {
  
  sub_data <- plot_data %>% filter(p_type == target_p_type)
  y_max <- max(0.95, max(sub_data$empirical_size, na.rm = TRUE) + 0.1)
  
  p_plot <- ggplot() +
    geom_hline(yintercept = 0.05, linetype = "dotted", color = "black", linewidth = 0.8) +
    
  geom_line(data = filter(sub_data, df_t != "Inf"), 
            aes(x = n, y = empirical_size, color = df_t, linetype = df_t, group = df_t), 
            linewidth = 0.9) +
    geom_point(data = filter(sub_data, df_t != "Inf"), 
               aes(x = n, y = empirical_size, color = df_t, shape = df_t, group = df_t), 
               size = 3.5, stroke = 1.2, fill = "white") +
    
    scale_color_manual(
      name = expression(italic(t)~distribution),
      breaks = c("3", "5", "10"),
      values = c("#D73027", "#F46D43", "#4DAF4A", "#984EA3"),
      labels = expression(italic(d) == 3, italic(d) == 5, italic(d) == 10)
    ) +
    scale_linetype_manual(
      name = expression(italic(t)~distribution),
      breaks = c("3", "5", "10"),
      values = c("dotdash", "longdash", "dashed", "twodash"),
      labels = expression(italic(d) == 3, italic(d) == 5, italic(d) == 10)
    ) +
    scale_shape_manual(
      name = expression(italic(t)~distribution),
      breaks = c("3", "5", "10"),
      values = c(24, 25, 23, 22), 
      labels = expression(italic(d) == 3, italic(d) == 5, italic(d) == 10)
    ) +
    
  new_scale_color() +
    new_scale("linetype") +
    new_scale("shape") +
    
  geom_line(data = filter(sub_data, df_t == "Inf"), 
            aes(x = n, y = empirical_size, color = df_t, linetype = df_t, group = df_t), 
            linewidth = 0.9) +
    geom_point(data = filter(sub_data, df_t == "Inf"), 
               aes(x = n, y = empirical_size, color = df_t, shape = df_t, group = df_t), 
               size = 3.5, stroke = 1.2, fill = "white") +
    

    scale_color_manual(
      name = NULL, 
      breaks = "Inf",
      values = "#377EB8",
      labels = expression(normal~distribution)
    ) +
    scale_linetype_manual(
      name = NULL,
      breaks = "Inf",
      values = "solid",
      labels = expression(normal~distribution)
    ) +
    scale_shape_manual(
      name = NULL,
      breaks = "Inf",
      values = 21, 
      labels = expression(normal~distribution)
    ) +
    

  scale_y_continuous(breaks = c(0.05, 0.30, 0.60, 0.90), 
                     limits = c(0, y_max)) +
    scale_x_continuous(breaks = unique(sub_data$n)) +
    
    theme_bw(base_size = 15, base_family = "serif") + 
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", linewidth = 1.2),
      
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      
      legend.position = c(0.05, 0.95), 
      legend.justification = c("left", "top"),
      
      legend.background = element_blank(), 
      
      legend.box.background = element_rect(fill = alpha("white", 0.85), 
                                           color = "black", linewidth = 0.5),
      
      legend.box.just = "left", 

      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.width = unit(1.5, "cm"),
      legend.spacing.y = unit(0.1, "cm")
    ) +
    labs(x = expression(n), y = NULL)
  
  ggsave(filename = file_name, plot = p_plot, 
         width = 5.5, height = 5, device = cairo_pdf)
  
  cat(sprintf("save to: %s\n", file_name))
  return(p_plot)
}

p1 <- create_and_save_plot("n",      "./result/plot_for_size/plot_size_n_50400.pdf")
p2 <- create_and_save_plot("n^2",    "./result/plot_for_size/plot_size_n2_50400.pdf")