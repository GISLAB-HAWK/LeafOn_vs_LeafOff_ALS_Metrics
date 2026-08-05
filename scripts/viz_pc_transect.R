library(lidR)
library(ggplot2)

# source setup script
source('src/setup.R', local = TRUE)

seasons <- c("pc_leafon_2023_angle20", "pc_leafoff_2024_angle20")
ppms <- c("ppm20", "ppm10", "ppm4")

season <- seasons[2]
ppm <- ppms[1]



# Transekt definieren (100m)
x1 <- 548000
x2 <- 548150
y1 <- 5730000
y2 <- 5730000  
breite <- 5  


# Alle Kombinationen durchlaufen
plot_list <- list()

for (season in seasons) {
  for (ppm in ppms) {
    dir <- file.path(processed_data_dir, season, ppm)
    
    las <- readLAScatalog(dir)
    transekt <- clip_transect(las, p1 = c(x1, y1), p2 = c(x2, y2), width = breite)
    df <- transekt@data
    
    if (season == "pc_leafon_2023_angle20"){
      s_annot <- "leaf-on"
    } else
      s_annot <- "leaf-off"
    
    p <- ggplot(df, aes(X, Z, color = factor(ReturnNumber))) +
      geom_point(size = 0.25, alpha = 0.7) +
  
      scale_color_manual(
        values = c(
          "#000000",
          "#0072B2",
          "#009E73",
          "#E69F00",
          "#56B4E9",
          "#CC79A7",
          "#D55E00"
        ),
        name = "Return"
      ) +
      
      coord_equal() +
      
      labs(
        title = paste(s_annot, ppm, sep = " | "),
        x = "Easting (m)",
        y = "Height ASL (m)"
      ) +
      
      theme_minimal() +
      
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 10
        ),
        
        panel.grid.minor = element_blank(),
        legend.position = "none"
      )
    
    plot_list[[paste(season, ppm, sep = "_")]] <- p
  }
}


plot_list[[1]] <- plot_list[[1]] +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

plot_list[[2]] <- plot_list[[2]] +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

plot_list[[3]] <- plot_list[[3]] +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# untere Reihe -> nur links Y-Achse entfernen
plot_list[[5]] <- plot_list[[5]] +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)
    )
  )

plot_list[[6]] <- plot_list[[6]] +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )



combined <- wrap_plots(plot_list, ncol = 3) +
  plot_layout(guides = "collect") 

combined


ggsave(
  filename = paste0(output_dir,"pc_plot.png"),
  plot = combined,
  width = 14,
  height = 8,
  dpi = 300
)
