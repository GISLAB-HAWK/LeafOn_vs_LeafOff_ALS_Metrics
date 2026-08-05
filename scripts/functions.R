library(managelidar)
library(lasR)
library(ggplot2)
library(units)
library(dplyr)
library(fs)

print_fileinfo <- function(vpc){
  filenames <- vpc |> get_names(full.names = TRUE)
  header <- get_header(filenames[1])
  h <- header[[1]]$header
  lasrinfo <- info(filenames[1])
  cat(glue::glue("
  
  additional Information:
  ----------------
  LAS Version          : {h$`Version Major`}.{h$`Version Minor`}
  Point Data Format ID : {h$`Point Data Format ID`}
  
  
  "))
}



print_pulseinfo <- function(vpc) {
  density <- vpc |> get_density()
  avg_pulse <- round(mean(density$pulsedensity, na.rm = TRUE), 1)
  avg_point <- round(mean(density$pointdensity, na.rm = TRUE), 1)
  
  penetration <- vpc |> get_penetration()
  to_pct <- function(x) round(mean(x, na.rm = TRUE) * 100, 1)
  
  avg_single   <- to_pct(penetration$single)
  avg_multiple <- to_pct(penetration$multiple)
  avg_two      <- to_pct(penetration$two)
  avg_three    <- to_pct(penetration$three)
  avg_four     <- to_pct(penetration$four)
  avg_five     <- to_pct(penetration$five)
  avg_six      <- to_pct(penetration$six)
  
  cat(glue::glue("
  Density (⌀):
  ----------------
  Pulse Density : {avg_pulse} pulses/m²
  Point Density : {avg_point} points/m²
  
  Pulse Penetration Rate (⌀):
  ----------------
  Single Returns   : {avg_single} %
  Multiple Returns : {avg_multiple} %
    Two Returns    : {avg_two} %
    Three Returns  : {avg_three} %
    Four Returns   : {avg_four} %
    Five Returns   : {avg_five} %
    Six Returns    : {avg_six} %
  "))
}



plot_density_comparison <- function(vpc_leafof, vpc_leafon,
                         y_label = "Density") {
  density_leafof <- get_density(vpc_leafof)
  density_leafof$group <- "Leaf-off"
  density_leafon <- get_density(vpc_leafon)
  density_leafon$group <- "Leaf-on"
  
  dplyr::bind_rows(density_leafof, density_leafon) |>
    tidyr::pivot_longer(
      cols = c(pointdensity, pulsedensity),
      names_to = "variable",
      values_to = "value"
    ) |>
    dplyr::mutate(
      variable = dplyr::recode(variable,
                               pointdensity = "Point density",
                               pulsedensity = "Pulse density")
    ) |>
    ggplot(aes(x = variable, y = value, fill = group)) +
    geom_boxplot(color = "black", linewidth = 0.3) +
    scale_fill_grey(start = 0.5, end = 1) +
    labs(x = NULL, y = y_label) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}





plot_penetration_comparison <- function(vpc_leafof, vpc_leafon){
  penetration_leafof <- get_penetration(vpc_leafof)
  penetration_leafof$group <- "leafof"
  penetration_leafon <- get_penetration(vpc_leafon)
  penetration_leafon$group <- "leafon"
  
  penetration_leafof |> 
    rbind(penetration_leafon) |> 
    tidyr::pivot_longer(cols = !c(filename, group),
                        names_to = "variable",
                        values_to = "value") |> 
    mutate(variable = factor(variable, levels = c("single", "multiple", "two", "three", "four", "five", "six"))) |> 
    ggplot(aes(x = variable, y = value, fill = group)) +
    geom_boxplot() +
    scale_fill_grey() +
    theme_minimal()
}

print_penetration_depth <- function(las){
  data <- las@data
  avg_first <- round(mean(data$HAG[data$ReturnNumber==1]), 1)
  avg_last <- round(mean(data$HAG[data$ReturnNumber==data$NumberOfReturns]), 1)
  avg_two <- round(mean(data$HAG[data$ReturnNumber==2]), 1)
  avg_three <- round(mean(data$HAG[data$ReturnNumber==3]), 1)
  avg_four <- round(mean(data$HAG[data$ReturnNumber==4]), 1)
  avg_five <- round(mean(data$HAG[data$ReturnNumber==5]), 1)
  avg_six <- round(mean(data$HAG[data$ReturnNumber==6]), 1)
  
  cat(glue::glue("
  
  Pulse Penetration Depth (⌀ m HAG):
  ----------------
  First Returns     : {avg_first} 
  Last Returns      : {avg_last} 
    Second Returns  : {avg_two} 
    Third Returns   : {avg_three} 
    Fourth Returns  : {avg_four}
    Fifth Returns   : {avg_five} 
    Sixth Returns   : {avg_six} 
  "))
}
