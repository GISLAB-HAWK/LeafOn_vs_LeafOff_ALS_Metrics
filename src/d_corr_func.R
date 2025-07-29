#' Assign Specific IDs Based on Tree Species Group
#'
#' Assigns an ID based on the tree species group.
#'
#' @param bagr Data frame indexing the column with tree species groups.
#'
#' @return Numeric vector containing the IDs.
#'
#' @details
#' The function assigns an ID based on the tree species group.
#'
#' @examples
#' df$tree_species_group <- input_d_korr(df$tree_species_groups)
#'

input_d_korr = function(bagr){

  bagr_d_korr <- ifelse(bagr == "EI", 110,
	               ifelse(bagr == "BU", 211,
	               ifelse(bagr == "ALH", 300,
	               ifelse(bagr == "ALN", 400,
	               ifelse(bagr == "FI", 500,
	               ifelse(bagr == "DGL", 611,
	               ifelse(bagr == "KI", 700, 
	                      800)
	               )))))) 	
}

#' Calculate Corrected DBH (Diameter at Breast Height)
#'
#' Calculates the corrected DBH based on diameter, measurement height, and tree species group.
#'
#' @param du Numeric vector indicating the diameter in centimeters.
#' @param abwmh Numeric vector indicating the measurement height in centimeters.
#' @param ba Numeric vector indicating the tree species group.
#'
#' @return Numeric vector containing the corrected DBH values.
#'
#' @details
#' The function creates a data frame with input data (du, abwmh, ba) and calculates the corrected DBH using specific formulas for each tree species group.
#'
#' @examples
#' d <- d_korr(du = df$dbh, abwmh = df$dbh_diff, ba = df$tree_species)

d_korr = function(du, abwmh, ba){
  
  data  <- data.frame(du, abwmh + 130, ba)
  data$bhd_korr <- NA
  data$ba <- as.numeric(as.character(data$ba))

  try(data[data$ba==110 & is.na(data$ba)==F,]$bhd_korr<-(1-0.250973942811556*((data[data$ba==110 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==110 & is.na(data$ba)==F,]$du^(1+0.0960473724548668*((data[data$ba==110 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Eiche
  try(data[data$ba==113 & is.na(data$ba)==F,]$bhd_korr<-(1-0.250973942811556*((data[data$ba==113 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==113 & is.na(data$ba)==F,]$du^(1+0.0960473724548668*((data[data$ba==113 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #REiche
  try(data[data$ba==211 & is.na(data$ba)==F,]$bhd_korr<-(1+0.000000000000000*((data[data$ba==211 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==211 & is.na(data$ba)==F,]$du^(1+0.0266729142937414*((data[data$ba==211 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Buche
  try(data[data$ba==300 & is.na(data$ba)==F,]$bhd_korr<-(1+0.099714148602944*((data[data$ba==300 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==300 & is.na(data$ba)==F,]$du^(1+0.0037810019407898*((data[data$ba==300 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #ALH
  try(data[data$ba==400 & is.na(data$ba)==F,]$bhd_korr<-(1-0.022313178325787*((data[data$ba==400 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==400 & is.na(data$ba)==F,]$du^(1+0.0419980094357119*((data[data$ba==400 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #ALN
  try(data[data$ba==500 & is.na(data$ba)==F,]$bhd_korr<-(1-0.130077988347608*((data[data$ba==500 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==500 & is.na(data$ba)==F,]$du^(1+0.0751116090656222*((data[data$ba==500 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Fichte
  try(data[data$ba==520 & is.na(data$ba)==F,]$bhd_korr<-(1-0.226435038499499*((data[data$ba==520 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==520 & is.na(data$ba)==F,]$du^(1+0.0918404787867948*((data[data$ba==520 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Tanne
  try(data[data$ba==611 & is.na(data$ba)==F,]$bhd_korr<-(1+0.195680347916216*((data[data$ba==611 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==611 & is.na(data$ba)==F,]$du^(1-0.0230574154637055*((data[data$ba==611 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Douglasie
  try(data[data$ba==700 & is.na(data$ba)==F,]$bhd_korr<-(1+0.140928222424959*((data[data$ba==700 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==700 & is.na(data$ba)==F,]$du^(1-0.0081348607807616*((data[data$ba==700 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T) #Kiefer
  try(data[data$ba==800 & is.na(data$ba)==F,]$bhd_korr<-(1-0.257968477710732*((data[data$ba==800 & is.na(data$ba)==F,]$abwmh/100)-1.3))*data[data$ba==800 & is.na(data$ba)==F,]$du^(1+0.11497179651835*((data[data$ba==800 & is.na(data$ba)==F,]$abwmh/100)-1.3)),silent=T)   #Laerche

  return(data$bhd_korr)
}