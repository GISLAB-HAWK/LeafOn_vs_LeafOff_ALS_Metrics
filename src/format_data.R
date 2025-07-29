#' Format Data Frame by Lowercasing Column Names and Extracting Substrings
#'
#' This function takes a data frame as input and performs data formatting steps on it.
#' It converts column names to lowercase and extracts substrings from the column names
#' based on the position of the last underscore ("_") character. The extracted substrings
#' are then set as the new column names.
#'
#' @param df A data frame to be formatted.
#'
#' @return The formatted data frame with modified column names.
#'
#' @examples
#' df_formatted <- format_data(df)
#'

format_data <- function(df) {
  
  # column names to lowercase
  names(df) <- tolower(names(df))
  
  dat <- data.frame(name = names(df), pos_sep = NA)
  
  # get position of last separator in column name strings
  for (i in 1:dim(dat)[1]) {
    
    y <- dat[i,]
    sep <- unlist(gregexpr('_', y, fixed = T))
    dat[i,]$pos_sep <- max(sep)
    
  }
  
  # replace column names
  name_n <- substr(dat$name, dat$pos_sep + 1, nchar(dat$name))
  names(df) <- name_n
  
  # remove intermediate variables
  rm(dat, y, sep, i, name_n)
  
  return(df)
  
}