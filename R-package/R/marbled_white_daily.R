#' Marbled White
#'
#' Daily site, year, and count data for marbled white 
#' 
#' @name marbled_white_daily 
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(marbled_white_daily)
#'
#' @rdname marbled_white_daily 
#'
#' @format A dataframe with 5 columns and 1048575 rows. The columns are DAYNO (day number), WEEKNO (week number), SITENO (site number), YEAR (year), SPECIES (species identify), COUNT (count), and DATE (date).
#' 
#' @examples
#'
#' library(ButterflyLifespan)
#'
#' data(marbled_white_daily)
#' head(marbled_white_daily)
#' plot(marbled_white_daily$COUNT)
#'
#'
NULL