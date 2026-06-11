#' Orange Tip Daily
#'
#' Daily site, year, and count data for orange tip 
#' 
#' @name orange_tip_daily
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(orange_tip_daily)
#'
#' @rdname orange_tip_daily
#'
#' @format A dataframe with 5 columns and 1048575 rows. The columns are DAYNO (day number), WEEKNO (week number), SITENO (site number), YEAR (year), SPECIES (species identify), COUNT (count), and DATE (date).
#' 
#' @examples
#'
#' library(ButterflyLifespan)
#'
#' data(orange_tip_daily)
#' head(orange_tip_daily)
#' plot(orange_tip_daily$COUNT)
#'
#'
NULL