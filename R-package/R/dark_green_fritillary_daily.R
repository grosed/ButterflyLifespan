#' Dark Green Fritillary
#'
#' Daily site, year, and count data for dark green fritillary 
#' 
#' @name dark_green_fritillary_daily 
#' 
#' @docType data
#'
#' @keywords datasets
#'
#' @usage data(dark_green_fritillary_daily )
#'
#' @rdname dark_green_fritillary_daily 
#'
#' @format A dataframe with 5 columns and 581129 rows. The columns are DAYNO (day number), WEEKNO (week number), SITENO (site number), YEAR (year), SPECIES (species identify), COUNT (count), and DATE (date).
#' 
#' @examples
#'
#' library(ButterflyLifespan)
#'
#' data(dark_green_fritillary_daily )
#' head(dark_green_fritillary_daily )
#' plot(dark_green_fritillary_daily $COUNT)
#'
#'
NULL