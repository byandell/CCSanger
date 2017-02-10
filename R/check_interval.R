#' Check chr, start_bp, end_bp for validity
#'
#' Simple checks.
#'
#' @param chr chromosome number
#' @param start_bp starting position in bp
#' @param end_bp ending position in bp
#'
#' @return none
check_interval <- function(chr, start_bp, end_bp) {
  if(missing(chr) | missing(start_bp) | missing(end_bp)) {
    stop("chr, start_bp and end_bp positions must be specified.")
  }
  if(start_bp < 0 | end_bp < 0 | end_bp < start_bp) {
    stop("end must be greater than start.")
  }
}
#' Convert to bp if in Mb
#'
#' Assumes positive value.
#'
#' @param value value in bp or Mbp
#' @param to_bp convert to bp if TRUE, Mbp if FALSE (default TRUE)
#'
#' @return value in bp
#' @export
#' 
convert_bp <- function(value, to_bp=TRUE) {
  if(max(value) < 200 & to_bp)
    value <- value * 1e6
  if(min(value) > 200 & !to_bp)
    value <- value * 1e-6
  value
}
