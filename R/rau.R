#' Rational arcsine transformed units
#'
#' Converts proportion into arcsine transformed units.
#' @param p Proportion of outcome.
#' @param N Total number of observations
#' @keywords proportion, arcsine, rational
#' @export
#' @examples
#' rau(.5, 12)

rau = function(p, N) {
  AU = asin(sqrt(p / (N + 1))) + asin(sqrt(( p + 1) / (N + 1)))
  RAU =  146 / pi *AU - 23

  return(RAU)
}
