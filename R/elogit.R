#' Empirical logits
#'
#' Converts proportion into empirical logits (elogit) and corresponding weights (elogit_weights) for regression analysis.
#' @param p Proportion of outcome.
#' @param N Total number of observations
#' @keywords proportion, empirical logit
#' @export
#' @examples
#' elogit(.5, 12)

elogit = function(p, N) {
  y = p * N
  return(log( (y+.5)/ (N-y + .5) ))
}

elogit_weight = function(p, N) {
  y = p * N
  return(1 / ( 1 / ( y + .5) + 1 / (N - y + .5) ))
}
