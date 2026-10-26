# The reporting scales, reachable from this package.
#
# `nsca_thresholds()` reports its two thresholds on a scale chosen by
# `threshold.x` and `threshold.y`, and the vocabulary for that choice lived
# only in SCAtools. An NSCA user who wanted to know what "percentile,
# directional" meant, or who wanted to put a number of their own on the same
# scale so that it could be compared with a reported threshold, had to know
# that the other package existed and attach it.
#
# So these delegate, in the same way and for the same reason as
# nsca_reference(): one implementation, one answer, one place to fix it. What
# they add is reachability -- and the guarantee that the scale an NSCA
# threshold is on and the scale a caller converts to are the same scale,
# because they are computed by the same function rather than by two that agree
# today.

#' Reporting scales available for NSCA threshold tables
#'
#' Lists the scales and conventions understood by [nsca_analysis()] and
#' [nsca_thresholds()], with a note on how each behaves when the direction
#' refers to a low level.
#'
#' Delegated to [SCAtools::sca_scales()]. Both thresholds in an NSCA table are
#' converted by that package's rule, so the list of scales is that package's
#' list; restating it here would create a second copy that would need to be
#' kept consistent.
#'
#' @return A data frame describing each scale.
#' @seealso [nsca_thresholds()], [nsca_rescale()], [SCAtools::sca_scales()]
#' @export
#' @examples
#' nsca_scales()
nsca_scales <- function() {
  SCAtools::sca_scales()
}

#' Convert values to an NSCA reporting scale
#'
#' The conversion that [nsca_thresholds()] uses, exported so that a quantity
#' obtained elsewhere can be placed on the same scale before it is compared
#' with a reported threshold.
#'
#' A necessity threshold on `"percentage.range"` and a percentage computed
#' elsewhere are comparable only if the latter was converted by the same rule,
#' including the direction-aware mirroring and the percentile strictness that
#' follows the inequality. Both components of an NSCA model share one axis and
#' one direction letter per axis, so one conversion serves both sides.
#'
#' Delegated to [SCAtools::sca_rescale()], which is the function that produced
#' the thresholds being compared against.
#'
#' @param values Numeric vector in actual units.
#' @param reference Observed values of the same variable, used by the
#'   `percentile` and `sd` scales.
#' @param bounds Length-two axis bounds, including any theoretical scope.
#'   Defaults to the range of `reference`.
#' @param scale One of the names returned by [nsca_scales()].
#' @param letter `"H"` or `"L"`, the direction letter for this axis.
#' @param convention `"absolute"` or `"directional"`.
#' @param inequality `"strict"` or `"inclusive"`; affects the `percentile`
#'   scale only.
#' @return A numeric vector on the requested scale. Non-finite inputs return
#'   `NA`.
#' @seealso [nsca_thresholds()], [nsca_scales()], [SCAtools::sca_rescale()]
#' @export
#' @examples
#' x <- c(1, 3, 4, 7, 9)
#' nsca_rescale(4, x, scale = "percentile", letter = "H")
#' nsca_rescale(4, x, scale = "percentile", letter = "L")
nsca_rescale <- function(values, reference,
                         bounds = range(reference, na.rm = TRUE),
                         scale = "percentage.range", letter = c("H", "L"),
                         convention = c("absolute", "directional"),
                         inequality = c("strict", "inclusive")) {
  SCAtools::sca_rescale(
    values = values, reference = reference, bounds = bounds, scale = scale,
    letter = match.arg(letter), convention = match.arg(convention),
    inequality = match.arg(inequality)
  )
}
