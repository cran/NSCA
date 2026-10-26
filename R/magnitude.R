# The conventional bands, reachable from this package, and the line this
# package will not cross.
#
# The two *component* effect sizes are ordinary NCA effect sizes, so Dul's
# bands apply to them and citing those bands is citing somebody else's
# published convention. Delegated for that reason, the way nsca_reference()
# and nsca_scales() are: one implementation, one answer.
#
# The joint indices are a different matter. `weakest_effect`,
# `balanced_joint_effect` and `joint_empty_space_coverage` are constructions of
# this package. Nobody has published a convention for them, so a band would be
# invented rather than cited -- and it would be invented by the package that
# also defines the index, which is the shape of a claim that cannot be checked.
# The print method has said so since 0.4.0 and it still does.

#' Conventional magnitude bands for a component effect size
#'
#' This function calls [SCAtools::sca_magnitude()]. It is provided here so
#' that NSCA tables can be interpreted without attaching the other package.
#'
#' @section Which numbers this applies to:
#' `d_nec` and `d_suf`, and `effect_size` in the component tables, under the
#' default `scale = "effect_size"`. These are ordinary NCA effect sizes, and
#' Dul's bands are a published convention that can be cited for them.
#' `scale = "beta"` labels the standardized coefficient that [nsca_compare()]
#' reports beside them, following Cohen's convention. The two conventions
#' share their edges but shift their labels by one band, so a value of 0.35
#' is labeled "medium" as a coefficient and "large" as an effect size.
#'
#' The bands do not apply to `weakest_effect`, `balanced_joint_effect` or
#' `joint_empty_space_coverage`. Those are this package's own constructions;
#' no convention has been published for them, and bands supplied by the same
#' package that defines the indices would lack an independent basis. The
#' `print()` method states the same, and `nsca_power()` allows a design to be
#' calibrated by simulation instead.
#'
#' @section Read the band beside the floor:
#' A band label indicates where a value lies on a conventional scale. It does
#' not indicate how much of that value the geometry would produce under
#' random pairing alone. `d_nec_null` and `d_suf_null` from [nsca_table()]
#' address that question, and we suggest reporting the two together: a
#' `d_nec` of 0.10 labeled "medium" at a floor of 0.055 is a different finding
#' from the same 0.10 at a floor of 0.007.
#'
#' @param value Optional numeric vector. If omitted, the band table is
#'   returned.
#' @param scale The convention to use: `"effect_size"` for Dul's bands on a
#'   component effect size, `"beta"` for Cohen's bands on the standardized
#'   coefficient that `nsca_compare()` reports. Both use .10, .30 and .50 as
#'   edges, with the labels shifted by one band, which is why `beta_cut` and
#'   `d_cut` default to different numbers.
#' @return With `value`, a character vector of labels. Without it, a data
#'   frame of bands.
#' @references
#' Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
#' Sciences* (2nd ed.). Lawrence Erlbaum.
#'
#' Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology
#' of "necessary but not sufficient" causality. *Organizational Research
#' Methods*, 19(1), 10-52. \doi{10.1177/1094428115584005}
#' @seealso [nsca_table()], [nsca_compare()], [SCAtools::sca_magnitude()]
#' @export
#' @examples
#' nsca_magnitude()
#' nsca_magnitude(c(0.04, 0.18, 0.35, 0.62))
#' nsca_magnitude(c(0.04, 0.18, 0.35, 0.62), scale = "beta")
nsca_magnitude <- function(value = NULL,
                           scale = c("effect_size", "beta")) {
  SCAtools::sca_magnitude(value, scale = match.arg(scale))
}
