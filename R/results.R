# Joint results: component effects, the three joint summaries, the admissible region,
# the geometry diagnostics and the joint-support decision.

# Conjunction rule for the two directional component tests. Requiring both to
# pass is an intersection-union combination, which is level-alpha without any
# assumption that the two are independent. The underlying null is random
# pairing, so this remains an empirical screen rather than a direct test of
# the full logical null "not necessary or not sufficient".
.nsca_verdict <- function(p_nec, p_suf, alpha, tested) {
  if (!isTRUE(tested)) {
    return("not tested")
  }
  if (!is.finite(p_nec) || !is.finite(p_suf)) {
    return("not tested")
  }
  if (p_nec <= alpha && p_suf <= alpha) {
    return("supported")
  }
  if (p_nec > alpha && p_suf > alpha) {
    return("neither component")
  }
  if (p_nec > alpha) "necessity fails" else "sufficiency fails"
}

# Names retired in 0.3.0, 0.4.0 and 0.5.x, mapped to their replacements. Kept
# reachable rather than merely documented, because a renamed column that
# silently returns NULL turns a working script into a wrong one instead of a
# broken one.
.nsca_legacy_columns <- c(
  # 0.3.0: name the geometry rather than the interpretation.
  d_nsca = "weakest_effect",
  j_nsca = "balanced_joint_effect",
  determinacy = "joint_empty_space_coverage",
  weaker_side = "weaker_component",
  p_nsca = "p_nsca_iut",
  # 0.4.0: adopt the vocabulary of condition analysis in degree.
  undetermined_share = "admissible_region_share",
  data_zone_share = "admissible_region_share",
  data_zone_width = "necessity_sufficiency_interval",
  undetermined_width = "necessity_sufficiency_interval",
  nsca_supported = "joint_support",
  conjunction = "significance_status",
  # `region` was renamed in 0.3.0 along with the rest of this block, and the
  # rename is listed in NEWS and in both READMEs, but the entry was never added
  # here -- so the one promise the retirement policy makes, that an old name
  # still leads somewhere, did not hold for it until 0.4.3.
  region = "threshold_status",
  boundary = "inequality",
  # 0.5.0: the status column looked only at the two p-values while the column
  # beside it, joint_support, also required relevance and acceptable geometry,
  # so one row could read "supported" and FALSE at once. The rename says which
  # question it answers. "empty zone" gives way to the framework's "empty
  # space" in the same pass.
  joint_support_status = "significance_status",
  # 0.5.1: `ceiling` named the estimator in the joint table while
  # `necessity_boundary` in the same row used "ceiling" for the boundary
  # itself. The column is `frontier` now, the name the component tables have
  # always used for it, and "ceiling" means only the boundary type.
  ceiling = "frontier",
  joint_empty_zone_coverage = "joint_empty_space_coverage"
)

# `empty_zone_area` was renamed to `empty_space_area` in 0.5.0 as well, but it
# is not listed here: it is a column of the component tables, which come from
# SCAtools and do not pass through .nsca_append_legacy(). Listing it would put
# a name in this map that no table served by the map can supply.
# SCAtools::sca_legacy_names() documents that one.

#' Names retired in NSCA 0.3.0, 0.4.0 and 0.5.x
#'
#' The 0.3.0 vocabulary followed the geometry rather than the interpretation:
#' the region between the two frontiers was named for what it is, and the two
#' joint indices for how they treat a weak component rather than for a
#' verdict.
#'
#' 0.4.0 adopted the vocabulary of *condition analysis in degree*. The region
#' compatible with both components is the **admissible region**; the distance
#' on the X scale between the two thresholds at one outcome target is the
#' **necessity-sufficiency interval**; evidence for both components is
#' **joint support**; and *boundary* is reserved for the theoretical line
#' that separates an expected empty space, so the argument controlling the
#' strictness of a reported inequality became `inequality`.
#'
#' 0.5.x replaced the shorthand "empty zone" with the framework's "empty
#' space", renamed `joint_support_status` to `significance_status` because it
#' considers only the two p-values, and named the technique column of the
#' joint tables `frontier`, as in the component tables, so that "ceiling"
#' refers only to the boundary type.
#'
#' Retired column names remain available. [nsca_table()] and
#' [nsca_thresholds()] append them as duplicate columns when called with
#' `legacy = TRUE`, each an exact copy of the current column. [nsca_extract()]
#' no longer accepts them: it stops with an error that names the replacement.
#' The retired argument `boundary` was removed in 0.5.0; its replacement is
#' `inequality`.
#'
#' @return A named character vector: names are the retired names, values the
#'   current ones.
#' @seealso [nsca_terms()]
#' @export
#' @examples
#' nsca_legacy_names()
nsca_legacy_names <- function() {
  .nsca_legacy_columns
}

#' Terminology of condition analysis in degree
#'
#' A machine-readable glossary linking the vocabulary of the condition analysis
#' in degree framework to the names this package uses. Names retired along the
#' way are listed by [nsca_legacy_names()].
#'
#' @return A data frame with one row per term.
#' @references
#' Chan, Y. (in press). *From Necessity to Sufficiency and
#' Necessity-and-Sufficiency: Extending Condition Analysis in Degree*.
#' University of Leeds.
#' @seealso [nsca_legacy_names()], [nsca_direction_map()], [nsca_table()]
#' @export
#' @examples
#' nsca_terms()
nsca_terms <- function() {
  data.frame(
    term = c(
      "scope", "expected empty space", "boundary", "frontier",
      "ceiling", "floor", "necessity threshold", "sufficiency threshold",
      "admissible region", "necessity-sufficiency interval",
      "threshold correspondence", "exact correspondence", "joint support",
      "component effect size", "frontier crossing", "central tendency",
      "component sum", "critical value"
    ),
    meaning = c(
      "Region of the X-Y space bounded by the theoretical or observed extremes of X and Y",
      "Part of the scope expected to contain no observations if the relation holds",
      "Theoretical line separating the expected empty space from the compatible region",
      "The boundary as estimated from data",
      "A boundary that limits how high Y can be for a given X",
      "A boundary that limits how low Y can be for a given X",
      "X level below which a stated outcome target cannot be attained",
      "X level at which the condition is sufficient for that target",
      "Part of the scope compatible with both components",
      "Distance on the X scale between the two thresholds at one outcome target",
      "How close the two thresholds are; assessed only after both components are supported",
      "The limiting case in which one X value is both required and sufficient",
      "Credible evidence for both components under their respective criteria",
      "One component's empty space as a proportion of the shared scope",
      "The two frontiers cross, so the component effect sizes can sum to more than one",
      "How the expected outcome changes with the condition; an average effect, not a boundary",
      "The sum of the two component effect sizes: the share of the scope that the two corners leave empty together",
      "The value that random pairing exceeds in only a proportion alpha of permutations; a statistic needs to exceed it for its permutation test to reject"
    ),
    package_name = c(
      "scope", "empty_space_area", "(theoretical; not estimated)",
      "frontier",
      "boundary_type == \"ceiling\"", "boundary_type == \"floor\"",
      "necessity_threshold", "sufficiency_threshold",
      "admissible_region_share", "necessity_sufficiency_interval",
      "threshold_status", "threshold_status == \"exact correspondence\"",
      "joint_support", "d_nec, d_suf", "frontiers_crossed, overlap_share",
      "reference = \"ols\", nsca_reference()",
      "component_sum, p_sum_perm",
      "d_nec_critical, d_suf_critical, component_sum_critical, weakest_effect_critical"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Central-tendency reference lines
#'
#' Fits the central-tendency lines for each condition, on the same observations
#' both frontiers were fitted to. `"ols"` is the ordinary least-squares
#' regression; `"qr"` is a quantile regression at `reference.tau`, which
#' defaults to `0.5` -- the median line, a central-tendency summary that is
#' less sensitive to individual observations than the least-squares line.
#'
#' Neither is a component of a condition claim. Both describe how the
#' expected outcome changes with the condition; the two frontiers describe
#' which condition-outcome combinations are absent. A reference line
#' therefore has a slope and an intercept but no empty space, no effect size,
#' no permutation p-value and no threshold.
#'
#' The computation is delegated to [SCAtools::sca_reference()] rather than
#' reimplemented. The sufficiency layer of an NSCA result is an `sca_result`
#' fitted to the same observations, so the line computed there is the line
#' this package would compute, and a single implementation is maintained.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x Optional condition names or positions. All conditions by default.
#' @param line Which reference lines to report. Defaults to those requested
#'   when the model was fitted, or to `"ols"` if none were requested.
#' @return A data frame with one row per condition and reference line.
#' @seealso [nsca_analysis()], [SCAtools::sca_reference()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
#'                      reference = c("ols", "qr"))
#' nsca_reference(fit)
nsca_reference <- function(model, x = NULL, line = NULL) {
  # Validates `x` against this package's condition list, so the error names
  # NSCA's own argument rather than surfacing from the other package.
  conditions <- .nsca_select_conditions(model, x)
  SCAtools::sca_reference(model$sufficiency, x = conditions, line = line)
}

.nsca_global_value <- function(fit, condition, parameter) {
  item <- fit$summaries[[condition]]
  if (is.null(item) || is.null(item$global) ||
      !(parameter %in% rownames(item$global))) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(
    item$global[match(parameter, rownames(item$global)), 1L]
  ))
}

.nsca_component_table <- function(model, component) {
  metadata <- .nsca_metadata(model)
  if (identical(component, "sufficiency")) {
    source <- SCAtools::sca_table(model$sufficiency)
    if (nrow(source) == 0L) {
      return(source)
    }
    source$component <- "sufficiency"
    source$nsca_number <- unname(.nsca_direction_number[source$direction])
    source$statement <- vapply(
      seq_len(nrow(source)),
      function(i) .nsca_sufficiency_statement(
        source$direction[[i]], source$condition[[i]], source$outcome[[i]]
      ),
      character(1L)
    )
    return(source[c(
      "component", "condition", "outcome", "direction", "relationship",
      "nsca_number", "empty_corner", "location", "empty_space",
      "boundary_type", "statement", "frontier",
      "observations", "empty_space_area", "effect_size",
      "effect_size_null", "effect_size_excess", "effect_size_critical",
      "p_value",
      "p_accuracy", "frontier_accuracy", "fit", "slope", "intercept"
    )])
  }

  rows <- list()
  index <- 0L
  for (condition in metadata$conditions) {
    direction <- unname(metadata$direction[[condition]])
    corner <- unname(metadata$necessity_corner[[condition]])
    for (ceiling in metadata$ceilings) {
      index <- index + 1L
      rows[[index]] <- data.frame(
        component = "necessity",
        condition = condition,
        outcome = metadata$outcome,
        direction = direction,
        relationship = .nsca_relationship_text(
          direction, condition, metadata$outcome
        ),
        nsca_number = unname(.nsca_direction_number[[direction]]),
        empty_corner = corner,
        location = SCAtools::sca_corner_location(corner),
        empty_space = .nsca_empty_space_text(
          corner, condition, metadata$outcome
        ),
        boundary_type = .nsca_boundary_type(corner),
        statement = .nsca_necessity_statement(
          direction, condition, metadata$outcome
        ),
        frontier = ceiling,
        observations = .nsca_global_value(
          model$necessity, condition, "Number of observations"
        ),
        empty_space_area = .nsca_nca_param(
          model$necessity, condition, ceiling, "Ceiling zone"
        ),
        effect_size = .nsca_nca_param(
          model$necessity, condition, ceiling, "Effect size"
        ),
        effect_size_null = .nsca_null_mean(
          model$necessity, condition, ceiling
        ),
        effect_size_excess = .nsca_excess(
          .nsca_nca_param(model$necessity, condition, ceiling, "Effect size"),
          .nsca_null_mean(model$necessity, condition, ceiling)
        ),
        effect_size_critical = .nsca_critical(
          .nsca_null_draws(model, NULL, "necessity", condition, ceiling),
          metadata$alpha
        ),
        p_value = .nsca_nca_param(
          model$necessity, condition, ceiling, "p-value"
        ),
        p_accuracy = .nsca_nca_param(
          model$necessity, condition, ceiling, "p-accuracy"
        ),
        frontier_accuracy = .nsca_as_proportion(.nsca_nca_param(
          model$necessity, condition, ceiling, "Ceiling accuracy"
        )),
        fit = .nsca_as_proportion(
          .nsca_nca_param(model$necessity, condition, ceiling, "Fit")
        ),
        slope = .nsca_nca_param(model$necessity, condition, ceiling, "Slope"),
        intercept = .nsca_nca_param(
          model$necessity, condition, ceiling, "Intercept"
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.nsca_band_for <- function(metadata, model, condition, ceiling) {
  raw <- metadata$raw$conditions[[condition]]
  x_bounds <- raw$x_bounds
  y_bounds <- raw$y_bounds
  grid <- seq(x_bounds[[1L]], x_bounds[[2L]], length.out = .nsca_grid_size)

  nec_corner <- unname(metadata$necessity_corner[[condition]])
  suf_corner <- unname(metadata$sufficiency_corner[[condition]])

  nec <- .nsca_frontier(
    raw$x, raw$y, nec_corner, ceiling, grid,
    slope = .nsca_nca_param(model$necessity, condition, ceiling, "Slope"),
    intercept = .nsca_nca_param(model$necessity, condition, ceiling, "Intercept"),
    peers = model$necessity$peers[[ceiling]][[condition]]
  )
  suf <- .nsca_frontier(
    raw$x, raw$y, suf_corner, ceiling, grid,
    slope = .nsca_sca_param(model$sufficiency, condition, ceiling, "Slope"),
    intercept = .nsca_sca_param(model$sufficiency, condition, ceiling, "Intercept"),
    peers = model$sufficiency$peers[[ceiling]][[condition]]
  )

  if (identical(.nsca_corner_side(nec_corner), "upper")) {
    band <- .nsca_band(nec, suf, grid, y_bounds)
  } else {
    band <- .nsca_band(suf, nec, grid, y_bounds)
  }
  scope_area <- diff(x_bounds) * diff(y_bounds)
  share <- function(value) {
    if (is.finite(value) && scope_area > 0) value / scope_area else NA_real_
  }
  list(
    share = share(band$area),
    overlap = share(band$overlap),
    # A step frontier always overlaps itself by roughly one tread; only a
    # smooth frontier that inverts is reporting a problem.
    #
    # The threshold is the caller's geometry.max_overlap, not a hard-coded
    # 0.01. Until 0.5.0 this line ignored the setting while geometry_acceptable
    # honored it, so a caller who widened the tolerance still got the warning
    # at 1% and then a geometry_acceptable of TRUE -- two statements about the
    # same row, disagreeing.
    suspect = !(ceiling %in% .nsca_step_ceilings) &&
      isTRUE(share(band$overlap) > metadata$geometry[["max_overlap"]])
  )
}

# The overlap a step frontier produces by construction, which is not a fault
# and must not count against the geometry. Between two consecutive
# observations the upper staircase still holds the earlier point's value while
# the lower staircase already holds the later one, so the two empty zones
# overlap on every tread. Summed over the n - 1 treads this is exactly
# 1 / (n - 1) of the scope, and it is attained by a perfect diagonal.
# Numerical slack for the "component effect sizes sum above one" check. Both
# shares are integrals over a grid, so an exactly saturating pair can land a
# few ulps above one; flagging that as crossed frontiers would be a false
# alarm about the one configuration the check is meant to bless.
.nsca_sum_tolerance <- 1e-9

# Columns that belong to nsca_thresholds() rather than nsca_table(). They are
# defined per outcome level, so nsca_extract() cannot return one of them for a
# condition and frontier and says so instead of calling the name unknown.
.nsca_threshold_columns <- c(
  "convention", "condition_scale", "outcome_scale", "outcome_level",
  "outcome_level_actual", "necessity_threshold", "sufficiency_threshold",
  "necessity_threshold_actual", "sufficiency_threshold_actual",
  "necessity_status", "sufficiency_status", "threshold_gap_actual",
  "necessity_sufficiency_interval", "overlap_width",
  "necessity_sufficiency_interval_actual", "overlap_width_actual",
  "thresholds_compatible",
  "threshold_status", "sufficiency_rule", "inequality"
)

.nsca_step_allowance <- function(ceiling, observations) {
  if (!(ceiling %in% .nsca_step_ceilings)) {
    return(0)
  }
  if (!is.finite(observations) || observations < 2) {
    return(0)
  }
  1 / (observations - 1)
}

#' Joint necessity and sufficiency results
#'
#' One row per condition and frontier technique. Both sides of a row use the
#' same technique, the same observations and the same scope.
#'
#' @section What is primary:
#' The component effect sizes `d_nec` and `d_suf` are primary. The remaining
#' indices combine them, and three such combinations are reported rather than
#' one, because the geometry of two empty zones does not single out any one
#' of them.
#'
#' \describe{
#'   \item{`weakest_effect`}{`min(d_nec, d_suf)`, the fully non-compensatory
#'   summary. It is zero whenever either component is zero and remains on the
#'   components' own scale, so its maximum is about 0.5 rather than 1.}
#'   \item{`balanced_joint_effect`}{`2 * sqrt(d_nec * d_suf)`, normalized onto
#'   `[0, 1]`. Partially compensatory: a larger component raises it at a
#'   diminishing rate, and it is zero whenever either component is zero. At
#'   an equal component sum of 0.80, `(0.40, 0.40)` gives 0.80 and
#'   `(0.70, 0.10)` gives 0.529. It is not a measure of balance:
#'   `(0.50, 0.125)` and `(0.25, 0.25)` both score 0.50. The `balance` column
#'   measures balance separately.}
#'   \item{`joint_empty_space_coverage`}{The share of the scope that at least
#'   one claim rules out, `1 - admissible_region_share`. Curvature-neutral: it
#'   reaches 1 whenever the admissible region closes, regardless of the shape
#'   of the relation.}
#' }
#'
#' `min` and the geometric mean are the power means at `p = -Inf` and `p = 0`,
#' so choosing between them amounts to choosing how much compensation to
#' allow; the component sum corresponds to `p = 1`, the fully compensatory
#' end, and is less suited to summarizing a conjunction because one strong
#' component can dominate it. [nsca_joint()] exposes the family directly.
#'
#' @section The area identity:
#' The two empty spaces and the admissible region tile the scope, so by
#' inclusion and exclusion
#'
#' \deqn{\mathrm{admissible\_region\_share} = 1 - d_{nec} - d_{suf} + O}
#'
#' where \eqn{O} is `overlap_share`. The remaining properties follow from
#' this identity rather than being asserted separately. Because
#' `admissible_region_share` is a share of the scope, it lies in `[0, 1]`
#' unconditionally, and therefore so does `joint_empty_space_coverage`. The
#' other two indices are bounded only up to the overlap:
#'
#' \deqn{d_{nec} + d_{suf} \le 1 + O, \quad
#'       \mathrm{weakest\_effect} \le 0.5 + O/2, \quad
#'       \mathrm{balanced\_joint\_effect} \le 1 + O}
#'
#' Both bounds are tight and are attained exactly by a perfect diagonal, where
#' `O` is `1 / (n - 1)`: at `n = 40` the observed values are 0.5128 and 1.0256
#' against bounds of 0.5128 and 1.0256. The excess reflects discretization
#' rather than evidence, and the indices are left unclamped so that it
#' remains visible.
#'
#' @section Vocabulary:
#' The column names follow the condition analysis in degree framework.
#' `relationship` states the direction as a claim about the theorized X-Y
#' relationship, as in Table 1 of that framework, not about the shape of the
#' estimated frontier. `necessity_empty_space` and `sufficiency_empty_space`
#' name what each expected empty space would contain, and
#' `necessity_boundary` and `sufficiency_boundary` state whether the boundary
#' separating it is a ceiling, which limits how high the outcome can be for a
#' given condition value, or a floor, which limits how low. One of the two is
#' a ceiling and the other a floor, because the two corners are diagonally
#' opposite.
#'
#' @section The disjointness check:
#' Whenever the two frontiers do not cross, the floor lies at or below the
#' ceiling throughout the scope, so no point belongs to both expected empty
#' spaces. Because each component effect size is that space's share of the
#' same scope, the two shares cannot sum to more than one. Component effect
#' sizes that sum to more than one therefore indicate that the frontiers have
#' crossed and that the joint result requires diagnosis before
#' interpretation, not that it is unusually strong. `component_sum` reports
#' `d_nec + d_suf` so that this rule can be applied directly.
#'
#' `frontiers_crossed` addresses the same question from the measurement
#' rather than from the sum. The sum serves as a proxy when the crossing
#' itself is not available; here it is available, because `overlap_share` is
#' the area on which the two empty spaces overlap, integrated from the
#' reconstructed frontiers. Using the measurement avoids a limitation of the
#' proxy: the sum and the allowance below are computed by different routes
#' and can disagree in the last few digits on precisely the configuration the
#' check is intended to accept.
#'
#' That allowance exists because a step frontier overlaps itself by
#' construction. Between two consecutive observations the upper staircase
#' still holds the earlier value while the lower one has already moved, so
#' the two spaces overlap on every tread; summed over the `n - 1` treads this
#' is exactly `1 / (n - 1)` of the scope, attained by a perfect diagonal,
#' which is an ideal relation rather than a failure. For every smooth
#' frontier the allowance is zero, and any overlap indicates a crossing.
#'
#' `frontiers_crossed` asks whether the frontiers crossed at all;
#' `geometry_acceptable` asks whether they crossed by more than
#' `geometry.max_overlap`. A row can be flagged as crossed and still have
#' acceptable geometry.
#'
#' @section Geometry diagnostics:
#' `reconstruction_error` is the residual of the identity above. The two
#' sides come from different sources: `admissible_region_share` and
#' `overlap_share` are integrated from the frontiers this package
#' reconstructs, while `d_nec` and `d_suf` are reported by the engines from
#' their own internal frontiers. A residual near zero indicates that the two
#' agree; a large one indicates that the reconstruction differs from what the
#' engine fitted. A plot may not reveal this, because an incorrectly
#' reconstructed frontier still looks like a frontier.
#'
#' `degenerate` marks a condition whose axis does not vary inside the declared
#' scope. `geometry_acceptable` combines both diagnostics with the excess
#' frontier overlap, allowing a step frontier the `1 / (n - 1)` it produces by
#' construction.
#'
#' @section Decisions:
#' Each p-value addresses a different question, and `joint_support` is based
#' on the first of them.
#'
#' * `p_nsca_iut` is `max(p_nec, p_suf)`: whether each of the two corners is
#'   emptier than random pairing produces. It passes only when both
#'   components reject at `test.p_threshold`. Against random pairing it is
#'   level-alpha without any assumption that the two component tests are
#'   independent (the intersection-union rule), so no correction is needed
#'   merely for combining two pre-specified tests. It is typically
#'   conservative, so an NSCA design tends to require more observations than
#'   either component alone. Its reliability depends on that of the component
#'   tests it combines; `?NSCA` describes a situation in which one of them may
#'   reject too readily.
#' * `p_weakest_perm` asks whether the smaller of the two empty areas,
#'   `weakest_effect`, is larger than random pairing produces. Its null is
#'   random pairing of both sides at once, so it is calibrated against that,
#'   and not against a state in which one component holds and the other does
#'   not. It concerns the smaller area, whereas `p_nsca_iut` concerns the
#'   weaker evidence, so the two can disagree when one side's random-pairing
#'   floor is high. On one shared sequence it is at most the p-value of the
#'   component whose observed area is the smaller, so it can reject where
#'   that component's own test does not; it is therefore not suited to serve
#'   as evidence that both corners are empty.
#' * `p_sum_perm` asks whether the two corners together, `component_sum =
#'   d_nec + d_suf`, are emptier than random pairing produces. A single
#'   strongly empty corner can be enough to make it small, so it describes the
#'   two corners jointly and complements `p_nsca_iut`, which asks the same of
#'   each corner separately. `joint_empty_space_coverage` differs from
#'   `component_sum` only by `overlap_share`, which is zero unless the two
#'   frontiers cross.
#'
#' `p_weakest_perm` and `p_sum_perm` require `shared.test.rep`, because the
#' null distribution of a statistic involving both components depends on how
#' the two co-vary under one shuffle. We suggest choosing the p-value that
#' corresponds to the research question before examining the results;
#' reporting several p-values as findings raises a multiple-testing question
#' (see `p.adjust`). Multiplicity across several conditions or frontiers is a
#' separate matter.
#'
#' `p_source` reports whether the p-values came from the engines' own
#' separate runs or from the shared sequence.
#'
#' @section Null distributions:
#' Each null summary in a row is taken from the distribution that produced
#' its p-value: the shared sequence when there is one, otherwise each
#' engine's own test. `d_nec_null` and `d_suf_null` are the means: what random
#' pairing alone produces at this sample size, scope and technique.
#' `d_nec_critical`, `d_suf_critical`, `component_sum_critical` and
#' `weakest_effect_critical` are the `1 - alpha` quantiles: the value that
#' random pairing exceeds in only a proportion `alpha` of permutations, which
#' each statistic needs to exceed for a single test at `alpha` to reject.
#' Each statistic has its own critical value; that of the sum is not the sum
#' of the two. The last two require a shared sequence. Adjustment for
#' multiple testing does not change them; see the section on multiple
#' testing.
#'
#' `p_nec_accuracy`, `p_suf_accuracy`, `p_weakest_perm_accuracy` and
#' `p_sum_perm_accuracy` give the Monte Carlo margin of each p-value, as NCA
#' reports its `p-accuracy`: the half-width of an interval at
#' `test.p_confidence`. They are zero when every ordering of the outcome was
#' evaluated, since the p-value is then exact.
#'
#' `necessity_supported` and `sufficiency_supported` require significance and,
#' when a `relevance` threshold was pre-specified, an effect at least that
#' large. `joint_support` additionally requires `geometry_acceptable`, so a
#' degenerate or poorly reconstructed row cannot be reported as supported.
#'
#' @section Blind spot:
#' All three joint summaries are areas over the declared scope and share one
#' blind spot. A constant outcome at the center of an imposed scope leaves the
#' upper and lower halves both empty, so both components reach 0.5, the
#' admissible region closes, and every index reaches its maximum for data
#' that carry no information. `balanced_joint_effect` is the most affected,
#' since `weakest_effect` at least follows the weaker side once the constant
#' lies off center. The `degenerate` flag and the permutation tests address
#' this case: permuting a constant outcome reproduces the same geometry each
#' time, so every p-value is 1.
#'
#' Dul's conventional bands apply to `d_nec` and `d_suf`, which are ordinary
#' NCA effect sizes; [nsca_magnitude()] returns them as a convention to cite
#' rather than as values calibrated for a particular design. We suggest
#' reading a band together with `d_nec_null`, which indicates how much of the
#' value the geometry produces at that sample size, scope and technique.
#'
#' No benchmarks are supplied for the joint indices. Conventions for a single
#' NCA effect size do not carry over to them: the scales differ, and under
#' independence both components carry a positive finite-sample floor that
#' `balanced_joint_effect` amplifies rather than cancels, by an amount
#' depending on the sample size, the frontier technique and the scope. The
#' joint indices are this package's own constructions, no convention has been
#' published for them, and one supplied here would come from the same package
#' that defines the indices. Calibration by simulation on the design at hand
#' is suggested instead; [nsca_power()] and [nsca_resample()] support this.
#'
#' @section Multiple testing:
#' Each condition involves two component tests, combined into one decision by
#' the intersection-union test, so a model with `k` conditions produces
#' `2k` component p-values and `k` joint ones.
#'
#' Within a condition, the multiplicity is already accounted for. The IUT
#' rejects only when both components reject, so `p_nsca_iut` is the larger of
#' the two p-values and the level of the test is at most `alpha`; the test is
#' conservative, and its rejection rate is usually below the nominal rate.
#'
#' Across conditions, no adjustment is made by default. If the `k` joint tests
#' were independent, the probability of at least one false positive at
#' `alpha` would be at most \eqn{1 - (1 - \alpha)^k}; it is an upper bound
#' because each joint test is itself conservative. The `Settings` block
#' printed by [print()] reports both points together, since the
#' conservativeness of each joint test alone gives an incomplete picture.
#'
#' Adjustment is available as an argument rather than applied by default,
#' for the reason given in [SCAtools::sca_table()]: the appropriate correction
#' depends on the design. When the family is known, `p.adjust` applies the
#' correction; it adjusts the two components and rebuilds the joint p-value
#' from them, as described under that argument.
#'
#' Adjustment changes the p-values and the decisions, not the `*_critical`
#' columns. These remain the critical values of single tests at `alpha`:
#' reference lines rather than the threshold an adjusted test needs to
#' exceed, so an effect size above its critical value may still fall short
#' after adjustment. An adjusted counterpart exists only for Bonferroni, as the
#' `1 - alpha / k` quantile of each null, which requires many more
#' permutations to estimate than the p-value does; the step-wise methods have
#' none, because the level at which each condition is tested depends on the
#' p-values of the others. The table prints a message to this effect whenever
#' it adjusts.
#'
#' Adjustment also raises the smallest p-value a condition can reach. A
#' permutation p-value cannot fall below `1 / test.rep` (`1 / (shared.test.rep
#' + 1)` from the shared sequence), and for a condition that is the only
#' strong one among `k` it is multiplied by at least `k`, so when that product
#' exceeds `alpha` the condition cannot be significant, regardless of the
#' data. The table issues a warning in that case and states the number of
#' permutations that would be required.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param legacy Whether to append the column names retired in 0.3.0, 0.4.0
#'   and 0.5.x as duplicates. See [nsca_legacy_names()].
#' @param p.adjust Multiple-testing correction applied across the conditions
#'   within each frontier technique, as a method name understood by
#'   [stats::p.adjust()]. `"none"`, the default, leaves the table unchanged.
#'
#'   Any other method adjusts each component family separately and rebuilds
#'   the joint p-value from the adjusted pair, so `p_nsca_iut_adjusted` is
#'   `max(p_nec_adjusted, p_suf_adjusted)` rather than
#'   `p.adjust(p_nsca_iut)`. These are different numbers, and only the first
#'   agrees with the component columns printed beside it. Every significance
#'   and support column is then re-derived; the relevance half of
#'   `necessity_supported` and `sufficiency_supported` is unchanged, since an
#'   effect size does not change when more tests are run. The `*_critical`
#'   columns are not adjusted; see the section on multiple testing.
#' @return A data frame of joint results.
#' @seealso [nsca_results()], [nsca_joint()], [nsca_extract()],
#'   [nsca_legacy_names()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_table(fit)
#'
#' # Column names retired in earlier releases, appended as duplicates
#' names(nsca_table(fit, legacy = TRUE))
nsca_table <- function(model, legacy = FALSE, p.adjust = "none") {
  metadata <- .nsca_metadata(model)
  out <- .nsca_joint_rows(model, metadata$conditions, metadata$ceilings)
  .nsca_append_legacy(.nsca_apply_p_adjust(out, metadata, p.adjust), legacy)
}

# The body of the joint table, over a chosen subset.
#
# Every row costs a frontier reconstruction on a 2001-point grid for each side,
# so nsca_extract() asking for one number used to pay for the whole table --
# every condition, every technique -- and to re-emit the frontiers-crossed
# warning each time. print() too. Splitting the loop out lets a caller that
# wants one cell compute one cell.
#
# 
.nsca_joint_rows <- function(model, conditions, ceilings) {
  metadata <- .nsca_metadata(model)
  shared <- metadata$shared
  relevance <- metadata$relevance
  tolerances <- metadata$geometry
  rows <- list()
  index <- 0L
  crossed_any <- character(0L)

  for (condition in conditions) {
    observations <- .nsca_global_value(
      model$necessity, condition, "Number of observations"
    )
    for (ceiling in ceilings) {
      d_nec <- .nsca_nca_param(model$necessity, condition, ceiling, "Effect size")
      d_suf <- .nsca_sca_param(model$sufficiency, condition, ceiling, "Effect size")
      weakest <- suppressWarnings(min(d_nec, d_suf))
      balanced <- nsca_joint(d_nec, d_suf, p = 0)
      component_sum <- d_nec + d_suf
      # The geometric mean over the arithmetic mean: 1 exactly when the two
      # components are equal, and the part of balanced_joint_effect that the
      # name promises but the index alone does not deliver.
      balance <- if (is.finite(component_sum) && component_sum > 0) {
        balanced / component_sum
      } else {
        NA_real_
      }

      band <- .nsca_band_for(metadata, model, condition, ceiling)
      if (isTRUE(band$suspect)) {
        crossed_any <- c(crossed_any, sprintf("%s/%s", condition, ceiling))
      }
      coverage <- if (is.finite(band$share)) 1 - band$share else NA_real_
      reconstruction <- if (is.finite(band$share) && is.finite(band$overlap) &&
                            is.finite(d_nec) && is.finite(d_suf)) {
        abs(1 - band$share - d_nec - d_suf + band$overlap)
      } else {
        NA_real_
      }
      allowance <- .nsca_step_allowance(ceiling, observations)
      excess_overlap <- if (is.finite(band$overlap)) {
        max(0, band$overlap - allowance)
      } else {
        NA_real_
      }

      # The disjointness check of the framework, section 5.2. Whenever the two
      # frontiers do not cross, the floor lies at or below the ceiling
      # throughout the scope, so no point belongs to both expected empty
      # spaces, and because each component effect size is that space's share of
      # the same scope the two cannot sum to more than one. The framework
      # therefore proposes d_nec + d_suf > 1 as a consistency check.
      #
      # That form of the check is a proxy: it infers crossing from the sum
      # because the crossing itself is not otherwise available. Here it is
      # available. overlap_share is the area on which the two empty spaces
      # actually overlap, integrated from the reconstructed frontiers, so the
      # crossing is measured rather than inferred. Reporting component_sum
      # leaves the framework's own rule applicable by any reader; deciding
      # frontiers_crossed from the measurement avoids the proxy's failure mode,
      # which is that the sum and the idealized allowance are computed by
      # different routes and can disagree in the last few digits exactly on the
      # perfect diagonal the check is meant to bless.
      #
      # The step allowance is already subtracted in excess_overlap. A step
      # frontier overlaps itself by construction: between two consecutive
      # observations the upper staircase still holds the earlier value while
      # the lower one has already moved, so the two spaces overlap on every
      # tread, exactly 1 / (n - 1) of the scope summed over the n - 1 treads.
      # For every smooth frontier the allowance is zero and any overlap at all
      # is a genuine crossing.
      frontiers_crossed <- if (is.finite(excess_overlap)) {
        excess_overlap > .nsca_sum_tolerance
      } else {
        NA
      }

      degenerate <- isTRUE(metadata$degeneracy[[condition]]$degenerate)
      geometry_acceptable <- !degenerate &&
        is.finite(excess_overlap) &&
        excess_overlap <= tolerances[["max_overlap"]] &&
        is.finite(reconstruction) &&
        reconstruction <= tolerances[["max_reconstruction"]]

      if (is.null(shared)) {
        p_nec <- .nsca_nca_param(model$necessity, condition, ceiling, "p-value")
        p_suf <- .nsca_sca_param(model$sufficiency, condition, ceiling, "p-value")
        p_nec_accuracy <- .nsca_nca_param(
          model$necessity, condition, ceiling, "p-accuracy"
        )
        p_suf_accuracy <- .nsca_sca_param(
          model$sufficiency, condition, ceiling, "p-accuracy"
        )
        p_weakest <- NA_real_
        p_sum <- NA_real_
        p_weakest_accuracy <- NA_real_
        p_sum_accuracy <- NA_real_
        p_source <- if (isTRUE(metadata$test.rep > 0)) "engine" else "none"
      } else {
        p_nec <- .nsca_permutation_p(shared$necessity[, condition, ceiling], d_nec)
        p_suf <- .nsca_permutation_p(shared$sufficiency[, condition, ceiling], d_suf)
        p_weakest <- .nsca_permutation_p(
          shared$weakest[, condition, ceiling], weakest
        )
        p_sum <- .nsca_permutation_p(
          shared$sum[, condition, ceiling], component_sum
        )
        reps <- dim(shared$necessity)[[1L]]
        exhaustive <- isTRUE(attr(shared, "exhaustive"))
        accuracy <- function(p) {
          .nsca_p_accuracy(p, reps, metadata$p_confidence, exhaustive)
        }
        p_nec_accuracy <- accuracy(p_nec)
        p_suf_accuracy <- accuracy(p_suf)
        p_weakest_accuracy <- accuracy(p_weakest)
        p_sum_accuracy <- accuracy(p_sum)
        p_source <- "shared"
      }

      # Every null summary is read from the distribution the p-values beside
      # it came from, so a floor, a critical value and a p-value in one row
      # always describe the same shuffles.
      null_nec <- .nsca_null_draws(model, shared, "necessity", condition, ceiling)
      null_suf <- .nsca_null_draws(model, shared, "sufficiency", condition, ceiling)
      null_weakest <- .nsca_null_draws(model, shared, "weakest", condition, ceiling)
      null_sum <- .nsca_null_draws(model, shared, "sum", condition, ceiling)
      tested <- !identical(p_source, "none")
      p_iut <- if (is.finite(p_nec) && is.finite(p_suf)) {
        max(p_nec, p_suf)
      } else {
        NA_real_
      }

      nec_significant <- if (tested && is.finite(p_nec)) {
        p_nec <= metadata$alpha
      } else {
        NA
      }
      suf_significant <- if (tested && is.finite(p_suf)) {
        p_suf <= metadata$alpha
      } else {
        NA
      }
      nsca_significant <- if (is.na(nec_significant) || is.na(suf_significant)) {
        NA
      } else {
        nec_significant && suf_significant
      }

      relevant <- function(value, threshold) {
        if (is.na(threshold)) {
          return(TRUE)
        }
        isTRUE(is.finite(value) && value >= threshold)
      }
      nec_supported <- if (is.na(nec_significant)) {
        NA
      } else {
        nec_significant && relevant(d_nec, relevance[["necessity"]])
      }
      suf_supported <- if (is.na(suf_significant)) {
        NA
      } else {
        suf_significant && relevant(d_suf, relevance[["sufficiency"]])
      }
      joint_support <- if (is.na(nec_supported) || is.na(suf_supported)) {
        NA
      } else {
        nec_supported && suf_supported && geometry_acceptable
      }

      index <- index + 1L
      rows[[index]] <- data.frame(
        condition = condition,
        outcome = metadata$outcome,
        direction = unname(metadata$direction[[condition]]),
        equivalence_class = unname(
          .nsca_equivalence_class[[metadata$direction[[condition]]]]
        ),
        relationship = .nsca_relationship_text(
          unname(metadata$direction[[condition]]), condition, metadata$outcome
        ),
        statement = unname(metadata$statement[[condition]]),
        frontier = ceiling,
        necessity_corner = unname(metadata$necessity_corner[[condition]]),
        necessity_empty_space = .nsca_empty_space_text(
          unname(metadata$necessity_corner[[condition]]),
          condition, metadata$outcome
        ),
        necessity_boundary = .nsca_boundary_type(
          unname(metadata$necessity_corner[[condition]])
        ),
        sufficiency_corner = unname(metadata$sufficiency_corner[[condition]]),
        sufficiency_empty_space = .nsca_empty_space_text(
          unname(metadata$sufficiency_corner[[condition]]),
          condition, metadata$outcome
        ),
        sufficiency_boundary = .nsca_boundary_type(
          unname(metadata$sufficiency_corner[[condition]])
        ),
        observations = observations,
        # The denominator both components share. d_nec and d_suf are each an
        # empty area over this one number, which is why they are on the same
        # scale and why every joint index below can combine them at all.
        # scope_adjusted says whether it is the scope as declared or the scope
        # the engines widened to hold the data; see nsca_scope().
        scope_area = .nsca_scope_area(metadata, condition),
        scope_adjusted = .nsca_scope_is_adjusted(metadata, condition),
        d_nec = d_nec,
        d_suf = d_suf,
        # What random pairing alone produces on each side at this n, scope
        # and technique, and how much of each component it does not account
        # for. Both frontiers are built from sample extremes, so both are
        # overstated, and the joint indices multiply the two overstatements
        # rather than canceling them -- which is exactly why the two floors
        # are printed beside the two components and not only in the docs.
        d_nec_null = .nsca_draw_mean(null_nec),
        d_suf_null = .nsca_draw_mean(null_suf),
        d_nec_excess = .nsca_excess(d_nec, .nsca_draw_mean(null_nec)),
        d_suf_excess = .nsca_excess(d_suf, .nsca_draw_mean(null_suf)),
        # The bar each statistic has to clear: the 1 - alpha quantile of its
        # own null distribution. Each statistic has its own; the critical
        # value of the sum is not the sum of the two critical values.
        d_nec_critical = .nsca_critical(null_nec, metadata$alpha),
        d_suf_critical = .nsca_critical(null_suf, metadata$alpha),
        component_sum = component_sum,
        component_sum_critical = .nsca_critical(null_sum, metadata$alpha),
        # Two points on one power-mean family, p = -Inf and p = 0. The minimum
        # is left as a literal min() rather than routed through nsca_joint():
        # the two agree on every value an effect size can take, but nsca_joint()
        # returns NA for a negative input where min() returns the negative
        # number, and an engine that ever reported a negative area should
        # surface it here rather than have it become a missing value.
        weakest_effect = weakest,
        weakest_effect_critical = .nsca_critical(null_weakest, metadata$alpha),
        weaker_component = if (!is.finite(d_nec) || !is.finite(d_suf)) {
          NA_character_
        } else if (d_nec <= d_suf) "necessity" else "sufficiency",
        balanced_joint_effect = balanced,
        balance = balance,
        admissible_region_share = band$share,
        joint_empty_space_coverage = coverage,
        overlap_share = band$overlap,
        frontiers_crossed = frontiers_crossed,
        degenerate = degenerate,
        geometry_acceptable = geometry_acceptable,
        reconstruction_error = reconstruction,
        p_source = p_source,
        p_nec = p_nec,
        p_nec_accuracy = p_nec_accuracy,
        p_suf = p_suf,
        p_suf_accuracy = p_suf_accuracy,
        p_nsca_iut = p_iut,
        p_weakest_perm = p_weakest,
        p_weakest_perm_accuracy = p_weakest_accuracy,
        p_sum_perm = p_sum,
        p_sum_perm_accuracy = p_sum_accuracy,
        necessity_significant = nec_significant,
        sufficiency_significant = suf_significant,
        nsca_significant = nsca_significant,
        necessity_supported = nec_supported,
        sufficiency_supported = suf_supported,
        joint_support = joint_support,
        significance_status = .nsca_verdict(p_nec, p_suf, metadata$alpha, tested),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  if (length(crossed_any) > 0L) {
    warning(
      sprintf(
        paste0(
          "The two frontiers cross for %s by more than %g%% of the scope, so ",
          "the expected empty spaces overlap rather than partition the scope, ",
          "and d_nec + d_suf can exceed 1. This is more than the small ",
          "overlap inherent in a step frontier; see frontiers_crossed. It is ",
          "advisable to inspect the plot before interpreting these areas."
        ),
        paste(unique(crossed_any), collapse = ", "),
        100 * tolerances[["max_overlap"]]
      ),
      call. = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  # No legacy handling here. This function builds the rows; the retired-name
  # duplication belongs to the exported wrapper, which is the thing that has a
  # `legacy` argument. The call was left behind when nsca_table() was split so
  # that nsca_extract() could build one row without rebuilding the whole
  # geometry, and it referred to a variable this function never had.
  out
}

# Append retired names as duplicate columns. Only names whose replacement this
# particular table carries are added, so one helper serves both the joint table
# and the dual-threshold table.
.nsca_append_legacy <- function(out, legacy) {
  if (!isTRUE(legacy)) {
    return(out)
  }
  for (old in names(.nsca_legacy_columns)) {
    new <- .nsca_legacy_columns[[old]]
    if (new %in% names(out)) {
      out[[old]] <- out[[new]]
    }
  }
  out
}

#' All three NSCA result layers
#'
#' Returns the necessity analysis, the sufficiency analysis, and their joint
#' necessary-and-sufficient analysis as three explicit tables. The third is
#' not an independent analysis: its `p_nsca_iut` is the conjunction
#' `max(p_nec, p_suf)` and can pass only when both component tests pass.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param object An object returned by [nsca_analysis()].
#' @param x An object returned by [nsca_analysis()], for the print method.
#' @param ... Additional arguments reserved for methods.
#' @return A named list containing `necessity`, `sufficiency`, and
#'   `necessary_and_sufficient` data frames. The methods print, or return a
#'   summary object.
#' @seealso [nsca_table()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' # With the two component tests, so that the joint decision can be reached.
#' # A report would use 1000 or more permutations.
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
#'                      test.rep = 99, test.seed = 1)
#' layers <- nsca_results(fit)
#' names(layers)
#' layers$necessity
#' fit
#' summary(fit)
nsca_results <- function(model) {
  .nsca_metadata(model)
  list(
    necessity = .nsca_component_table(model, "necessity"),
    sufficiency = .nsca_component_table(model, "sufficiency"),
    necessary_and_sufficient = nsca_table(model)
  )
}

#' Extract one value from an NSCA result
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x Condition name. Defaults to the first condition.
#' @param ceiling Frontier technique. Defaults to the first requested.
#' @param param A column of [nsca_table()], or a parameter name understood by
#'   the underlying engines, prefixed with `"nec:"` or `"suf:"`. Retired column
#'   names, listed by [nsca_legacy_names()], are rejected with an error that
#'   names their replacement.
#' @return The selected value.
#' @seealso [nsca_table()], [nsca_legacy_names()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_extract(fit, param = "weakest_effect")
#' nsca_extract(fit, param = "nec:Ceiling accuracy")
#' nsca_extract(fit, param = "suf:Effect size")
nsca_extract <- function(model, x = NULL, ceiling = NULL,
                         param = "weakest_effect") {
  metadata <- .nsca_metadata(model)
  if (is.null(x)) {
    x <- metadata$conditions[[1L]]
  }
  if (!(x %in% metadata$conditions)) {
    stop(sprintf("Unknown condition '%s'.", x), call. = FALSE)
  }
  if (is.null(ceiling)) {
    ceiling <- metadata$ceilings[[1L]]
  }

  # Checked before any route rather than left to the row builder. The `nec:`
  # and `suf:` routes read the engines directly and returned NA for a frontier
  # that was never estimated, which reads as a missing measurement rather
  # than a wrong request; the joint route passes whatever the caller asked
  # for, where nsca_table() only ever passes the techniques the model holds.
  if (length(ceiling) != 1L || !(ceiling %in% metadata$ceilings)) {
    stop(
      sprintf(
        "Frontier '%s' was not estimated. Available: %s.",
        paste(ceiling, collapse = ", "), paste(metadata$ceilings, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (grepl("^nec:", param)) {
    return(.nsca_nca_param(model$necessity, x, ceiling, sub("^nec:", "", param)))
  }
  if (grepl("^suf:", param)) {
    return(.nsca_sca_param(model$sufficiency, x, ceiling, sub("^suf:", "", param)))
  }

  if (param %in% names(.nsca_legacy_columns)) {
    .nsca_retired_stop(param)
  }

  # One cell, not the whole table. The geometry behind a row costs a 2001-point
  # reconstruction per side, and the frontiers-crossed warning belongs to
  # nsca_table(), which the caller can run when they want the table.
  row <- suppressWarnings(.nsca_joint_rows(model, x, ceiling))
  if (is.null(row) || nrow(row) == 0L) {
    stop(sprintf("No result for condition '%s' with frontier '%s'.", x, ceiling),
         call. = FALSE)
  }
  if (!(param %in% names(row))) {
    # A dual-threshold column is defined per outcome level, not per condition
    # and frontier, so there is no single value to return. Say that rather than
    # reporting the name as unknown, which it is not.
    if (param %in% .nsca_threshold_columns) {
      stop(
        sprintf(
          paste0(
            "'%s' is a column of nsca_thresholds(), which has one row per ",
            "outcome level rather than one value per condition and frontier. ",
            "The column can be selected from the output of nsca_thresholds()."
          ),
          param
        ),
        call. = FALSE
      )
    }
    stop(
      sprintf(
        paste0(
          "Unknown parameter '%s'. Expected a column of nsca_table(), or an ",
          "engine parameter prefixed with 'nec:' or 'suf:'."
        ),
        param
      ),
      call. = FALSE
    )
  }
  row[[param]][[1L]]
}

#' @rdname nsca_results
#' @export
print.nsca_result <- function(x, ...) {
  metadata <- .nsca_metadata(x)
  tables <- nsca_results(x)
  joint <- tables$necessary_and_sufficient
  cat("\nNecessary and Sufficient Condition Analysis (NSCA)\n")
  cat("Engines: NCA ", metadata$engines[["NCA"]],
      " (necessity), SCAtools ", metadata$engines[["SCAtools"]],
      " (sufficiency)\n", sep = "")

  mapping <- data.frame(
    condition = metadata$conditions,
    direction = unname(metadata$direction[metadata$conditions]),
    nec_corner = unname(metadata$necessity_corner[metadata$conditions]),
    suf_corner = unname(metadata$sufficiency_corner[metadata$conditions]),
    statement = unname(metadata$statement[metadata$conditions]),
    stringsAsFactors = FALSE
  )
  print(mapping, row.names = FALSE)

  cat("\n1. Necessary Condition Analysis\n")
  print(
    tables$necessity[c(
      "condition", "frontier", "effect_size", "p_value",
      "frontier_accuracy", "fit"
    )],
    row.names = FALSE
  )

  cat("\n2. Sufficient Condition Analysis\n")
  print(
    tables$sufficiency[c(
      "condition", "frontier", "effect_size", "p_value",
      "frontier_accuracy", "fit"
    )],
    row.names = FALSE
  )

  cat("\n3. Necessary and Sufficient Condition Analysis\n")
  cat("Three summaries combining the same two components; none replaces them.\n")
  cat("  weakest_effect            = min(d_nec, d_suf)       non-compensatory, max ~0.5 + O/2\n")
  cat("  balanced_joint_effect     = 2*sqrt(d_nec*d_suf)     partially compensatory, max ~1 + O\n")
  cat("  joint_empty_space_coverage = 1 - admissible_region_share     curvature-neutral, within [0, 1]\n")
  print(
    joint[c(
      "condition", "frontier", "d_nec", "d_suf", "component_sum",
      "weakest_effect", "balanced_joint_effect", "balance",
      "admissible_region_share", "joint_empty_space_coverage",
      "overlap_share", "frontiers_crossed"
    )],
    row.names = FALSE
  )
  # Split, because the two halves of that table are in different positions.
  # d_nec and d_suf are ordinary NCA effect sizes and Dul's bands are a
  # convention to cite for them -- the package was already applying his
  # lower edge as the default d_cut, so refusing to name the bands was
  # disclosing less than it was already doing. The joint indices are this
  # package's own constructions and nobody has published a convention for
  # them; one supplied here would be invented by the same package that
  # defines the index.
  cat(
    "\nd_nec and d_suf are ordinary NCA effect sizes; nsca_magnitude()",
    "\nreturns the conventional bands of Dul (2016) for them, to be cited",
    "\nrather than treated as calibrated. A band is best read beside",
    "\nd_nec_null, which indicates how much of the value the geometry",
    "\nproduces at this n, scope and technique.",
    "\n\nThe joint indices have no such convention and no established",
    "\nbenchmarks; they are this package's own constructions. We suggest",
    "\nreporting them alongside both components; it is also advisable to",
    "\ncalibrate by simulation before naming a value large or small.\n",
    sep = ""
  )

  cat("\n4. Decision\n")
  cat("Each p-value addresses a different question; joint_support uses p_nsca_iut.\n")
  cat("  p_nsca_iut      each corner is emptier than random pairing: max(p_nec, p_suf)\n")
  if (identical(joint$p_source[[1L]], "shared")) {
    cat("  p_weakest_perm  the smaller empty area is larger than random pairing produces\n")
    cat("  p_sum_perm      the two corners together are emptier than random pairing\n")
  }
  cat("Support also requires any pre-specified relevance threshold to be met\n")
  cat("and the geometry to be acceptable.\n")
  print(
    joint[c(
      "condition", "frontier", "p_source", "p_nec", "p_suf", "p_nsca_iut",
      "p_weakest_perm", "p_sum_perm", "necessity_supported",
      "sufficiency_supported", "geometry_acceptable", "joint_support"
    )],
    row.names = FALSE
  )

  if (identical(joint$p_source[[1L]], "none")) {
    cat("\nNo permutation test was run, so no joint-support decision is available.\n")
    cat("Setting test.rep tests both components; joint support requires both.\n")
  }
  if (all(is.na(metadata$relevance))) {
    cat("\nNo practical-relevance threshold was pre-specified, so support rests\n")
    cat("on significance alone. Setting 'relevance' adds an effect-size criterion.\n")
  }
  if (any(joint$degenerate)) {
    cat(
      "\nDegenerate row(s) present: an axis does not vary inside the declared",
      "\nscope, so every index above reflects the scope rather than the data,",
      "\nand overstates the evidence. The p-values are more informative here.\n",
      sep = ""
    )
  }
  bad <- joint$reconstruction_error >
    metadata$geometry[["max_reconstruction"]]
  if (isTRUE(any(bad, na.rm = TRUE))) {
    cat(
      "\nThe area identity does not close for at least one row. The frontiers",
      "\nreconstructed by this package differ from those the engines fitted,",
      "\nso the admissible region and the plots may be inaccurate even if they",
      "\nlook plausible; see reconstruction_error.\n",
      sep = ""
    )
  }
  .nsca_settings_block(metadata)

  if (isTRUE(attr(x, "show.plots"))) {
    plot(x)
  }
  cat(
    "
Inference note: component and joint p-values screen against the\n",
    "random-pairing/independence null; they do not by themselves establish\n",
    "a causal or deterministic necessary-and-sufficient relation.\n",
    sep = ""
  )
  invisible(x)
}

#' @rdname nsca_results
#' @export
summary.nsca_result <- function(object, ...) {
  metadata <- .nsca_metadata(object)
  tables <- nsca_results(object)
  out <- list(
    call = metadata$call,
    map = nsca_direction_map(),
    necessity = tables$necessity,
    sufficiency = tables$sufficiency,
    necessary_and_sufficient = tables$necessary_and_sufficient,
    results = tables$necessary_and_sufficient,
    alpha = metadata$alpha,
    relevance = metadata$relevance,
    test.rep = metadata$test.rep,
    shared.test.rep = metadata$shared.test.rep
  )
  class(out) <- "summary_nsca_result"
  out
}

#' @rdname nsca_results
#' @export
print.summary_nsca_result <- function(x, ...) {
  cat("\nNecessary and Sufficient Condition Analysis summary\n\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("1. Necessary Condition Analysis:\n")
  print(x$necessity, row.names = FALSE)
  cat("\n2. Sufficient Condition Analysis:\n")
  print(x$sufficiency, row.names = FALSE)
  cat("\n3. Necessary and Sufficient Condition Analysis:\n")
  print(x$necessary_and_sufficient, row.names = FALSE)
  cat("\nThe joint-support screen uses p_nsca_iut = max(p_nec, p_suf) at\nalpha = ",
      format(x$alpha), ", so both directional component tests need to reject.\n",
      sep = "")
  cat("This is not a direct test of a causal necessary-and-sufficient relation.\n")
  cat("weakest_effect, balanced_joint_effect and joint_empty_space_coverage are\n")
  cat("three summaries of the same two components, differing in how much\n")
  cat("compensation they allow and in how they treat curvature. None has\n")
  cat("established magnitude benchmarks.\n")
  invisible(x)
}

# Retired names are no longer usable. They stop with the replacement named,
# because a name that quietly returns NULL turns a working script into a wrong
# one rather than a broken one. The column duplicates behind `legacy = TRUE`
# remain: that mechanism copies the current column rather than translating a
# name, so it cannot give a different answer.
.nsca_retired_stop <- function(name) {
  replacement <- .nsca_legacy_columns[[name]]
  stop(
    sprintf(
      paste0(
        "'%s' was renamed to '%s' and is no longer accepted here. Retired ",
        "column names are still available as duplicates from nsca_table(legacy ",
        "= TRUE) and nsca_thresholds(legacy = TRUE); see nsca_legacy_names()."
      ),
      name, replacement
    ),
    call. = FALSE
  )
}

# Same convention as SCAtools: every share is reported on [0, 1]. Defined here
# rather than reached for with ::: , which R CMD check flags and which would
# tie NSCA to another package's internals.
.nsca_as_proportion <- function(value) {
  if (!is.finite(value)) {
    return(NA_real_)
  }
  value / 100
}

# What the analysis was actually run with, and which of it the caller chose.
#
# `nsca()` takes four arguments and reports two engines, three joint indices,
# an admissible region, a permutation test and a geometry check. Almost every
# number it prints therefore rests on a choice the package made: which
# frontiers, how many replications, what alpha, which scope, what counts as
# acceptable geometry. Nothing in the tables distinguishes those from choices
# the researcher made, and a methods section has to. The recorded call names
# what was passed, so everything absent from it is a default.
#
# Deliberately kept in the same shape and order as SCAtools' block. The two
# packages are read side by side, and a reader comparing a d_suf here with an
# effect size there needs to see at a glance that the same conventions
# produced both.
.nsca_supplied <- function(metadata) {
  call <- metadata$call
  if (is.null(call)) {
    return(NULL)
  }
  setdiff(names(as.list(call)), c("", "data", "x", "y"))
}

.nsca_settings_block <- function(metadata) {
  supplied <- .nsca_supplied(metadata)
  # NULL is "cannot tell", not "nothing was supplied". Print no marks rather
  # than wrong ones.
  mark <- if (is.null(supplied)) {
    function(...) ""
  } else {
    function(...) if (any(c(...) %in% supplied)) "" else "  (default)"
  }

  scope_text <- if (identical(metadata$scope_source, "observed")) {
    "observed range of each condition and the outcome"
  } else {
    paste0(
      "declared: ",
      paste(format(unlist(metadata$scope_declared)), collapse = ", ")
    )
  }

  relevance <- metadata$relevance
  relevance_text <- if (is.null(relevance) || all(is.na(relevance))) {
    "none pre-specified"
  } else if (length(relevance) == 2L) {
    sprintf("d_nec >= %s, d_suf >= %s",
            format(relevance[[1L]]), format(relevance[[2L]]))
  } else {
    paste(format(relevance), collapse = ", ")
  }

  cat("\nSettings\n")
  cat(sprintf("  frontiers       %s%s\n",
              paste(metadata$ceilings, collapse = ", "), mark("ceilings")))
  if (length(metadata$reference) > 0L) {
    cat(sprintf("  reference       %s%s\n",
                paste(metadata$reference, collapse = ", "), mark("reference")))
    if ("qr" %in% metadata$reference) {
      cat(sprintf("  reference.tau   %s%s\n",
                  format(if (is.null(metadata$reference.tau)) 0.5
                         else metadata$reference.tau),
                  mark("reference.tau")))
    }
  }
  cat(sprintf("  test.rep        %s%s\n",
              format(metadata$test.rep), mark("test.rep")))
  cat(sprintf("  test.seed       %s%s\n",
              if (is.null(metadata$test.seed)) "not fixed" else
                format(metadata$test.seed),
              mark("test.seed")))
  if (isTRUE(metadata$test.rep > 0) && is.null(metadata$test.seed)) {
    cat("                the p-values above may differ on the next run
")
  }
  cat(sprintf("  shared.test.rep %s%s\n",
              format(metadata$shared.test.rep), mark("shared.test.rep")))
  cat(sprintf("  shared.seed     %s%s\n",
              if (is.null(metadata$shared.seed)) "not fixed" else
                format(metadata$shared.seed),
              mark("shared.seed")))
  cat(sprintf("  alpha           %s%s\n",
              format(metadata$alpha), mark("test.p_threshold")))
  .nsca_multiplicity_block(metadata)

  cat(sprintf("  relevance       %s%s\n", relevance_text, mark("relevance")))
  cat(sprintf("  scope           %s%s\n", scope_text, mark("scope")))
  # Not a default and not a choice: something the engines did to the choice.
  # The declared scope is what a methods section reports, the widened one is
  # what every area above was divided by.
  adjustment <- metadata$scope_adjustment
  if (!is.null(adjustment) && nrow(adjustment) > 0L) {
    cat(sprintf(
      "                widened by the engines to hold every observation (%s); see nsca_scope()\n",
      paste(unique(adjustment$condition), collapse = ", ")
    ))
  }
  # Not a setting either, and for the same reason as the line above: it says
  # whether the areas describe the data or describe the scope. Printed only
  # when an axis is thin, so an ordinary analysis shows nothing here.
  degeneracy <- metadata$degeneracy
  if (length(degeneracy) > 0L) {
    worst <- vapply(degeneracy, function(item) {
      min(c(item$x_coverage, item$y_coverage))
    }, numeric(1L))
    thin <- is.finite(worst) & worst < 0.5
    if (any(thin)) {
      cat(sprintf(
        "                observations span %s of it\n",
        paste(sprintf("%s %.0f%%", names(worst)[thin], 100 * worst[thin]),
              collapse = ", ")
      ))
    }
  }
  cat(sprintf("  na.action       %s%s\n",
              metadata$na.action, mark("na.action")))
  cat(sprintf("  thresholds      x on %s, y on %s%s\n",
              metadata$threshold.x, metadata$threshold.y,
              mark("threshold.x", "threshold.y")))
  cat(sprintf("  convention      %s%s\n",
              metadata$convention, mark("convention")))
  # The two geometry tolerances decide geometry_acceptable, which decides
  # joint_support. A verdict that depends on a number the reader has never
  # seen is not a reportable verdict.
  # Zero is the default and means "no tolerance declared", which is a
  # reportable fact rather than an absence: it says the correspondence column
  # marks only exact zeros.
  cat(sprintf("  correspondence  %s%s\n",
              .nsca_tolerance_text(metadata$correspondence.tolerance),
              mark("correspondence.tolerance")))
  cat(sprintf("  geometry        max_overlap %s, max_reconstruction %s%s\n",
              format(metadata$geometry[["max_overlap"]]),
              format(metadata$geometry[["max_reconstruction"]]),
              mark("geometry.max_overlap", "geometry.max_reconstruction")))
  invisible(NULL)
}

# One line for a tolerance that is one number per condition.
#
# Stored per condition because it is a width in the units of X and the
# conditions need not share a scale. Printed as one value when they all agree,
# which is the common case and the only one where a bare number is not
# ambiguous.
.nsca_tolerance_text <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return("0")
  }
  if (length(unique(value)) == 1L) {
    return(format(unname(value[[1L]])))
  }
  paste(sprintf("%s %s", names(value), format(unname(value))), collapse = ", ")
}

# The counterpart of SCAtools' `.sca_multiplicity_block()`, and different in
# the one way that matters: a condition here costs two component tests, not
# one, and they are combined before anything is decided.
#
# That combination is already conservative. The intersection-union test
# rejects only when both components reject, so its level is at most alpha and
# usually below it -- the within-condition multiplicity is not merely
# disclosed, it is paid for by construction. Across conditions nothing is
# paid, which is why the two facts are printed together: a reader who sees
# "conservative" without seeing "across conditions, not adjusted" has been
# told the reassuring half.
#
# Disclosure rather than adjustment, for the same reason as the other package:
# which correction fits depends on what the family is, and that is a property
# of the design rather than of the data.
.nsca_multiplicity_block <- function(metadata) {
  tested <- isTRUE(metadata$test.rep > 0) ||
    isTRUE(metadata$shared.test.rep > 0)
  if (!tested) {
    return(invisible(NULL))
  }
  k <- length(metadata$conditions)
  m <- length(metadata$ceilings)
  alpha <- suppressWarnings(as.numeric(metadata$alpha))

  if (k > 1L) {
    cat(sprintf(
      "  multiplicity  %d conditions%s: %d component p-values, %d joint tests\n",
      k,
      if (m > 1L) sprintf(" x %d techniques", m) else "",
      2L * k * m, k * m
    ))
    cat("                the intersection-union test is conservative within a\n")
    cat("                condition; nothing is adjusted across conditions\n")
    if (isTRUE(is.finite(alpha)) && alpha > 0 && alpha < 1) {
      cat(sprintf(
        paste0("                were the %d joint tests independent, P(at least one",
               "\n                false positive) at alpha %s would be at most %.2f\n"),
        k, format(alpha), 1 - (1 - alpha)^k
      ))
    }
    .nsca_supplementary_line(metadata, k * m)
  } else {
    cat(sprintf(
      "  multiplicity  1 condition%s: %d component p-values, %d joint test%s\n",
      if (m > 1L) sprintf(" x %d techniques", m) else "",
      2L * m, m, if (m > 1L) "s" else ""
    ))
    cat("                the intersection-union test is conservative: it rejects\n")
    cat("                only when both components do\n")
    if (m > 1L) {
      cat("                it is advisable to report the technique chosen in advance\n                rather than the one with the smallest p\n")
    }
    .nsca_supplementary_line(metadata, m)
  }
  invisible(NULL)
}

# The supplementary tests a shared sequence adds, listed apart from the joint
# decisions. They use the same shuffles as the component tests and answer
# different questions; reported as findings alongside the joint decision they
# add to the number of chances at a false positive, so the line says which
# one to choose and when.
.nsca_supplementary_line <- function(metadata, joint_tests) {
  if (!isTRUE(metadata$shared.test.rep > 0)) {
    return(invisible(NULL))
  }
  cat(sprintf(
    paste0(
      "                plus %d supplementary tests on the same shuffles",
      " (p_weakest_perm,\n                p_sum_perm); we suggest choosing the",
      " p-value that answers the\n                research question before",
      " reading the results\n"
    ),
    2L * joint_tests
  ))
  invisible(NULL)
}

# Multiple-testing adjustment for the joint table.
#
# The counterpart of SCAtools' `.sca_apply_p_adjust()`, and different in the
# one way that matters: a condition here is decided by two component tests
# combined, so there are two families, not one.
#
# Each component's p-values are adjusted over the conditions, and the joint
# p-value is then rebuilt from the adjusted pair. It is *not* obtained by
# handing `p_nsca_iut` to `stats::p.adjust()`, which would be a different and
# wrong number: the intersection-union test rejects only when both components
# reject, so the coherent adjusted decision is
# `max(p_nec_adjusted, p_suf_adjusted) <= alpha`, exactly as the unadjusted one
# is `max(p_nec, p_suf) <= alpha`. Adjusting the maximum instead would control
# nothing in particular and would not agree with the component columns printed
# beside it.
#
# Controlling each component family at alpha leaves the joint decision at most
# alpha, since it is the more conservative of the two. That is the same
# direction the IUT was already conservative in, so the adjustment costs power
# and never buys false rejections.
#
# As in SCAtools, the family is the conditions within one frontier technique.
.nsca_adjust_within <- function(p, group, method) {
  out <- rep(NA_real_, length(p))
  for (g in unique(group)) {
    i <- group == g
    out[i] <- stats::p.adjust(p[i], method = method)
  }
  out
}

.nsca_apply_p_adjust <- function(out, metadata, method) {
  method <- .nsca_validate_p_adjust(method)
  if (identical(method, "none") || nrow(out) == 0L) {
    return(out)
  }

  # The floor of the p-values the table actually holds: the shared sequence's
  # when there is one, the engines' otherwise.
  if (is.null(metadata$shared)) {
    .nsca_warn_adjust_floor(
      length(unique(out$condition)), metadata$test.rep, metadata$alpha, method,
      if (is.null(metadata$observations)) Inf else metadata$observations
    )
  } else {
    .nsca_warn_adjust_floor(
      length(unique(out$condition)), dim(metadata$shared$necessity)[[1L]],
      metadata$alpha, method, shared = TRUE
    )
  }
  group <- out$frontier
  nec <- .nsca_adjust_within(out$p_nec, group, method)
  suf <- .nsca_adjust_within(out$p_suf, group, method)
  iut <- suppressWarnings(pmax(nec, suf))
  weakest <- .nsca_adjust_within(out$p_weakest_perm, group, method)
  total <- .nsca_adjust_within(out$p_sum_perm, group, method)

  out <- .nsca_insert_after(out, "p_nec_accuracy", "p_nec_adjusted", nec)
  out <- .nsca_insert_after(out, "p_suf_accuracy", "p_suf_adjusted", suf)
  out <- .nsca_insert_after(out, "p_nsca_iut", "p_nsca_iut_adjusted", iut)
  out <- .nsca_insert_after(out, "p_weakest_perm_accuracy",
                            "p_weakest_perm_adjusted", weakest)
  out <- .nsca_insert_after(out, "p_sum_perm_accuracy", "p_sum_perm_adjusted",
                            total)
  out <- .nsca_insert_after(out, "p_sum_perm_adjusted", "p_adjust_method",
                            method)

  # Every decision re-derived from the adjusted pair. `necessity_supported`
  # and `sufficiency_supported` keep their relevance half untouched -- an
  # effect size does not become smaller because more tests were run -- so the
  # relevance verdict is read back off the existing column rather than
  # recomputed.
  alpha <- metadata$alpha
  tested <- is.finite(nec) & is.finite(suf)

  nec_sig <- ifelse(is.finite(nec), nec <= alpha, NA)
  suf_sig <- ifelse(is.finite(suf), suf <= alpha, NA)

  nec_relevant <- .nsca_relevance_half(out$necessity_supported,
                                       out$necessity_significant)
  suf_relevant <- .nsca_relevance_half(out$sufficiency_supported,
                                       out$sufficiency_significant)

  out$necessity_significant <- nec_sig
  out$sufficiency_significant <- suf_sig
  out$nsca_significant <- ifelse(is.na(nec_sig) | is.na(suf_sig), NA,
                                 nec_sig & suf_sig)
  out$necessity_supported <- ifelse(is.na(nec_sig), NA, nec_sig & nec_relevant)
  out$sufficiency_supported <- ifelse(is.na(suf_sig), NA,
                                      suf_sig & suf_relevant)
  out$joint_support <- ifelse(
    is.na(out$necessity_supported) | is.na(out$sufficiency_supported), NA,
    out$necessity_supported & out$sufficiency_supported
  )
  out$significance_status <- vapply(
    seq_len(nrow(out)),
    function(i) .nsca_verdict(nec[[i]], suf[[i]], alpha, tested[[i]]),
    character(1L)
  )

  # The critical values stay where they were, and nothing in the table can say
  # so, so the output does. An adjusted counterpart exists only for Bonferroni
  # -- the 1 - alpha/k quantile of each null, which takes far more permutations
  # to estimate than the p-value does -- and not at all for the step-wise
  # methods, where the level a condition is tested at depends on the p-values
  # of the others.
  if (any(tested)) {
    message(sprintf(
      paste0(
        "p-values adjusted with '%s' across the conditions within each ",
        "frontier; every significance and support column, including ",
        "joint_support, follows the adjusted p-values. The *_critical ",
        "columns are not adjusted: they remain single-test reference lines ",
        "at alpha = %s."
      ),
      method, format(alpha)
    ))
  }
  out
}

# Recover "was this component relevant?" from the two columns already in the
# table, rather than re-reading `relevance` and the effect size.
#
# `supported` is `significant && relevant`, so where the component was
# significant the two columns differ exactly when relevance failed. Where it
# was not significant, `supported` is FALSE whatever relevance said and the
# question cannot be answered from these two columns -- but it does not need
# to be, because an adjusted p can only be larger, so a component that was not
# significant before cannot become significant now and its `supported` stays
# FALSE either way.
.nsca_relevance_half <- function(supported, significant) {
  ifelse(is.na(significant), NA,
         ifelse(significant, supported, TRUE))
}

.nsca_insert_after <- function(out, after, name, value) {
  out[[name]] <- value
  keep <- setdiff(names(out), name)
  before <- keep[seq_len(match(after, keep))]
  out[c(before, name, setdiff(keep, before))]
}

# Mirrors SCAtools' validator rather than reaching for it: it is private
# there, and a four-line check is not worth an export. The message lists what
# `stats::p.adjust` accepts, so the two can never disagree about the set.
.nsca_validate_p_adjust <- function(method) {
  method <- as.character(method)[[1L]]
  if (!method %in% stats::p.adjust.methods) {
    stop(
      sprintf(
        "Unknown p.adjust method '%s'. Available: %s.",
        method, paste(stats::p.adjust.methods, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  method
}

# What a random pairing of the same data produces, per component.
#
# The counterpart of SCAtools' `.sca_null_mean()`, and needed separately
# because that one is private and because NSCA reads two engines. Both hold
# an NCA-shaped `tests` list, so one reader serves both sides.
#
# Both frontiers are fitted to sample extremes and both are therefore
# overstated. That matters more here than in a one-sided analysis: the joint
# indices are built from the two components, so a small sample inflates both
# and the product, geometric mean and coverage inherit the inflation rather
# than canceling it. The two floors are printed beside the two components so
# that a reader can see how much of a joint index is sample size.
#
# Reported, never subtracted: the difference is a comparison, not a corrected
# estimator, since the bias under independence is not the bias under the
# alternative.
.nsca_null_mean <- function(fit, condition, ceiling) {
  item <- fit$tests[[condition]]
  if (is.null(item) || is.null(item[[ceiling]])) {
    return(NA_real_)
  }
  draws <- suppressWarnings(as.numeric(item[[ceiling]]$data))
  draws <- draws[is.finite(draws)]
  if (length(draws) == 0L) {
    return(NA_real_)
  }
  mean(draws)
}

.nsca_excess <- function(observed, null_mean) {
  if (!is.finite(observed) || !is.finite(null_mean)) {
    return(NA_real_)
  }
  observed - null_mean
}

# The smallest p-value a permutation test in this package can return.
#
# The engines report max((count + 1) / (R + 1), 1 / R), so nothing below
# 1 / R; the shared sequence reports (count + 1) / (R + 1), so nothing below
# 1 / (R + 1). Both cap R at the n! orderings a sample of n has. A test whose
# smallest attainable p-value is above alpha rejects nothing, whatever the
# data, and its silence looks exactly like a null result.
.nsca_smallest_p <- function(reps, observations = Inf, shared = FALSE) {
  reps <- min(reps, .nsca_orderings(observations))
  if (shared) 1 / (reps + 1) else 1 / reps
}

# Why a test at this replication count and sample size cannot reach `alpha`,
# or NULL when it can. The two causes have different remedies -- more
# permutations, or more observations -- so they are told apart.
.nsca_test_floor_reason <- function(reps, alpha, observations = Inf,
                                    shared = FALSE, label = "test.rep",
                                    threshold = "alpha") {
  if (!isTRUE(reps > 0)) {
    return(NULL)
  }
  free <- .nsca_smallest_p(reps, Inf, shared)
  if (free > alpha) {
    return(sprintf(
      paste0(
        "%s = %d is too few for %s = %s: the permutation test cannot return ",
        "a p-value below 1/%d = %.4f"
      ),
      label, as.integer(reps), threshold, format(alpha),
      as.integer(round(1 / free)), free
    ))
  }
  p_min <- .nsca_smallest_p(reps, observations, shared)
  if (p_min > alpha) {
    return(sprintf(
      paste0(
        "%d observations are too few for %s = %s: the outcome has only %d ",
        "orderings, and the permutation test cannot return a p-value below ",
        "1/%d = %.4f"
      ),
      as.integer(observations), threshold, format(alpha),
      as.integer(.nsca_orderings(observations)), as.integer(round(1 / p_min)),
      p_min
    ))
  }
  NULL
}

# The smallest replication count whose test can reach `alpha`.
.nsca_test_rep_needed <- function(alpha, shared = FALSE) {
  needed <- as.integer(ceiling(round(1 / alpha, 6)))
  if (shared) needed - 1L else needed
}

# Warn, at fitting time, when no component can be significant.
#
# The p-values in the tables come from the shared sequence when there is one
# and from the engines otherwise, so that is the count checked. SCAtools
# raises its own warning for the sufficiency side; nsca_analysis() muffles it
# in favor of this one, which speaks of both components.
.nsca_warn_test_floor <- function(test_rep, shared_rep, alpha,
                                  observations, shared_observations) {
  shared <- isTRUE(shared_rep > 0)
  reps <- if (shared) shared_rep else test_rep
  label <- if (shared) "shared.test.rep" else "test.rep"
  reason <- .nsca_test_floor_reason(
    reps, alpha,
    if (shared) shared_observations else observations,
    shared, label, "test.p_threshold"
  )
  if (is.null(reason)) {
    return(invisible(NULL))
  }
  remedy <- if (.nsca_smallest_p(reps, Inf, shared) > alpha) {
    sprintf(" %s >= %d is required; 1000 or more is common in practice.", label,
            .nsca_test_rep_needed(alpha, shared))
  } else {
    ""
  }
  warning(
    paste0(
      reason, ", so neither component can be significant and joint support ",
      "cannot be reached, regardless of the data.", remedy
    ),
    call. = FALSE
  )
  invisible(NULL)
}

# The same arithmetic for a power simulation, where its consequence is a power
# of exactly 0 that reads as a property of the design.
.nsca_warn_power_floor <- function(test_rep, alpha, n, shared = FALSE) {
  reason <- .nsca_test_floor_reason(test_rep, alpha, min(n), shared,
                                    "test.rep", "alpha")
  if (is.null(reason)) {
    return(invisible(NULL))
  }
  consequence <- if (.nsca_smallest_p(test_rep, Inf, shared) > alpha) {
    sprintf(
      paste0(
        ", so no simulated replication can reject and the estimated power is 0 ",
        "at every sample size. This reflects the number of permutations rather ",
        "than the design. test.rep >= %d is required, and 100 or more is ",
        "advisable for a table that will be reported."
      ),
      .nsca_test_rep_needed(alpha, shared)
    )
  } else {
    short <- n[vapply(n, function(size) {
      .nsca_smallest_p(test_rep, size, shared) > alpha
    }, logical(1L))]
    sprintf(
      paste0(
        ", so the estimated power at n = %s is 0. This reflects the number of ",
        "permutations rather than the design."
      ),
      paste(short, collapse = ", ")
    )
  }
  warning(paste0(reason, consequence), call. = FALSE)
  invisible(NULL)
}

# Warn when the adjustment cannot reject a condition, whatever the data says.
#
# The permutation p-value has a floor (.nsca_smallest_p()), and adjustment
# raises it. By how much depends on the method and on how many conditions are
# strong at once. A condition that is the only strong one among k has its
# p-value multiplied by k under every method stats::p.adjust() offers -- by k
# times the k-th harmonic number under BY -- while under BH, hochberg and
# hommel several equally strong conditions can share a smaller value. Both
# cases are computed with stats::p.adjust() itself, so the warning cannot drift
# from the adjustment it describes. Each component family is adjusted this
# way, and the joint test needs both, so the floor of one is the floor of the
# joint decision.
#
# The lone strong condition is the case the warning is about: it is the common
# one, and a missing rejection there is the easiest to misread as a null
# finding. When even the most favorable case cannot reach alpha, the warning
# says that instead.
#
# A warning rather than a refusal: the adjusted p-values are still correct,
# and a caller may want them for their own reasons. What they must not do is
# read the absence of a rejection as evidence.
.nsca_adjust_floor <- function(k, p_min, method) {
  c(
    lone = stats::p.adjust(c(p_min, rep(1, k - 1L)), method = method)[[1L]],
    best = stats::p.adjust(rep(p_min, k), method = method)[[1L]]
  )
}

.nsca_warn_adjust_floor <- function(k, reps, alpha, method,
                                    observations = Inf, shared = FALSE) {
  if (!isTRUE(is.finite(reps) && reps > 0) || k < 2L) {
    return(invisible(NULL))
  }
  label <- if (shared) "shared.test.rep" else "test.rep"
  p_min <- .nsca_smallest_p(reps, observations, shared)
  floor <- .nsca_adjust_floor(k, p_min, method)
  if (floor[["lone"]] <= alpha) {
    return(invisible(NULL))
  }
  # The count at which a lone strong condition can reach alpha: the lone floor
  # is a fixed multiple of the unadjusted one.
  multiple <- .nsca_adjust_floor(k, 1e-8, method)[["lone"]] / 1e-8
  needed <- as.integer(ceiling(round(multiple / alpha, 6)))
  if (shared) {
    needed <- needed - 1L
  }
  outcome <- if (floor[["best"]] > alpha) {
    sprintf(
      paste0(
        "after '%s' adjustment it is at least %.4f, so no condition can reach ",
        "alpha = %s regardless of the strength of the evidence"
      ),
      method, floor[["best"]], format(alpha)
    )
  } else {
    sprintf(
      paste0(
        "after '%s' adjustment a condition that is the only strong one cannot ",
        "fall below %.4f, so it cannot reach alpha = %s regardless of the ",
        "strength of its evidence"
      ),
      method, floor[["lone"]], format(alpha)
    )
  }
  remedy <- if (.nsca_orderings(observations) < needed) {
    sprintf(
      paste0(
        "A sample of %d has only %d orderings of the outcome, so more ",
        "permutations cannot lower it."
      ),
      as.integer(observations), as.integer(.nsca_orderings(observations))
    )
  } else {
    sprintf("%s >= %d would be required.", label, needed)
  }
  warning(
    sprintf(
      paste0(
        "With %d conditions and %s = %d, the smallest p-value the permutation ",
        "test can return is %.4f, and %s. An absence of joint support here ",
        "would therefore reflect the resolution of the test rather than the ",
        "data. %s"
      ),
      k, label, as.integer(reps), p_min, outcome, remedy
    ),
    call. = FALSE
  )
  invisible(NULL)
}
