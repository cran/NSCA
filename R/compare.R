# Average effect beside a condition claim, on both sides at once.
#
# SCAtools has sca_compare() and sca_locationshift(); NSCA had neither. A
# reader who wanted the regression comparison for a necessary-and-sufficient
# model therefore had to reach into model$sufficiency and run the SCAtools
# function on it, which answers the question for one of the two components and
# returns a table with nothing on it saying so.
#
# That is the whole reason these are not one-line delegations. The average
# effect is a single line through the cloud and there is one of it; the
# boundary it is being compared against is one of two, and which one changes
# the reading completely. A flat slope with an empty upper-left corner and a
# flat slope with an empty lower-right corner are different findings, and a
# table that reports "configuration II" without naming the corner has thrown
# away the half of the answer this package exists to supply.
#
# So both functions return one row per condition *and component*, and the
# comparison columns repeat down the pair: same beta, same slope, two
# boundaries, two readings.

# The four configurations of the condition analysis in degree framework, with
# the component's own name in the label.
#
# The numbering is the framework's, where it is stated for sufficiency. The
# necessity rows ask the same two questions of the other corner, which is the
# older argument rather than a new typology -- NCA's case for itself has always
# been that a condition can be necessary without raising the expected outcome,
# which is configuration II read on the necessity side. The numerals are
# borrowed all the same, and the documentation says so.
.nsca_configuration_label <- list(
  necessity = c(
    I = "neither",
    II = "necessity-dominant",
    III = "mean-dominant",
    IV = "mean and necessity"
  ),
  sufficiency = c(
    I = "neither",
    II = "sufficiency-dominant",
    III = "mean-dominant",
    IV = "mean and sufficiency"
  )
)

.nsca_configuration_reading <- list(
  necessity = c(
    I = paste0(
      "Neither the average effect nor the empty corner reaches its stated ",
      "cut. There is no evidence that %s is necessary for %s, nor of an ",
      "average effect."
    ),
    II = paste0(
      "The average slope is flat, but the corner excluded by necessity is ",
      "largely empty: below some level of %s, cases with high %s no longer ",
      "occur, although higher levels of the condition do not raise the ",
      "expected outcome. This is the situation for which NCA was proposed, ",
      "and a regression-only design would report it as no effect."
    ),
    III = paste0(
      "Higher %s is associated with higher expected %s on average, but cases ",
      "with a low level of the condition and a high outcome remain. The ",
      "average effect is supported; necessity is not."
    ),
    IV = paste0(
      "Higher %s is associated with higher expected %s on average, and the ",
      "corner excluded by necessity is largely empty. These are two distinct ",
      "findings."
    )
  ),
  sufficiency = c(
    I = paste0(
      "Neither the average effect nor the empty corner reaches its stated ",
      "cut. There is no evidence that %s is sufficient for %s, nor of an ",
      "average effect."
    ),
    II = paste0(
      "The average slope is flat, but the corner excluded by sufficiency is ",
      "largely empty: beyond some level of %s, cases with low %s no longer ",
      "occur. A regression-only design would report this relationship as ",
      "absent. nsca_locationshift() can be used to examine whether the ",
      "frontier shifts the whole distribution or marks a genuine threshold."
    ),
    III = paste0(
      "Higher %s is associated with higher expected %s on average, but ",
      "counterexamples remain. The average effect is supported; sufficiency ",
      "is not."
    ),
    IV = paste0(
      "Higher %s is associated with higher expected %s on average, and the ",
      "corner excluded by sufficiency is largely empty. These are two ",
      "distinct findings."
    )
  )
)

# Cohen (1988): .10 small, .30 medium, .50 large, for a standardized
# coefficient. Dul (2016): <.10 small, .10-.30 medium, for an effect size.
# Two numbers because two traditions, each the bottom of its own "medium".
.nsca_beta_cut_default <- 0.30
.nsca_d_cut_default <- 0.10

.nsca_configuration <- function(average_present, boundary_present, beta, d) {
  if (!is.finite(beta) || !is.finite(d)) {
    return(NA_character_)
  }
  if (!average_present && !boundary_present) {
    "I"
  } else if (!average_present && boundary_present) {
    "II"
  } else if (average_present && !boundary_present) {
    "III"
  } else {
    "IV"
  }
}

# The two per-component cuts, validated the way `relevance` is.
#
# `d_cut` decides a reported classification, exactly as `relevance` decides a
# reported support verdict, so it is worth the same protection: a name that is
# not one of the two components is an error rather than a hint to fall back on
# position, because c(sufficiency = 0.2, nec = 0.1) would otherwise be applied
# in reverse with no warning.
.nsca_d_cuts <- function(d_cut, relevance) {
  if (is.null(d_cut)) {
    stated <- relevance
    if (is.null(stated) || length(stated) != 2L) {
      stated <- stats::setNames(c(NA_real_, NA_real_), .nsca_components)
    }
    values <- vapply(.nsca_components, function(component) {
      value <- suppressWarnings(as.numeric(stated[[component]]))
      if (is.finite(value)) value else .nsca_d_cut_default
    }, numeric(1L))
    return(stats::setNames(unname(values), .nsca_components))
  }
  .nsca_validate_relevance(d_cut)
}

#' Average effect beside both condition claims
#'
#' Reports the ordinary least-squares coefficient beside the necessity and the
#' sufficiency effect size for the same relationship, and classifies each pair
#' into the four configurations of the condition analysis in degree framework.
#'
#' @section Why two rows per condition:
#' There is one average effect and there are two boundaries. A flat slope with
#' an empty necessity corner and a flat slope with an empty sufficiency corner
#' are different findings, so a single row per condition would have to select
#' one of them without indicating which. Each condition therefore produces a
#' `necessity` row and a `sufficiency` row; `beta` and the slope columns repeat
#' down the pair, because they describe the same line.
#'
#' @section What is being compared:
#' `beta` is the standardized regression coefficient, which describes how far
#' the expected outcome moves with the condition. `effect_size` is the share of
#' the scope left empty in the corner that the component's claim excludes.
#' Neither implies the other. A relationship can have a strong average effect
#' and no empty corner, or an empty corner and almost no slope.
#'
#' @section The four configurations:
#' \describe{
#'   \item{I}{Neither the average effect nor the empty corner reaches its cut.}
#'   \item{II}{A flat slope with a largely empty corner, often reflecting a
#'   level shift or a threshold. This configuration is not detectable by a
#'   design that reports only a regression coefficient, because configurations
#'   I and II have the same flat slope.}
#'   \item{III}{A slope with counterexamples still present. The condition
#'   shifts the expected outcome without bounding it.}
#'   \item{IV}{Both an average effect and a largely empty corner.}
#' }
#'
#' The numbering follows the framework, where it is stated for sufficiency.
#' The necessity rows apply the same two questions to the other corner. This
#' reflects the argument on which NCA was proposed, namely that a condition can
#' be necessary without raising the expected outcome, rather than a new
#' typology; the numerals are nevertheless borrowed, as acknowledged here.
#'
#' @section The cuts:
#' `beta_cut` defaults to 0.30, the lower edge of the medium band in Cohen
#' (1988). `d_cut` defaults to the `relevance` given to [nsca_analysis()],
#' per component, when one was pre-specified, because the researcher has then
#' already stated a criterion that this function does not override; otherwise
#' it defaults to 0.10, the lower edge of the medium band in Dul (2016). One
#' number may be supplied for both components, or
#' `c(necessity = , sufficiency = )` for separate values. The cuts used are
#' returned as columns, so that each table records the values that produced
#' its classification.
#'
#' These are conventions, not benchmarks calibrated for this design.
#'
#' `effect_size_null`, `effect_size_excess` and `effect_size_critical` are
#' reported beside the classification for that reason. A condition can reach
#' `d_cut` with an effect size that largely reflects the envelope fitting
#' sample extremes, and the label would read "dominant" in either case. The cut
#' is kept at its conventional value; the table adds the values that qualify
#' it. See [nsca_table()] for how the floor is measured.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x Optional condition names or positions. All conditions by default.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param beta_cut Absolute standardized coefficient at or above which the
#'   average effect counts as present.
#' @param d_cut Effect size at or above which the empty corner counts as
#'   present: one number, or `c(necessity = , sufficiency = )`. Defaults to the
#'   model's `relevance`, then to 0.10.
#' @return A data frame with one row per condition and component.
#' @references
#' Chan, Y. (in press). *Mean Effects Are Not Sufficiency Effects: How
#' Sufficient Condition Analysis Complements Regression Analysis*.
#' University of Leeds.
#'
#' Chan, Y. (in press). *From Necessity to Sufficiency and
#' Necessity-and-Sufficiency: Extending Condition Analysis in Degree*.
#' University of Leeds.
#' @seealso [nsca_locationshift()], [nsca_reference()], [nsca_table()],
#'   [SCAtools::sca_compare()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_compare(fit)
nsca_compare <- function(model, x = NULL, ceiling = NULL,
                         beta_cut = .nsca_beta_cut_default, d_cut = NULL) {
  metadata <- .nsca_metadata(model)
  conditions <- .nsca_select_conditions(model, x)
  ceiling <- .nsca_resolve_ceiling(metadata, ceiling)

  beta_cut <- abs(suppressWarnings(as.numeric(beta_cut))[[1L]])
  if (!is.finite(beta_cut)) {
    stop("'beta_cut' must be a finite number.", call. = FALSE)
  }
  cuts <- .nsca_d_cuts(d_cut, metadata$relevance)

  # One line per condition, not one per row: the two components share it.
  reference <- nsca_reference(model, x = conditions, line = "ols")
  components <- .nsca_component_rows(model, conditions, ceiling)

  rows <- list()
  index <- 0L
  for (condition in conditions) {
    line <- reference[reference$condition == condition, , drop = FALSE][1L, ]
    beta <- line$beta
    for (component in .nsca_components) {
      row <- components[
        components$component == component & components$condition == condition, ,
        drop = FALSE
      ]
      if (nrow(row) == 0L) {
        stop(
          sprintf("No %s result for condition %s with frontier %s.",
                  component, condition, ceiling),
          call. = FALSE
        )
      }
      row <- row[1L, ]
      d <- row$effect_size
      d_cut_here <- unname(cuts[[component]])
      average_present <- isTRUE(is.finite(beta) && abs(beta) >= beta_cut)
      boundary_present <- isTRUE(is.finite(d) && d >= d_cut_here)
      configuration <- .nsca_configuration(
        average_present, boundary_present, beta, d
      )

      index <- index + 1L
      rows[[index]] <- data.frame(
        condition = condition,
        outcome = metadata$outcome,
        direction = unname(metadata$direction[[condition]]),
        component = component,
        empty_corner = row$empty_corner,
        boundary_type = row$boundary_type,
        frontier = ceiling,
        observations = row$observations,
        beta = beta,
        beta_magnitude = SCAtools::sca_magnitude(beta, scale = "beta"),
        slope = line$slope,
        slope_se = line$slope_se,
        slope_p = line$slope_p,
        r_squared = line$r_squared,
        effect_size = d,
        magnitude = SCAtools::sca_magnitude(d),
        # See the note in SCAtools::sca_compare(): the classification
        # thresholds the raw effect size, and the floor says how much of
        # that effect size random pairing alone accounts for.
        effect_size_null = row$effect_size_null,
        effect_size_excess = row$effect_size_excess,
        effect_size_critical = row$effect_size_critical,
        p_value = row$p_value,
        beta_cut = beta_cut,
        d_cut = d_cut_here,
        average_effect_present = average_present,
        boundary_present = boundary_present,
        configuration = configuration,
        configuration_label = if (is.na(configuration)) {
          NA_character_
        } else {
          unname(.nsca_configuration_label[[component]][[configuration]])
        },
        reading = if (is.na(configuration)) {
          NA_character_
        } else {
          sprintf(
            unname(.nsca_configuration_reading[[component]][[configuration]]),
            condition, metadata$outcome
          )
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Location-shift diagnostic for both frontiers
#'
#' Places the ordinary least-squares slope, the median (quantile regression at
#' `tau = 0.5`) slope and the fitted frontier slope side by side, once for each
#' component.
#'
#' When all three agree, the frontier is a *location shift*: the central
#' tendency has shifted as a whole to one edge of the data, so the component
#' adds a level to what the regression already describes but does not add a
#' slope. When the frontier slope departs from the other two, the boundary has
#' a shape of its own and the relationship is not a shift.
#'
#' The two components can give different answers, and may frequently do so. A
#' band whose upper edge tracks the regression while its lower edge does not
#' describes a relationship with a floor of its own shape and a ceiling that is
#' only the mean shifted upward, and this asymmetry is not visible in a table
#' with one row per condition.
#'
#' This is a natural follow-up to configuration II of [nsca_compare()], in
#' which the average slope is flat but an excluded corner is empty.
#'
#' Only frontier techniques that fit a line report a slope. A step frontier
#' (`ce_fdh`, `ce_vrs`) has none by construction, so `frontier_slope` is `NA`
#' for those techniques and the comparison requires a regression-type frontier
#' such as `cr_fdh`.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x Optional condition names or positions. All conditions by default.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param tolerance Relative spread among the three slopes at or below which
#'   they are reported as agreeing.
#' @return A data frame with one row per condition and component.
#' @seealso [nsca_compare()], [SCAtools::sca_locationshift()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "cr_fdh")
#' nsca_locationshift(fit)
nsca_locationshift <- function(model, x = NULL, ceiling = NULL,
                               tolerance = 0.10) {
  metadata <- .nsca_metadata(model)
  conditions <- .nsca_select_conditions(model, x)
  ceiling <- .nsca_resolve_ceiling(metadata, ceiling)
  tolerance <- suppressWarnings(as.numeric(tolerance))[[1L]]
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("'tolerance' must be a non-negative number.", call. = FALSE)
  }

  components <- .nsca_component_rows(model, conditions, ceiling)

  # The two central-tendency slopes are SCAtools' own rather than a second
  # computation here, for the reason nsca_reference() delegates: the
  # sufficiency layer is an sca_result fitted on exactly these observations,
  # so the lines that package fits are the lines this one would. Two
  # implementations of a median slope can drift apart; one cannot. It also
  # keeps quantreg out of this package's dependencies, where it would have
  # been declared for a single call.
  #
  # Only the central tendency comes from there. Both frontier slopes are read
  # from the component tables below, so the necessity and sufficiency
  # boundaries reach this table by one route rather than two.
  central <- SCAtools::sca_locationshift(
    model$sufficiency, x = conditions, ceiling = ceiling,
    tolerance = tolerance
  )

  rows <- list()
  index <- 0L
  for (condition in conditions) {
    line <- central[central$condition == condition, , drop = FALSE][1L, ]
    ols_slope <- line$ols_slope
    median_slope <- line$median_slope

    for (component in .nsca_components) {
      row <- components[
        components$component == component & components$condition == condition, ,
        drop = FALSE
      ]
      frontier_slope <- if (nrow(row) == 0L) NA_real_ else row$slope[[1L]]
      corner <- if (nrow(row) == 0L) NA_real_ else row$empty_corner[[1L]]

      slopes <- c(ols_slope, median_slope, frontier_slope)
      spread <- if (all(is.finite(slopes))) {
        scale <- max(abs(slopes))
        if (scale > 0) (max(slopes) - min(slopes)) / scale else 0
      } else {
        NA_real_
      }

      index <- index + 1L
      rows[[index]] <- data.frame(
        condition = condition,
        outcome = metadata$outcome,
        component = component,
        empty_corner = corner,
        frontier = ceiling,
        observations = line$observations,
        ols_slope = ols_slope,
        median_slope = median_slope,
        frontier_slope = frontier_slope,
        relative_spread = spread,
        tolerance = tolerance,
        location_shift = if (is.finite(spread)) spread <= tolerance else NA,
        reading = if (!is.finite(spread)) {
          paste0(
            "Not assessable: a step frontier has no slope, so the comparison ",
            "requires a regression-type frontier such as cr_fdh."
          )
        } else if (spread <= tolerance) {
          sprintf(
            paste0(
              "The three slopes agree. The %s frontier corresponds to the ",
              "central tendency shifted to the edge of the data, so that ",
              "component adds a level rather than a slope."
            ),
            component
          )
        } else {
          sprintf(
            paste0(
              "The %s frontier slope departs from the central tendency, so ",
              "that boundary has a shape of its own. This is not a location ",
              "shift."
            ),
            component
          )
        },
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Both component tables, restricted to one technique and one set of conditions.
#
# Read from the component tables rather than from the engines, so that a number
# in one of these tables and the same number in nsca_results() cannot drift
# apart: there is one route from an engine parameter to a reported column, and
# it is the one nsca_results() already uses.
.nsca_component_rows <- function(model, conditions, ceiling) {
  out <- rbind(
    .nsca_component_table(model, "necessity"),
    .nsca_component_table(model, "sufficiency")
  )
  out[out$frontier == ceiling & out$condition %in% conditions, , drop = FALSE]
}

.nsca_resolve_ceiling <- function(metadata, ceiling) {
  available <- metadata$ceilings
  if (is.null(ceiling)) {
    if (length(available) == 0L) {
      stop("The model contains no estimated frontier technique.",
           call. = FALSE)
    }
    return(available[[1L]])
  }
  ceiling <- as.character(ceiling)
  if (length(ceiling) != 1L) {
    stop("This function requires exactly one frontier technique.",
         call. = FALSE)
  }
  if (!ceiling %in% available) {
    stop(
      sprintf(
        "The model did not estimate '%s'. It estimated: %s.",
        ceiling, paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  ceiling
}
