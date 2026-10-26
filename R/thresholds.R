# Dual thresholds.
#
# For a given outcome level, necessity and sufficiency answer two different
# questions, and together they cut the condition axis into three regions:
#
# For a high-X statement, low values are out of reach and high values are
# guaranteed. For a low-X statement the ordering reverses. In both cases the
# direction-oriented gap is positive when the two frontiers leave an admitted
# interval between them.
#
# The middle region is the admissible region at that outcome level: the
# interval the frontiers still admit. Its width is the necessity-sufficiency
# interval, the row-level counterpart of admissible_region_share, and it
# collapses to zero exactly where the relation is deterministic at that level.
# "Admits" is the right word: the interval says observations are not ruled out
# there, not that any were seen there.

# Read a threshold column from an engine bottleneck table in actual units.
# The exact values are carried on an attribute; the printed text is the
# fallback and is rounded to three decimals.
.nsca_actual_column <- function(table, condition, rows) {
  column <- table[[condition]]
  actual <- attr(column, "mpx.actual")
  if (is.null(actual)) {
    actual <- attr(table[condition], "mpx.actual")
  }
  if (!is.null(actual)) {
    actual <- as.numeric(actual)
    if (length(actual) == rows) {
      return(list(values = actual, exact = TRUE))
    }
  }
  text <- trimws(as.character(column))
  text[is.na(text)] <- "NA"
  parsed <- suppressWarnings(as.numeric(text))
  # "NN" is an infinite requirement in the engine's flipped coordinates and
  # "NA" a requirement outside the scope; both are non-finite here, and the
  # side-specific readings are attached by the caller.
  parsed[!is.finite(parsed)] <- ifelse(text[!is.finite(parsed)] == "NN", Inf, NA_real_)
  list(values = parsed, exact = FALSE)
}

# Necessity markers keep their ordinary necessity reading, because the
# necessity side is being read as necessity. Only a necessity computation read
# as sufficiency needs the meanings inverted, which is what SCAtools handles.
.nsca_necessity_status <- function(values, x_letter) {
  ifelse(
    is.na(values), "unattainable",
    ifelse(
      is.infinite(values),
      if (identical(x_letter, "H")) "no_minimum" else "no_maximum",
      "estimable"
    )
  )
}

#' Dual thresholds and the necessity-sufficiency interval
#'
#' For each outcome level, reports the necessity and sufficiency thresholds
#' and the direction-oriented distance between them.
#'
#' For a fixed outcome target, the *necessity-sufficiency interval* is the
#' distance on the condition axis between the two thresholds: if a target
#' requires at least 40 units of the condition and 70 units are sufficient for
#' it, the interval extends from 40 to 70. It is not a confidence interval and
#' does not express statistical uncertainty.
#'
#' `necessity_sufficiency_interval` and `overlap_width` are reported on
#' `scale`, the same scale as the two threshold columns whose distance they
#' measure, so that a row can be checked by subtraction.
#' `necessity_sufficiency_interval_actual` and `overlap_width_actual` report
#' the same two quantities in the units of the data, whatever `scale` is
#' requested, and `threshold_gap_actual` is the signed gap in actual units from
#' which both are derived.
#'
#' Under `scale = "percentile"` the reported column is a difference of ranks
#' rather than a distance, so it is interpreted as the share of cases lying
#' between the two thresholds. This is a useful quantity, but it is not a
#' width: equal percentile gaps in a dense and in a sparse part of the
#' distribution correspond to very different distances on the condition axis.
#' For a normalized width that is comparable across studies,
#' `"percentage.range"` or `"sd"` may be more suitable.
#'
#' The interval is also distinct from the component effect sizes. Effect sizes
#' aggregate exclusion over the whole scope, whereas the interval is evaluated
#' at one outcome target. Large component effects therefore do not guarantee a
#' narrow interval at any particular level, and a narrow interval at one target
#' does not imply anything about the others. Reporting both quantities, at
#' several targets, is not redundant.
#'
#' @section The three regions:
#' Necessity and sufficiency address different questions about the same
#' outcome level, and together they divide the condition axis into three parts.
#' For a high-X direction:
#'
#' \describe{
#'   \item{below `necessity_threshold`}{The outcome level is out of reach.}
#'   \item{between the two thresholds}{It is admitted but not guaranteed. This
#'   is the admissible interval, and `necessity_sufficiency_interval` is its
#'   width. It is the horizontal cross-section, at one outcome level, of the
#'   region whose area `admissible_region_share` reports.}
#'   \item{above `sufficiency_threshold`}{The fitted frontier guarantees it.}
#' }
#'
#' For a low-X direction these inequalities reverse: values above the necessity
#' threshold are out of reach and values below the sufficiency threshold are
#' guaranteed. `threshold_gap_actual` is direction-oriented, so a positive
#' value indicates an admitted interval in all four directions. A negative
#' value is reported in `threshold_status`, with its magnitude in
#' `overlap_width` and `overlap_width_actual`, rather than being truncated to
#' zero.
#'
#' The interpretation of a crossing depends on the technique. An exact envelope
#' cannot cross itself: CE-FDH's two frontiers are a running minimum and a
#' running maximum over the same observations, so the floor cannot rise above
#' the ceiling, and a crossing there reflects either the one-tread artifact of
#' a step function or a genuine specification problem. Fitted lines, by
#' contrast, can cross. The line techniques (CR-FDH, CR-VRS, `qr`, `cols`)
#' regress a straight line through the frontier peers, and the narrower the
#' band between the two frontiers, the more easily two separately fitted lines
#' cross inside it. A small crossing under a line technique may therefore
#' indicate a band narrower than the fitting error, which is the pattern an
#' almost-perfect necessary-and-sufficient relation would produce. Such a
#' crossing is consistent with strong correspondence rather than with a
#' misspecified model.
#'
#' The two interpretations are distinguished by magnitude rather than by sign.
#' Within `correspondence.tolerance` the two thresholds are treated as lying at
#' the same X, whichever side the estimate fell on, and `corresponds` is
#' `TRUE`. Beyond the tolerance in the negative direction, the floor lies above
#' the ceiling over an interval, which is incompatible with the relation;
#' `corresponds` is then `NA`, because the correspondence question cannot be
#' answered there.
#'
#' The interval collapses to zero exactly where the relation is deterministic
#' at that level, reported as `"exact correspondence"`: one X value is then both
#' the minimum required and the sufficient level. An estimated interval of zero
#' is consistent with exact correspondence but does not by itself establish
#' it, and a small non-zero interval is best described as near correspondence
#' only when a substantively justified tolerance has been defined. A positive
#' interval does not negate joint support; it indicates imperfect threshold
#' correspondence. Neither NCA nor SCA alone can produce this table, since each
#' supplies one edge of the interval.
#'
#' The two out-of-range markers are interpreted differently on the two sides,
#' because their meaning inverts under contraposition. On the necessity side
#' they retain their necessity interpretation, giving `"no_minimum"`,
#' `"no_maximum"` and `"unattainable"`; on the sufficiency side SCAtools
#' supplies the sufficiency interpretation.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param ceiling A single frontier technique, applied to both sides.
#' @param x Optional condition names or positions.
#' @param scale,outcome_scale Reporting scales, defaulting to those requested
#'   in [nsca_analysis()]. See [SCAtools::sca_scales()].
#' @param convention `"absolute"` or `"directional"`.
#' @param correspondence.tolerance The largest necessity-sufficiency interval,
#'   in the units of the condition, that still counts as threshold
#'   correspondence. `NULL` uses the tolerance with which the model was fitted,
#'   which is `0` unless one was declared. See [nsca_analysis()] for the
#'   reasons why no substantive default is provided.
#'
#'   A non-zero value adds two columns, `corresponds` and
#'   `correspondence_tolerance`, and enables the `"near correspondence"`
#'   category of `threshold_status`. At zero these columns are omitted, because
#'   they would carry no judgment, and a binary verdict placed beside a
#'   continuous interval tends to attract more attention than the interval
#'   itself. The default output therefore reports the interval; a verdict is
#'   added only when a width has been chosen and can be justified.
#' @param inequality Strict (`"strict"`, the default) or inclusive
#'   (`"inclusive"`) inequalities in the reported rules.
#' @param digits Significant digits in the formatted text.
#' @param legacy Whether to append the column names retired in 0.3.0, 0.4.0
#'   and 0.5.x as duplicates. See [nsca_legacy_names()].
#' @return A data frame with one row per condition and outcome level.
#' @seealso [nsca_table()], [nsca_analysis()], [nsca_legacy_names()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_thresholds(fit)
#' nsca_thresholds(fit, scale = "actual")
nsca_thresholds <- function(
    model,
    ceiling = NULL,
    x = NULL,
    scale = NULL,
    outcome_scale = NULL,
    convention = NULL,
    correspondence.tolerance = NULL,
    inequality = c("strict", "inclusive"),
    digits = 6L,
    legacy = FALSE) {
  metadata <- .nsca_metadata(model)
  # NULL means the tolerance the model was fitted with, so a reader who did
  # not set one gets zero and a reader who did does not have to restate it.
  if (is.null(correspondence.tolerance)) {
    correspondence.tolerance <- metadata$correspondence.tolerance
  }
  # Resolved against the conditions this table will cover, by the same rule
  # the fit uses: a width in the units of X may not be stated once for axes
  # that need not share a scale.
  correspondence.tolerance <- .nsca_validate_correspondence(
    correspondence.tolerance, metadata$conditions
  )
  inequality <- match.arg(inequality)
  if (is.null(ceiling)) {
    ceiling <- metadata$ceilings[[1L]]
  }
  if (length(ceiling) != 1L) {
    stop("A threshold table requires exactly one frontier technique.",
         call. = FALSE)
  }
  if (!(ceiling %in% metadata$ceilings)) {
    stop(
      sprintf(
        "Frontier '%s' was not estimated. Available: %s.",
        ceiling, paste(metadata$ceilings, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # The sufficiency side is produced by SCAtools, which already applies the
  # direction-aware scale conversion and the sufficiency reading of the
  # out-of-range markers.
  sufficiency <- SCAtools::sca_thresholds(
    model$sufficiency,
    ceiling = ceiling,
    x = x,
    scale = scale,
    outcome_scale = outcome_scale,
    convention = convention,
    inequality = inequality,
    digits = digits
  )
  condition_names <- unique(sufficiency$condition)

  nec_table <- model$necessity$bottlenecks[[ceiling]]
  if (is.null(nec_table)) {
    stop(
      sprintf("The necessity side produced no threshold table for '%s'.", ceiling),
      call. = FALSE
    )
  }
  nec_levels <- as.numeric(nec_table[[1L]])

  rows <- list()
  # Conditions whose thresholds came from the engine's printed text rather
  # than from its exact attribute. .nsca_actual_column() has always reported
  # this; nothing read it until 0.5.0, so NSCA fell back to three-decimal
  # values in silence while SCAtools warned in the same situation. It matters
  # here more than there: "exact correspondence" is decided against a tolerance
  # of about 1.5e-8, and a gap rounded to three decimals cannot meet it, so a
  # genuinely zero interval would be reported as an admissible one.
  degraded <- character(0L)
  for (i in seq_along(condition_names)) {
    condition <- condition_names[[i]]
    suf_rows <- sufficiency[sufficiency$condition == condition, , drop = FALSE]

    nec_column <- .nsca_actual_column(nec_table, condition, length(nec_levels))
    if (!isTRUE(nec_column$exact)) {
      degraded <- c(degraded, condition)
    }

    # Align the necessity rows to the sufficiency rows by actual outcome level.
    # Both sides were given the same levels, but the engine reverses them for a
    # lower corner (p_bottleneck_table.R: `if (flip.y) mpy <- rev(values)`) and
    # an NSCA model always has one corner of each, so this match does real work
    # on every call rather than being a defensive nicety.
    #
    # The values themselves are bit-identical -- the engine only sorts and
    # filters the steps it is given -- so exact matching is right and a
    # tolerance would be worse: two outcome levels close together could be
    # mis-paired. What was missing is the check that it worked. Before 0.5.0 a
    # failed match left NA in `position`, which made every threshold NA and
    # every row "not estimable", indistinguishable from a genuinely
    # unattainable target.
    position <- match(
      round(suf_rows$outcome_level_actual, 10L), round(nec_levels, 10L)
    )
    if (anyNA(position)) {
      stop(
        sprintf(
          paste0(
            "The necessity and sufficiency threshold tables could not be ",
            "aligned for '%s': %d of %d outcome levels have no counterpart. ",
            "Both sides are given the same levels, so this indicates that the ",
            "engine altered them. Without this alignment, the two thresholds ",
            "in a row would not refer to the same outcome level."
          ),
          condition, sum(is.na(position)), length(position)
        ),
        call. = FALSE
      )
    }
    nec_actual <- nec_column$values[position]
    direction <- unique(suf_rows$direction)
    x_letter <- substr(direction[[1L]], 1L, 1L)
    nec_status <- .nsca_necessity_status(nec_actual, x_letter)
    nec_actual[!is.finite(nec_actual)] <- NA_real_

    suf_actual <- suf_rows$sufficiency_threshold_actual
    gap <- if (identical(x_letter, "H")) {
      suf_actual - nec_actual
    } else {
      nec_actual - suf_actual
    }
    x_span <- diff(metadata$raw$conditions[[condition]]$x_bounds)
    tolerance <- .Machine$double.eps^0.5 * max(1, abs(x_span))
    # Two different tolerances, deliberately not merged. `tolerance` above
    # decides when two computed numbers are the same number -- a floating-point
    # question with a floating-point answer. `substantive` decides how wide an
    # interval a researcher is still willing to call correspondence, which is a
    # claim about the subject matter and has no numerical answer at all. It is
    # zero unless the caller set it, because a package-supplied default would
    # be exactly the unjustified tolerance the framework warns against: the
    # width at which a relation stops being necessary-and-sufficient is not
    # something a package can know.
    substantive <- max(correspondence.tolerance[[condition]], tolerance)
    width_actual <- ifelse(is.finite(gap), pmax(gap, 0), NA_real_)
    overlap_actual <- ifelse(is.finite(gap), pmax(-gap, 0), NA_real_)

    threshold_status <- ifelse(
      nec_status == "estimable" & suf_rows$status == "estimable",
      ifelse(
        is.finite(gap) & abs(gap) <= tolerance,
        "exact correspondence",
        ifelse(
          is.finite(gap) & abs(gap) <= substantive,
          "near correspondence",
          ifelse(
            is.finite(gap) & gap < 0,
            "frontiers overlap", "admissible interval"
          )
        )
      ),
      ifelse(nec_status == "unattainable", "out of reach",
        ifelse(suf_rows$status == "no_threshold", "never guaranteed",
          ifelse(nec_status == "no_minimum", "no minimum required",
            ifelse(nec_status == "no_maximum", "no maximum required",
                   "not estimable")
          )
        )
      )
    )

    rows[[i]] <- data.frame(
      condition = condition,
      outcome = suf_rows$outcome,
      direction = suf_rows$direction,
      relationship = .nsca_relationship_text(
        unname(metadata$direction[[condition]]), condition, metadata$outcome
      ),
      frontier = ceiling,
      convention = suf_rows$convention,
      inequality = inequality,
      condition_scale = suf_rows$condition_scale,
      outcome_scale = suf_rows$outcome_scale,
      outcome_level = suf_rows$outcome_level,
      outcome_level_actual = suf_rows$outcome_level_actual,
      necessity_threshold_actual = nec_actual,
      sufficiency_threshold_actual = suf_actual,
      necessity_status = nec_status,
      sufficiency_status = suf_rows$status,
      threshold_gap_actual = gap,
      necessity_sufficiency_interval_actual = width_actual,
      overlap_width_actual = overlap_actual,
      thresholds_compatible = ifelse(
        is.finite(gap), gap >= -tolerance, NA
      ),
      # Whether the two thresholds correspond at this level, at the tolerance
      # that was declared. With the default tolerance of zero this is TRUE
      # only where the interval is exactly zero, which is the same set of rows
      # threshold_status calls "exact correspondence".
      correspondence_tolerance = substantive,
      # Symmetric in the signed gap, which is the point.
      #
      # A crossing is not one thing. An exact envelope cannot cross itself:
      # CE-FDH's two frontiers are a running minimum and a running maximum
      # over the same points, so floor <= ceiling always, and a crossing there
      # is the one-tread artifact or a genuine specification problem. A fitted
      # line can. CR-FDH, CR-VRS, QR and COLS regress a straight line through
      # the peers, and the narrower the band the more easily two separately
      # fitted lines cross inside it -- so a small crossing under a line
      # technique is the signature of a band thinner than the fitting error,
      # which is what an almost-perfect necessary-and-sufficient relation
      # looks like. Calling that a modeling problem would throw away exactly
      # the case where correspondence is strongest.
      #
      # What separates the two is magnitude, not sign. Within the declared
      # tolerance the two thresholds are at the same X whichever way the
      # estimate fell; beyond it in the negative direction the floor really is
      # above the ceiling over a stretch, and the correspondence question has
      # no answer.
      corresponds = ifelse(
        !is.finite(gap), NA,
        ifelse(abs(gap) <= substantive, TRUE,
               ifelse(gap > 0, FALSE, NA))
      ),
      threshold_status = threshold_status,
      sufficiency_rule = suf_rows$rule,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out <- do.call(rbind, rows)

  # The correspondence columns appear only when a tolerance was declared.
  #
  # At the default of zero they carry no judgment: `corresponds` marks
  # exactly the rows `threshold_status` already calls "exact correspondence",
  # and `correspondence_tolerance` is a column of zeros. What they do carry is
  # a binary verdict sitting beside a continuous interval, and the binary one
  # is the column that gets read. The default output should report the
  # interval and let the reader judge it; a verdict belongs in the table only
  # once someone has chosen a width and can say why that width. Declaring a
  # tolerance is that choice, so it is what brings the columns.
  # Decided from what was declared, not from what was stored. The stored
  # value is max(declared, floating-point tolerance), so it is never exactly
  # zero and testing it would keep the columns on every call.
  if (all(correspondence.tolerance <= 0)) {
    out <- out[, setdiff(names(out),
                         c("corresponds", "correspondence_tolerance")),
               drop = FALSE]
  }

  # Express both thresholds on the requested scale with the very same converter
  # SCAtools used, so the two columns are directly comparable. Anything else
  # risks the two sides differing by a convention rather than by the evidence.
  out$necessity_threshold <- .nsca_rescale(
    metadata, out, "necessity_threshold_actual", scale, convention, inequality
  )
  out$sufficiency_threshold <- .nsca_rescale(
    metadata, out, "sufficiency_threshold_actual", scale, convention, inequality
  )

  # The interval and the overlap are distances between the two threshold
  # columns, so they are reported on the scale those columns are reported on.
  # A width in actual units printed beside two percentages cannot be checked by
  # subtraction, and the Rd told the reader to normalize without offering the
  # normalized number.
  #
  # The magnitude is taken from the reported columns and the sign is carried
  # over from threshold_gap_actual. Every scale is a monotone transformation of
  # X, so the magnitude survives it; but under the "directional" convention a
  # low-X direction is mirrored, and reading the sign off the reported
  # difference would then invert admissible interval and frontiers overlap.
  # The sign is geometry, not presentation, and is settled in actual units.
  reported_gap <- abs(out$sufficiency_threshold - out$necessity_threshold) *
    sign(out$threshold_gap_actual)
  out$necessity_sufficiency_interval <- ifelse(
    is.finite(reported_gap), pmax(reported_gap, 0), NA_real_
  )
  out$overlap_width <- ifelse(
    is.finite(reported_gap), pmax(-reported_gap, 0), NA_real_
  )

  # A reported width of zero does not always mean the two thresholds coincide.
  #
  # The comment above says every scale is a monotone transformation of X so the
  # magnitude survives it. That is true of the affine scales and false of
  # `percentile`, which counts observations: it is monotone but not
  # distance-preserving, so two distinct X values with no observation between
  # them land on the same percentile and the reported width collapses to zero
  # while the actual gap is positive. A row could then print
  # threshold_status = "admissible interval" beside an interval of 0, which
  # contradicts the documented reading that a width of zero is exact
  # correspondence.
  #
  # The status stays as it is: it is decided in actual units, and that is the
  # geometry. What was missing is a way to tell a real zero from a zero the
  # scale produced, so the two are now distinguishable instead of identical.
  if (length(degraded) > 0L) {
    warning(
      sprintf(
        paste0(
          "The frontier engine did not preserve exact necessity thresholds ",
          "for %s; NSCA used its printed output instead, which is rounded to ",
          "three decimals. 'exact correspondence' is assessed against a ",
          "tolerance near 1e-8, so at that precision an interval of zero may ",
          "be reported as an admissible one. Rescaling the affected variable ",
          "may help if that precision is insufficient."
        ),
        paste(unique(degraded), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  collapse_tolerance <- .Machine$double.eps^0.5
  out$interval_collapsed <- is.finite(out$threshold_gap_actual) &
    is.finite(reported_gap) &
    abs(out$threshold_gap_actual) > collapse_tolerance &
    abs(reported_gap) <= collapse_tolerance

  if (any(out$interval_collapsed, na.rm = TRUE)) {
    warning(
      sprintf(
        paste0(
          "%d row(s) report a necessity-sufficiency interval of zero on the ",
          "'%s' scale while the two thresholds differ in actual units. A ",
          "percentile is a step function, so thresholds with no observation ",
          "between them share a rank. For those rows, threshold_gap_actual ",
          "gives the gap in actual units; interval_collapsed identifies them."
        ),
        sum(out$interval_collapsed, na.rm = TRUE),
        if (is.null(scale)) metadata$threshold.x else scale
      ),
      call. = FALSE
    )
  }

  order_columns <- c(
    "condition", "outcome", "direction", "relationship",
    "frontier",
    "convention", "inequality",
    "condition_scale", "outcome_scale", "outcome_level", "outcome_level_actual",
    "necessity_threshold", "sufficiency_threshold",
    "necessity_threshold_actual", "sufficiency_threshold_actual",
    "necessity_status", "sufficiency_status", "threshold_gap_actual",
    "necessity_sufficiency_interval", "overlap_width", "interval_collapsed",
    "necessity_sufficiency_interval_actual", "overlap_width_actual",
    "thresholds_compatible",
    "correspondence_tolerance", "corresponds",
    "threshold_status", "sufficiency_rule"
  )
  # intersect(), not the bare list: the two correspondence columns are removed
  # above when no tolerance was declared, and a whitelist that named a column
  # the frame no longer has would fail rather than order what is there.
  out <- out[, intersect(order_columns, names(out)), drop = FALSE]
  rownames(out) <- NULL
  .nsca_append_legacy(out, legacy)
}

# Accept the retired 'boundary' argument for one more release cycle. The name
# was retired because the framework uses "boundary" for the theoretical line
# that separates an expected empty space from the compatible region -- the
# thing a frontier estimates -- not for the strictness of a reported
# inequality.
# Convert an actual-unit column to the reporting scale, one condition at a
# time, because the percentile and standard-deviation scales depend on that
# condition's own data.
.nsca_rescale <- function(metadata, out, column, scale, convention, inequality) {
  scale <- if (is.null(scale)) metadata$threshold.x else scale
  convention <- if (is.null(convention)) metadata$convention else convention
  values <- rep(NA_real_, nrow(out))
  for (condition in unique(out$condition)) {
    rows <- out$condition == condition
    raw <- metadata$raw$conditions[[condition]]
    letter <- substr(unname(metadata$direction[[condition]]), 1L, 1L)
    values[rows] <- SCAtools::sca_rescale(
      out[[column]][rows],
      reference = raw$x,
      bounds = raw$x_bounds,
      scale = scale,
      letter = letter,
      convention = convention,
      inequality = inequality
    )
  }
  values
}
