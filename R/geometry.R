# Geometry shared by the two frontiers.
#
# The admissible region is computed directly, by integrating the vertical gap
# the two fitted frontiers, rather than algebraically as 1 - d_nec - d_suf + O.
# The two routes are equal in theory: the empty spaces and the admissible region tile
# the scope, so inclusion and exclusion give the identity exactly. They are not
# equal in practice, because the algebraic route uses the effect sizes the
# engines computed from their own internal frontiers while the geometric route
# uses the frontiers this package reconstructs.
#
# That is precisely why the geometric route is the one used. The difference
# between the two is reported as reconstruction_error, and it is the only
# available check that the reconstruction still matches what the engine
# actually fitted. A wrongly reconstructed frontier still looks like a
# frontier on a plot and still produces plausible areas; only the identity
# refuses to close.

.nsca_grid_size <- 2001L

.nsca_numeric_column <- function(column) {
  suppressWarnings(as.numeric(as.character(column)))
}

# The scope the engine will use: the observed range, widened by any theoretical
# scope supplied by the caller.
.nsca_scope_pair <- function(scope, index, offset) {
  if (is.null(scope)) {
    return(c(NA_real_, NA_real_))
  }
  segment <- if (is.list(scope)) {
    if (index <= length(scope)) scope[[index]] else scope[[1L]]
  } else {
    values <- as.numeric(scope)
    if (length(values) == 4L) {
      values
    } else if (length(values) >= index * 4L) {
      values[((index - 1L) * 4L + 1L):(index * 4L)]
    } else {
      return(c(NA_real_, NA_real_))
    }
  }
  segment <- as.numeric(segment)
  if (length(segment) < offset + 2L) {
    return(c(NA_real_, NA_real_))
  }
  segment[c(offset + 1L, offset + 2L)]
}

.nsca_bounds <- function(observed, scope_pair) {
  observed <- observed[is.finite(observed)]
  if (length(observed) == 0L) {
    stop("A variable has no usable observations.", call. = FALSE)
  }
  c(
    min(c(observed, scope_pair), na.rm = TRUE),
    max(c(observed, scope_pair), na.rm = TRUE)
  )
}

# How much of the declared scope the observations actually span, per axis.
#
# Every NSCA quantity is an area measured against the scope, so as this falls
# toward zero the numbers stop describing the data and start describing the
# scope. The limiting case is instructive rather than exotic: a constant
# outcome placed at the center of an imposed unit scope leaves the whole upper
# half empty and the whole lower half empty, so d_nec = d_suf = 0.5, the band
# closes, and every joint index reports perfect joint support for data that
# carries no information whatsoever. Moving the constant to 0.8 changes the
# answer again, which is the tell: nothing about the data moved.
#
# Two different pictures are called a "horizontal line" and they land in
# opposite places. A flat *frontier* over a filled scope means no empty zone at
# all and every index is zero, correctly. A flat *data cloud* means the
# indices are at their maximum, wrongly. Only the second is flagged here.
.nsca_degenerate_tolerance <- 1e-8
.nsca_thin_coverage <- 0.05

.nsca_axis_coverage <- function(values, bounds) {
  span <- suppressWarnings(diff(as.numeric(bounds)))
  values <- values[is.finite(values)]
  if (length(span) != 1L || !is.finite(span) || span <= 0 ||
      length(values) == 0L) {
    return(NA_real_)
  }
  diff(range(values)) / span
}

# Complete pairs only, matching what the engine fits on.
#
# Both axes are bounded by the same rows. Until 0.5.0 the condition axis used
# that condition's complete pairs while the outcome axis used every row with a
# finite outcome, whether or not the condition was observed there. That matters
# more here than in SCAtools: `.nsca_band_for()` divides by
# diff(x_bounds) * diff(y_bounds), so with missing data the admissible region
# was measured against a taller box than the one the engines used for d_nec and
# d_suf, and the area identity could not close. reconstruction_error was
# reporting a fault in this function.
#
# `y_bounds` at the top level stays shared: it lays out the outcome levels, and
# the necessity and sufficiency threshold tables are only readable side by side
# if both were tabulated at the same levels. The per-condition `y_bounds` is
# the denominator.
.nsca_raw_data <- function(data, x_names, y_name, scope,
                           na.action = "pairwise") {
  frame <- as.data.frame(data)
  y_raw <- .nsca_numeric_column(frame[[y_name]])
  y_scope <- .nsca_scope_pair(scope, 1L, 2L)

  x_raw <- lapply(x_names, function(nm) .nsca_numeric_column(frame[[nm]]))
  names(x_raw) <- x_names
  usable <- lapply(x_raw, function(column) is.finite(column) & is.finite(y_raw))
  if (identical(na.action, "listwise")) {
    common <- Reduce(`&`, usable)
    usable <- lapply(usable, function(ignored) common)
  }
  covered <- Reduce(`|`, usable)

  conditions <- vector("list", length(x_names))
  names(conditions) <- x_names
  for (i in seq_along(x_names)) {
    keep <- usable[[i]]
    conditions[[i]] <- list(
      x = x_raw[[i]][keep],
      y = y_raw[keep],
      x_bounds = .nsca_bounds(x_raw[[i]][keep], .nsca_scope_pair(scope, i, 0L)),
      y_bounds = .nsca_bounds(y_raw[keep], y_scope)
    )
  }
  list(
    outcome = y_raw[covered],
    y_bounds = .nsca_bounds(y_raw[covered], y_scope),
    na.action = na.action,
    conditions = conditions
  )
}

# The exact free-disposal-hull frontier for one empty corner, evaluated on a
# grid. Which running extreme is taken depends only on which corner is empty.
.nsca_step_frontier <- function(x, y, corner, grid) {
  vapply(
    grid,
    function(t) {
      subset <- if (corner %in% c(1L, 3L)) y[x <= t] else y[x >= t]
      if (length(subset) == 0L) {
        return(NA_real_)
      }
      if (corner %in% c(1L, 2L)) max(subset) else min(subset)
    },
    numeric(1L)
  )
}

# Piecewise-linear frontier through the peers, for ce_vrs.
.nsca_peer_frontier <- function(peers, grid) {
  if (is.null(peers) || !is.matrix(peers) || nrow(peers) < 1L) {
    return(rep(NA_real_, length(grid)))
  }
  peers <- peers[stats::complete.cases(peers), , drop = FALSE]
  if (nrow(peers) < 1L) {
    return(rep(NA_real_, length(grid)))
  }
  if (nrow(peers) == 1L) {
    return(rep(peers[1L, 2L], length(grid)))
  }
  ordered <- peers[order(peers[, 1L]), , drop = FALSE]
  keep <- !duplicated(ordered[, 1L])
  stats::approx(
    ordered[keep, 1L], ordered[keep, 2L], xout = grid, rule = 2
  )$y
}

#' Frontier values on a grid
#'
#' @param x,y Condition and outcome values for one condition.
#' @param corner The empty corner this frontier bounds.
#' @param ceiling Frontier technique.
#' @param grid Condition values at which to evaluate.
#' @param slope,intercept Line coefficients, for regression-type frontiers.
#' @param peers Peer matrix, used by `ce_vrs`.
#' @return A numeric vector of frontier values, or `NA` where undefined.
#' @noRd
.nsca_frontier <- function(x, y, corner, ceiling, grid,
                           slope = NA_real_, intercept = NA_real_,
                           peers = NULL) {
  if (identical(ceiling, "ce_fdh")) {
    return(.nsca_step_frontier(x, y, corner, grid))
  }
  if (identical(ceiling, "ce_vrs")) {
    return(.nsca_peer_frontier(peers, grid))
  }
  if (is.finite(slope) && is.finite(intercept)) {
    return(intercept + slope * grid)
  }
  rep(NA_real_, length(grid))
}

# A step frontier overlaps itself by construction. Between two consecutive
# observations the upper staircase holds the lower observation's value while
# the lower staircase already holds the upper one, so the two empty zones
# overlap in every tread. The overlap is about one step, of order 1/n, and is
# discretization rather than a modeling failure. Only smooth frontiers can
# cross in a way that says something is wrong.
.nsca_step_ceilings <- "ce_fdh"

#' Undetermined band and overlap between two frontiers
#'
#' @param upper,lower Frontier values on a common grid.
#' @param grid The grid.
#' @param y_bounds Outcome axis bounds, used to clip to the scope.
#' @return A list with the band area, the overlap area, and the maximum depth
#'   to which the frontiers were inverted.
#' @noRd
.nsca_band <- function(upper, lower, grid, y_bounds) {
  usable <- is.finite(upper) & is.finite(lower)
  if (!any(usable)) {
    return(list(
      area = NA_real_, overlap = NA_real_, depth = NA_real_, crossed = NA
    ))
  }
  upper <- pmin(pmax(upper, y_bounds[[1L]]), y_bounds[[2L]])
  lower <- pmin(pmax(lower, y_bounds[[1L]]), y_bounds[[2L]])
  gap <- upper - lower
  gap[!usable] <- 0
  band <- gap
  band[band < 0] <- 0
  overlap <- -gap
  overlap[overlap < 0] <- 0
  list(
    area = .nsca_integrate(band, grid),
    overlap = .nsca_integrate(overlap, grid),
    depth = max(c(0, overlap)),
    crossed = any(gap[usable] < -.Machine$double.eps^0.5)
  )
}

# Trapezoidal integration, written out because utils::trapz does not exist and
# base R gained a trapezoid helper only recently.
.nsca_integrate <- function(values, grid) {
  n <- length(grid)
  if (n < 2L) {
    return(0)
  }
  sum((values[-1L] + values[-n]) / 2 * diff(grid))
}

# Where the engines widened the declared scope, and by how much.
#
# `p_loop_data.R` builds the scope it divides by as the declared scope widened
# to hold every observation, so a declared bound that excludes a data point is
# replaced by the data point. Every d_nec, every d_suf and therefore every
# joint index is an area over that widened scope, while the declared scope is
# what the researcher wrote down and will report.
#
# Computed but not warned about here. `SCAtools::sca_analysis()` runs first in
# every NSCA fit and raises the warning already, naming the same conditions
# and bounds; a second identical warning from this layer would say nothing
# new. What this layer adds is the record: stored in the metadata, printed in
# the settings block, and returned bound by bound by nsca_scope().
.nsca_scope_adjustment <- function(scope, raw, outcome_name) {
  empty <- data.frame(
    condition = character(0L), axis = character(0L),
    bound = character(0L), declared = numeric(0L), used = numeric(0L),
    stringsAsFactors = FALSE
  )
  if (is.null(scope)) {
    return(empty)
  }
  conditions <- names(raw$conditions)
  rows <- list()
  for (i in seq_along(conditions)) {
    condition <- conditions[[i]]
    item <- raw$conditions[[condition]]
    axes <- list(
      list(axis = condition,
           declared = .nsca_scope_pair(scope, i, 0L),
           used = item$x_bounds),
      list(axis = outcome_name,
           declared = .nsca_scope_pair(scope, 1L, 2L),
           used = item$y_bounds)
    )
    for (entry in axes) {
      for (side in 1:2) {
        declared <- entry$declared[[side]]
        used <- entry$used[[side]]
        if (!is.finite(declared) || !is.finite(used) ||
            isTRUE(all.equal(declared, used))) {
          next
        }
        rows[[length(rows) + 1L]] <- data.frame(
          condition = condition,
          axis = entry$axis,
          bound = if (side == 1L) "lower" else "upper",
          declared = declared,
          used = used,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L) empty else do.call(rbind, rows)
}

# The NCA engine writes a warning-suppression flag into the caller's global
# environment and does not remove it. SCAtools clears the one its own engine
# call created, but NSCA calls NCA::nca_analysis() directly afterward, which
# creates it again. Removed here unless it was already there before the fit.
.nsca_clear_engine_scope_flag <- function(existed) {
  if (existed) {
    return(invisible(NULL))
  }
  if (exists("scope_warnings_nca", envir = globalenv(), inherits = FALSE)) {
    rm("scope_warnings_nca", envir = globalenv())
  }
  invisible(NULL)
}

.nsca_scope_area <- function(metadata, condition) {
  item <- metadata$raw$conditions[[condition]]
  if (is.null(item)) {
    return(NA_real_)
  }
  area <- diff(as.numeric(item$x_bounds)) * diff(as.numeric(item$y_bounds))
  if (length(area) != 1L || !is.finite(area)) NA_real_ else area
}

.nsca_scope_is_adjusted <- function(metadata, condition) {
  adjustment <- metadata$scope_adjustment
  if (is.null(adjustment) || nrow(adjustment) == 0L) {
    return(FALSE)
  }
  condition %in% adjustment$condition
}

#' The scope by which each component effect size is divided
#'
#' Both components of an NSCA result are empty areas divided by the area of
#' the X-Y scope, so the scope is the denominator of every quantity the
#' analysis reports, including the joint indices. It is not necessarily the
#' scope that was requested: the NCA engine widens a declared scope until it
#' contains every observation, which changes the denominator without changing
#' the scope the researcher declared.
#'
#' This function reports both. `declared_*` is `NA` when no scope was declared
#' and the observed range was used. `used_*` is the scope by which the engines
#' divided, and `adjusted` is `TRUE` for a condition on which the two differ.
#' This is the counterpart of [SCAtools::sca_scope()].
#'
#' @param model An object returned by [nsca_analysis()].
#' @return A data frame with one row per condition.
#' @seealso [nsca_table()], [nsca_analysis()], [SCAtools::sca_scope()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
#'                      scope = c(0, 1, 0, 1))
#' nsca_scope(fit)
nsca_scope <- function(model) {
  metadata <- .nsca_metadata(model)
  conditions <- metadata$conditions
  declared <- metadata$scope_declared
  rows <- lapply(seq_along(conditions), function(i) {
    condition <- conditions[[i]]
    item <- metadata$raw$conditions[[condition]]
    x_declared <- .nsca_scope_pair(declared, i, 0L)
    y_declared <- .nsca_scope_pair(declared, 1L, 2L)
    data.frame(
      condition = condition,
      outcome = metadata$outcome,
      source = metadata$scope_source,
      declared_xmin = x_declared[[1L]],
      declared_xmax = x_declared[[2L]],
      declared_ymin = y_declared[[1L]],
      declared_ymax = y_declared[[2L]],
      used_xmin = as.numeric(item$x_bounds)[[1L]],
      used_xmax = as.numeric(item$x_bounds)[[2L]],
      used_ymin = as.numeric(item$y_bounds)[[1L]],
      used_ymax = as.numeric(item$y_bounds)[[2L]],
      scope_area = .nsca_scope_area(metadata, condition),
      adjusted = .nsca_scope_is_adjusted(metadata, condition),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

# Refuse a per-condition scope whose outcome half is not the same in every
# segment. The counterpart of SCAtools' .sca_validate_scope_y(), and the reason
# it matters more here.
#
# The two components are not a source of this. Necessity and sufficiency are
# two corners of one scatter plot: `nsca_analysis()` hands both engines the
# same `x`, the same `y` and the same `scope`, and only `corner` differs.
# `corner` drives `flip.y`, which reverses the order of the outcome levels
# (`matrix(rev(values))`) and leaves the box alone. The two components can
# therefore never disagree about the scope.
#
# Conditions are. With `x = c("A", "B")` and a per-condition scope whose Y
# halves differ, `d_nec` and `d_suf` for condition 2 are areas over condition
# 2's box, because the engines read `scope[[id.x]]`, while `.nsca_band_for()`
# clips and divides by `raw$conditions[[i]]$y_bounds`, which this package
# builds from the first segment. The admissible region is then measured against
# a different rectangle from the two components, so
#
#   1 - d_nec - d_suf + overlap_share  =  admissible_region_share
#
# stops holding. `reconstruction_error` catches it and `geometry_acceptable`
# goes FALSE, which is the check doing its job -- but the message a user sees
# is "the area identity does not close", which nothing connects back to a
# scope written two different ways. Refusing the input says it at the point
# where it can still be fixed.
.nsca_scope_segments <- function(scope, n) {
  if (is.null(scope)) {
    return(NULL)
  }
  if (is.list(scope)) {
    return(lapply(scope, function(segment) suppressWarnings(as.numeric(segment))))
  }
  values <- suppressWarnings(as.numeric(scope))
  # A single c(xmin, xmax, ymin, ymax) is replicated to every condition and
  # cannot disagree with itself.
  if (length(values) == 4L || length(values) < 4L * n) {
    return(NULL)
  }
  lapply(seq_len(n), function(i) values[((i - 1L) * 4L + 1L):(i * 4L)])
}

# Equal, or absent on both sides. A bound declared for one condition and left
# NA for another is the quiet form of this mistake: NA means "no theoretical
# bound, use the observed range", so the two end up on different boxes without
# anyone writing two different numbers.
.nsca_same_bound <- function(a, b) {
  if (is.na(a) && is.na(b)) {
    return(TRUE)
  }
  if (is.na(a) || is.na(b)) {
    return(FALSE)
  }
  isTRUE(all.equal(a, b))
}

.nsca_describe_bound <- function(value) {
  if (is.na(value)) "not declared" else format(value)
}

.nsca_validate_scope_y <- function(scope, x_names, y_name) {
  segments <- .nsca_scope_segments(scope, length(x_names))
  if (is.null(segments) || length(segments) < 2L) {
    return(invisible(scope))
  }
  first <- segments[[1L]]
  if (length(first) < 4L) {
    return(invisible(scope))
  }

  mismatches <- character(0L)
  for (i in seq_along(segments)[-1L]) {
    segment <- segments[[i]]
    if (length(segment) < 4L) {
      next
    }
    for (side in 3:4) {
      if (.nsca_same_bound(first[[side]], segment[[side]])) {
        next
      }
      mismatches <- c(mismatches, sprintf(
        "%s %s: %s (%s says %s)",
        y_name, if (side == 3L) "lower" else "upper",
        .nsca_describe_bound(segment[[side]]),
        x_names[[1L]], .nsca_describe_bound(first[[side]])
      ))
    }
  }
  if (length(mismatches) == 0L) {
    return(invisible(scope))
  }

  stop(
    sprintf(
      paste0(
        "'scope' gives '%s' more than one theoretical range. %s.\n",
        "\n",
        "Both components are areas over the same rectangle, which places ",
        "d_nec and d_suf on one scale and allows the joint indices to combine ",
        "them, and the admissible region between them is measured against ",
        "that rectangle as well. Two outcome ranges would place the ",
        "components and the region on different rectangles, so the ",
        "identity\n",
        "  1 - d_nec - d_suf + overlap_share = admissible_region_share\n",
        "would no longer hold, reconstruction_error would report the ",
        "discrepancy, and geometry_acceptable would be FALSE for reasons that ",
        "the resulting message would not trace back to this input.\n",
        "\n",
        "This concerns conditions, not the two components: necessity and ",
        "sufficiency are two corners of one scatter plot and share a scope ",
        "by design.\n",
        "\n",
        "Two possible remedies:\n",
        "  * Specify the same ymin and ymax in every segment. A per-condition ",
        "scope allows conditions on different measurement scales to have ",
        "different X ranges; the outcome half is a repetition required by the ",
        "format, not a separate setting.\n",
        "  * If a condition requires a different outcome scope, it can be ",
        "fitted in a separate call. No information is lost: conditions are ",
        "estimated independently, so one call with several conditions is ",
        "equivalent to several single-condition analyses."
      ),
      y_name, paste(mismatches, collapse = "; ")
    ),
    call. = FALSE
  )
}

# Make `na.action` reach both engines, not only this package's bookkeeping.
#
# The counterpart of SCAtools' `.sca_apply_na_action()`, and needed separately
# because NSCA calls `NCA::nca_analysis()` directly for the necessity side.
# That engine drops incomplete pairs per condition inside its own loop, so a
# full frame is fitted per condition however `na.action` was set. Under
# `"listwise"` that would leave the two components on different row sets from
# each other's bookkeeping and from `raw`, which is the disagreement the rule
# exists to prevent -- and here it would also put `d_nec` and `d_suf` on
# different samples, which is the one thing the joint indices cannot survive.
#
# Only the selected columns are considered: a missing value elsewhere in the
# data frame is not this analysis's business.
.nsca_apply_na_action <- function(frame, x_names, y_name, na.action) {
  if (!identical(na.action, "listwise")) {
    return(frame)
  }
  usable <- lapply(
    c(x_names, y_name),
    function(nm) is.finite(.nsca_numeric_column(frame[[nm]]))
  )
  keep <- Reduce(`&`, usable)
  if (!any(keep)) {
    stop(
      paste0(
        "No row has every selected variable observed, so listwise deletion ",
        "leaves no data to fit. With na.action = \"pairwise\", each ",
        "condition is fitted on its own complete pairs."
      ),
      call. = FALSE
    )
  }
  frame[keep, , drop = FALSE]
}

# One side's influential observations, tagged with which side they came from.
#
# `NCA::nca_outliers()` returns NULL when it finds nothing and prints a message
# rather than raising, so a missing side is ordinary and must not become an
# error. It also drops its `scope` column when a theoretical scope was declared,
# so the two sides can come back with different columns and cannot simply be
# rbind'ed.
.nsca_outliers_side <- function(side, condition, corner, data, x, y, ceiling,
                                scope, k, min.dif, max.results, condensed) {
  result <- NCA::nca_outliers(
    data = data, x = x, y = y, ceiling = ceiling, corner = corner,
    scope = scope, k = k, min.dif = min.dif, max.results = max.results,
    plotly = FALSE, condensed = condensed
  )
  if (is.null(result) || nrow(result) == 0L) {
    return(NULL)
  }
  result <- as.data.frame(result, stringsAsFactors = FALSE)
  cbind(
    data.frame(
      condition = condition, component = side, corner = corner,
      stringsAsFactors = FALSE
    ),
    result
  )
}

# The declared scope belonging to one condition. The counterpart of SCAtools'
# `.sca_scope_for()`, and separate for the same reason the na.action helper is:
# NSCA reaches the necessity engine directly and cannot borrow the other
# package's private function.
.nsca_scope_for <- function(scope, index, n) {
  if (is.null(scope)) {
    return(NULL)
  }
  if (is.list(scope)) {
    segment <- if (index <= length(scope)) scope[[index]] else scope[[1L]]
    return(suppressWarnings(as.numeric(segment)))
  }
  values <- suppressWarnings(as.numeric(scope))
  if (length(values) >= 4L * n && n > 1L) {
    return(values[((index - 1L) * 4L + 1L):(index * 4L)])
  }
  values
}

#' Influential observations on both sides
#'
#' Identifies the individual observations that support the necessity frontier
#' and those that support the sufficiency frontier. This is the counterpart of
#' [SCAtools::sca_outliers()], which addresses the question for one corner.
#'
#' A necessary-and-sufficient claim rests on two frontiers simultaneously, and
#' they are supported by different observations: the necessity frontier by the
#' cases at the top-left of the cloud, the sufficiency frontier by the cases at
#' the bottom-right. A case that is influential for one side is often not
#' influential for the other, so a single list would obscure which claim is
#' fragile. The `component` column identifies the side.
#'
#' @section Two ways to call it:
#' When a fitted model is passed, the conditions, outcome, direction and scope
#' are taken from it, so that the diagnostic cannot describe a different
#' analysis from the one being reported:
#'
#' ```
#' fit <- nsca_analysis(data, x = c("A", "B"), y = "Y", direction = "HH")
#' nsca_outliers(fit)
#' ```
#'
#' Every condition in the model is examined unless `x` selects a subset.
#' Passing `y`, `direction` or `scope` alongside a model is an error rather
#' than an override: the model already fixes all three, and a second
#' specification could place the two components on different corners, which
#' would invalidate a two-sided table.
#'
#' A data frame may be passed instead, as in earlier versions, when no fitted
#' model is available.
#'
#' @section Reading it:
#' The output is best interpreted as a fragility check rather than as a
#' data-cleaning step. An observation that moves a frontier behaves as the
#' method expects, since frontiers are determined by extreme cases rather than
#' by central tendency, so identifying one is not in itself grounds for
#' removing it. The result indicates how much of the claim rests on how few
#' cases, which is worth reporting in either case.
#'
#' @param data A fitted [nsca_analysis()] model, or a data frame.
#' @param x With a model, optional condition names or positions; all conditions
#'   by default. With a data frame, one condition column name.
#' @param y Outcome column name. Taken from the model when one is supplied;
#'   supplying it as well is an error.
#' @param direction One of `"HH"`, `"LH"`, `"HL"`, `"LL"`. Taken from the model
#'   when one is supplied; supplying it as well is an error.
#' @param ceiling Frontier technique, passed to [NCA::nca_outliers()]. With a
#'   model, defaults to the first technique that the model estimated.
#' @param scope Optional theoretical scope `c(xmin, xmax, ymin, ymax)`. Taken
#'   from the model when one is supplied; supplying it as well is an error.
#' @param k Number of observations that may be removed at once.
#' @param min.dif Smallest change in the effect size to report.
#' @param max.results Maximum number of rows to keep per side.
#' @param condensed Passed to [NCA::nca_outliers()].
#' @return A data frame whose first columns name the condition and the
#'   component, with one block per side, or `NULL` when neither side has
#'   anything to report. When only one side does, only that side is returned.
#' @seealso [SCAtools::sca_outliers()], [nsca_corners()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_outliers(fit)
nsca_outliers <- function(data, x = NULL, y = NULL, direction = "HH",
                          ceiling = NULL, scope = NULL, k = 1, min.dif = 1e-2,
                          max.results = 25, condensed = FALSE) {
  if (inherits(data, "nsca_result")) {
    given <- c(
      if (!is.null(y)) "y",
      if (!missing(direction)) "direction",
      if (!is.null(scope)) "scope"
    )
    if (length(given) > 0L) {
      stop(
        sprintf(
          paste0(
            "%s cannot be supplied alongside a fitted model: the model ",
            "already fixes %s, and a second specification could place the two ",
            "components on corners that the fit did not use. The argument can ",
            "be omitted, or nsca_outliers() can be called on the data frame ",
            "instead."
          ),
          paste(sprintf("'%s'", given), collapse = " and "),
          if (length(given) == 1L) "it" else "them"
        ),
        call. = FALSE
      )
    }

    metadata <- .nsca_metadata(data)
    conditions <- .nsca_select_conditions(data, x)
    available <- metadata$conditions
    ceiling <- if (is.null(ceiling)) metadata$ceilings[[1L]] else ceiling
    outcome <- metadata$outcome

    blocks <- list()
    for (condition in conditions) {
      item <- metadata$raw$conditions[[condition]]
      frame <- data.frame(
        as.numeric(item$x), as.numeric(item$y),
        stringsAsFactors = FALSE
      )
      names(frame) <- c(condition, outcome)
      box <- .nsca_scope_for(
        metadata$scope_declared,
        match(condition, available),
        length(available)
      )
      for (side in .nsca_components) {
        corner <- if (identical(side, "necessity")) {
          unname(metadata$necessity_corner[[condition]])
        } else {
          unname(metadata$sufficiency_corner[[condition]])
        }
        blocks[[length(blocks) + 1L]] <- .nsca_outliers_side(
          side, condition, corner, frame, condition, outcome, ceiling, box,
          k, min.dif, max.results, condensed
        )
      }
    }
    blocks <- Filter(Negate(is.null), blocks)
    if (length(blocks) == 0L) {
      return(NULL)
    }
    out <- .nsca_stack_outliers(blocks)
    attr(out, "nsca_direction") <- paste(
      unique(unname(metadata$direction[conditions])), collapse = ", "
    )
    attr(out, "necessity_corner") <- paste(
      unique(unname(metadata$necessity_corner[conditions])), collapse = ", "
    )
    attr(out, "sufficiency_corner") <- paste(
      unique(unname(metadata$sufficiency_corner[conditions])), collapse = ", "
    )
    class(out) <- unique(c("nsca_outliers", class(out)))
    return(out)
  }

  if (is.null(x) || is.null(y)) {
    stop(
      "'x' and 'y' are required when the first argument is a data frame.",
      call. = FALSE
    )
  }
  direction <- .nsca_validate_direction(direction, 1L)
  corners <- nsca_corners(direction)
  frame <- as.data.frame(data, stringsAsFactors = FALSE)
  condition <- if (is.numeric(x)) names(frame)[[x[[1L]]]] else as.character(x)[[1L]]

  necessity <- .nsca_outliers_side(
    "necessity", condition, corners$necessity_corner[[1L]], frame, x, y,
    ceiling, scope, k, min.dif, max.results, condensed
  )
  sufficiency <- .nsca_outliers_side(
    "sufficiency", condition, corners$sufficiency_corner[[1L]], frame, x, y,
    ceiling, scope, k, min.dif, max.results, condensed
  )

  sides <- Filter(Negate(is.null), list(necessity, sufficiency))
  if (length(sides) == 0L) {
    return(NULL)
  }
  if (length(sides) == 2L &&
      !identical(names(sides[[1L]]), names(sides[[2L]]))) {
    # This should be unreachable, which is why it is an error rather than a
    # repair. Both sides are given the same `scope`, `k`, `min.dif` and
    # `condensed`, so both take the same branch of `NCA::nca_outliers()` and
    # return the same columns -- including the `scope` column it drops when a
    # theoretical scope is declared, which it drops for both. If the two ever
    # disagree, the engine's interface has changed underneath this function,
    # and silently keeping the intersection would lose whichever column was
    # new. The same reasoning as the engine contract check: an interface that
    # moves is reported once and loudly.
    stop(
      sprintf(
        paste0(
          "The two sides were returned with different columns, which the ",
          "engine is not expected to produce for identical arguments: ",
          "necessity has %s, sufficiency has %s. This suggests a change in ",
          "the NCA interface rather than a property of the data; please ",
          "report it."
        ),
        paste(names(sides[[1L]]), collapse = ", "),
        paste(names(sides[[2L]]), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out <- do.call(rbind, sides)
  row.names(out) <- NULL
  attr(out, "nsca_direction") <- direction
  attr(out, "necessity_corner") <- corners$necessity_corner[[1L]]
  attr(out, "sufficiency_corner") <- corners$sufficiency_corner[[1L]]
  class(out) <- unique(c("nsca_outliers", class(out)))
  out
}


#' @export
print.nsca_outliers <- function(x, ...) {
  cat(
    "\nNSCA outlier sensitivity: direction ", attr(x, "nsca_direction"),
    "\nnecessity corner ", attr(x, "necessity_corner"),
    "; sufficiency corner ", attr(x, "sufficiency_corner"), "\n",
    sep = ""
  )
  NextMethod()
  invisible(x)
}
# Blocks that need not share columns. See `.sca_stack_outliers()` for why they
# can differ; the model form loops over conditions, and a declared scope may be
# given per condition, so one condition can sit on a declared box and another
# on the observed range.
#
# The single-condition data-frame path keeps its strict check instead: there,
# both sides are given exactly the same arguments, so a disagreement can only
# be an engine interface change and is reported rather than repaired.
.nsca_stack_outliers <- function(blocks) {
  columns <- unique(unlist(lapply(blocks, names), use.names = FALSE))
  filled <- lapply(blocks, function(block) {
    for (column in setdiff(columns, names(block))) {
      block[[column]] <- NA
    }
    block[columns]
  })
  out <- do.call(rbind, filled)
  row.names(out) <- NULL
  out
}

# Refuse an analysis in which no condition varies.
#
# A constant outcome has been refused since 0.4.x -- it leaves the scope with
# no height, so every empty area is measured against zero. A constant
# *condition* does the same thing on the other axis and was not checked, so it
# returned a fit in which `scope_area` was 0, every effect size was NaN and
# every p-value NA. The same degeneracy, refused on one axis and reported as a
# number on the other.
#
# The refusal is only for the case where *every* condition is flat, because
# that is the case where the analysis has nothing to say at all. One flat
# condition among several still leaves the others estimable, and stopping
# would throw those away; those rows keep the coverage warning and their NaN,
# which is the honest answer for that row and not for the call.
.nsca_require_varying_condition <- function(raw) {
  spans <- vapply(raw$conditions, function(item) {
    suppressWarnings(diff(as.numeric(item$x_bounds)))
  }, numeric(1L))
  flat <- !is.finite(spans) | spans <= 0
  if (!all(flat)) {
    return(invisible(raw))
  }
  stop(
    sprintf(
      paste0(
        "No condition varies inside its scope (%s), so both components would ",
        "be areas over a scope of zero width and neither is defined. If the ",
        "condition is bounded by design, a theoretical 'scope' with distinct ",
        "xmin and xmax can be supplied (the result then depends on where the ",
        "constant lies inside that scope); otherwise, it is advisable to ",
        "check the data for a column that did not vary."
      ),
      paste(names(spans), collapse = ", ")
    ),
    call. = FALSE
  )
}

# A declared box with no width, or with its bounds the wrong way round.
#
# The counterpart of SCAtools' `.sca_validate_scope_bounds()`. Both were
# accepted until 0.5.0: the engines repair them silently -- they widen whatever
# they are given until the data fits -- so the analysis was right while
# nsca_scope() reported a declared minimum above its own maximum, which is not
# a scope anybody could put in a methods section.
#
# Refused rather than repaired, on the same grounds as the per-condition
# outcome scope: no legitimate call reaches either shape.
.nsca_validate_scope_bounds <- function(scope, n) {
  if (is.null(scope)) {
    return(invisible(NULL))
  }
  segments <- if (is.list(scope)) {
    lapply(scope, function(s) suppressWarnings(as.numeric(s)))
  } else {
    values <- suppressWarnings(as.numeric(scope))
    if (length(values) == 4L) {
      list(values)
    } else if (length(values) >= 4L * n) {
      lapply(seq_len(n), function(i) values[((i - 1L) * 4L + 1L):(i * 4L)])
    } else {
      list()
    }
  }
  for (i in seq_along(segments)) {
    s <- segments[[i]]
    if (length(s) != 4L) next
    for (pair in list(list(1L, 2L, "x"), list(3L, 4L, "y"))) {
      lo <- s[[pair[[1L]]]]
      hi <- s[[pair[[2L]]]]
      if (is.na(lo) || is.na(hi)) next
      where <- if (length(segments) > 1L) sprintf(" (condition %d)", i) else ""
      if (lo > hi) {
        stop(
          sprintf(
            paste0("'scope' gives the %s axis a minimum of %s above its ",
                   "maximum of %s%s. The engines would swap them without ",
                   "warning."),
            pair[[3L]], format(lo), format(hi), where
          ),
          call. = FALSE
        )
      }
      if (isTRUE(all.equal(lo, hi))) {
        stop(
          sprintf(
            paste0("'scope' gives the %s axis no width (%s to %s)%s, so both ",
                   "components would be areas over a scope of zero width. The ",
                   "engines would widen it to the observed range without ",
                   "warning."),
            pair[[3L]], format(lo), format(hi), where
          ),
          call. = FALSE
        )
      }
    }
  }
  invisible(NULL)
}

# The shapes a declared scope may take, checked before the engine sees it.
#
# Until 0.5.0 a wrong length reached the engines, whose refusal reads "The
# length of scope needs to be equal to the length of x \n or" -- its own line
# break arriving mid-sentence, and no statement of what the legal lengths
# are. A caller with three conditions could not tell whether four numbers
# were wanted or twelve.
#
# Legal: four numbers for one box shared by every condition; four per
# condition; or a list with one entry per condition (or a single entry
# applied to all).
.nsca_validate_scope_length <- function(scope, n) {
  if (is.null(scope)) {
    return(invisible(NULL))
  }
  if (is.list(scope)) {
    widths <- vapply(scope, function(s) length(suppressWarnings(as.numeric(s))),
                     integer(1L))
    if (any(widths != 4L)) {
      stop(
        sprintf(
          paste0("Every entry of a 'scope' list must be four numbers ",
                 "c(xmin, xmax, ymin, ymax); entry %d has %d."),
          which(widths != 4L)[[1L]], widths[widths != 4L][[1L]]
        ),
        call. = FALSE
      )
    }
    if (!(length(scope) %in% c(1L, n))) {
      stop(
        sprintf(
          paste0("'scope' is a list of %d boxes but there %s %d condition%s. ",
                 "Either one box or one box per condition is expected."),
          length(scope), if (n == 1L) "is" else "are", n,
          if (n == 1L) "" else "s"
        ),
        call. = FALSE
      )
    }
    return(invisible(NULL))
  }
  values <- suppressWarnings(as.numeric(scope))
  if (length(values) %in% c(4L, 4L * n)) {
    return(invisible(NULL))
  }
  stop(
    sprintf(
      paste0("'scope' has %d value%s. With %d condition%s, it requires ",
             "either 4 values (c(xmin, xmax, ymin, ymax), one box for all ",
             "conditions) or %d values (four per condition, in the order in ",
             "which the conditions were given)."),
      length(values), if (length(values) == 1L) "" else "s",
      n, if (n == 1L) "" else "s", 4L * n
    ),
    call. = FALSE
  )
}
