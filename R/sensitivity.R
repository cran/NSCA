# How much of the joint result is the box.
#
# SCAtools has sca_scope_sensitivity(); NSCA had nsca_scope(), which reports
# the box, and nothing that reports what changing it would do. That gap is
# worse here than there, because in a one-sided analysis the scope is the
# denominator of one number and here it is the denominator of *both* -- and
# therefore of every joint index built from them.
#
# The consequence is not obvious until it is printed, and it runs the opposite
# way to the first guess. Widening the box does not dilute the empty areas: the
# frontiers stay where the data is, so every square unit added beyond the cloud
# lies on the empty side of one frontier or the other. The empty areas grow
# faster than the box, d_nec and d_suf both *rise*, and
# admissible_region_share -- one minus the two, plus the overlap -- falls
# toward 0. In the limit the data is a speck in a large box and almost all of
# that box is forbidden by one claim or the other.
#
# Measured on one 80-case sample, observed range to twice the span:
#
#   widen  d_nec  d_suf  overlap  reconstruction_error  geometry_acceptable
#   0.00    .353   .471  .00008                .00005  TRUE
#   0.25    .462   .515  .00003                .05555  FALSE
#   0.50    .526   .555  .00002                .12502  FALSE
#   1.00    .595   .608  .00001                .22221  FALSE
#
# The second half of the finding, and it is not what it looks like. The
# geometry check fails on `reconstruction_error`, not on overlap: the two
# frontiers never come close to crossing (the overlap *falls*, to one part in
# a hundred thousand). What breaks is the agreement between the engines and
# this package's own reconstruction of the frontiers -- outside the observed
# range of X there are no data to fix a step frontier, so the engine's
# extension and the 2001-point reconstruction here diverge, and the identity
#
#   1 - d_nec - d_suf + overlap_share = admissible_region_share
#
# stops closing. That is precisely what reconstruction_error exists to catch,
# working as designed.
#
# The distinction matters for how a flagged row should be read. Those rows are
# not known to be wrong -- they are *unverified*: the package can no longer
# confirm that the areas it reports describe the frontiers the engines fitted.
# So the series carries overlap_share and reconstruction_error beside the flag,
# because "geometry failed" without saying which criterion failed is the kind
# of summary that gets read as whichever failure the reader already expected.
#
# A relationship can therefore be written up as loosely bounded on one box and
# as almost entirely forbidden on another, with the same data and the same two
# frontiers. Nobody estimated that difference; somebody chose it.
#
# As in SCAtools: the series reports the numbers on each box and does not
# subtract any two of them, because which row is the baseline is the caller's
# question. The observed range answers "how much of this is box, counting from
# the smallest box that holds the data"; the box the model was fitted on
# answers "how sensitive is my number to the box I chose". Different questions,
# and the table cannot silently pick one.

.nsca_widened_scope <- function(x, y, fraction) {
  x_span <- diff(range(x))
  y_span <- diff(range(y))
  c(
    min(x) - fraction * x_span, max(x) + fraction * x_span,
    min(y) - fraction * y_span, max(y) + fraction * y_span
  )
}

.nsca_sensitivity_row <- function(metadata, condition, ceiling, scope, label) {
  item <- metadata$raw$conditions[[condition]]
  outcome <- metadata$outcome
  frame <- data.frame(item$x, item$y)
  names(frame) <- c(condition, outcome)

  fit <- tryCatch(
    suppressWarnings(suppressMessages(nsca_analysis(
      frame, condition, outcome,
      direction = unname(metadata$direction[[condition]]),
      ceilings = ceiling, scope = scope,
      # The engine settings that move a frontier have to come from the model
      # being tested, or the series compares boxes and techniques at once.
      cutoff = metadata$cutoff,
      qr.tau = if (is.null(metadata$qr.tau)) 0.95 else metadata$qr.tau,
      test.rep = 0
    ))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(NULL)
  }
  row <- tryCatch(suppressWarnings(nsca_table(fit)), error = function(e) NULL)
  if (is.null(row)) {
    return(NULL)
  }
  row <- row[row$frontier == ceiling, , drop = FALSE]
  if (nrow(row) == 0L) {
    return(NULL)
  }

  data.frame(
    condition = condition,
    outcome = outcome,
    frontier = ceiling,
    scope = label,
    xmin = scope[[1L]], xmax = scope[[2L]],
    ymin = scope[[3L]], ymax = scope[[4L]],
    scope_area = row$scope_area[[1L]],
    d_nec = row$d_nec[[1L]],
    d_suf = row$d_suf[[1L]],
    weakest_effect = row$weakest_effect[[1L]],
    balanced_joint_effect = row$balanced_joint_effect[[1L]],
    admissible_region_share = row$admissible_region_share[[1L]],
    # Beside the verdict, the two numbers that produce it. A widened box fails
    # on reconstruction_error while overlap_share falls, so the flag alone
    # would be read as the wrong failure.
    overlap_share = row$overlap_share[[1L]],
    reconstruction_error = row$reconstruction_error[[1L]],
    geometry_acceptable = row$geometry_acceptable[[1L]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Sensitivity of the joint result to the scope
#'
#' Refits one condition on a series of scopes and reports both components and
#' both joint indices for each. This is the counterpart of
#' [SCAtools::sca_scope_sensitivity()]; the scope has a greater influence in a
#' two-sided analysis, as described below.
#'
#' @section Why the scope matters more in a two-sided analysis:
#' In a one-sided analysis the scope is the denominator of one effect size. In
#' a two-sided analysis it is the denominator of both, and therefore of every
#' index built from them.
#'
#' The direction of the change may be contrary to expectation. Widening the box
#' does not dilute the empty areas: the frontiers remain where the data are, so
#' every square unit added beyond the cloud falls on the empty side of one
#' frontier or the other. The empty areas grow faster than the box, `d_nec`
#' and `d_suf` both rise, `weakest_effect` and `balanced_joint_effect` rise
#' with them, and `admissible_region_share` falls toward 0. In the limit the
#' data occupy a very small part of a large box, and almost all of that box is
#' excluded by one claim or the other.
#'
#' The same data and the same two frontiers can therefore be described as a
#' loosely bounded relationship or as an almost entirely constrained one,
#' depending on a box that was chosen rather than estimated. This table allows
#' researchers to assess which of these descriptions their own result is closer
#' to.
#'
#' @section Interpreting geometry_acceptable and the criterion that failed:
#' The flag can become `FALSE` at modest amounts of widening, and the reason may
#' differ from the one first assumed. On one 80-case sample:
#'
#' ```
#' widen  d_nec  d_suf  overlap_share  reconstruction_error  geometry
#'  0.00   .353   .471         .00008                .00005  TRUE
#'  0.25   .462   .515         .00003                .05555  FALSE
#'  0.50   .526   .555         .00002                .12502  FALSE
#'  1.00   .595   .608         .00001                .22221  FALSE
#' ```
#'
#' The two frontiers do not cross: `overlap_share` decreases, to one part in a
#' hundred thousand. The criterion that fails is `reconstruction_error`,
#' against its 0.02 tolerance. Outside the observed range of X there are no
#' data to determine a step frontier, so the engines' extension of it and this
#' package's own reconstruction diverge, and the identity
#' \eqn{1 - d_{nec} - d_{suf} + \mathrm{overlap} = \mathrm{admissible}}
#' no longer holds. Detecting this is the purpose of `reconstruction_error`.
#'
#' A flagged row should therefore be regarded as unverified rather than as
#' incorrect: the package can no longer confirm that the areas it reports
#' describe the frontiers the engines fitted. Both diagnostic columns are
#' included in the table so that the reader can see which criterion was not
#' met, rather than interpreting a geometry failure as whichever failure was
#' expected.
#'
#' @section What the boxes are:
#' By default, each box is the observed range widened by a fraction of the
#' observed span on each side, because such a quantity is easy for a reader to
#' judge: "a box half as wide again as the data" is more interpretable than a
#' pair of numbers. With the default `widen`, the first row is the observed
#' range itself, which is the box the engines use when no scope is declared and
#' the smallest box that contains the data, since any smaller box is widened
#' back to it.
#'
#' When the model was fitted on a declared scope, that box is added as a final
#' row labeled `"as fitted (declared)"`, so that the user can see where the
#' chosen scope lies among the alternatives without having to judge it by eye.
#'
#' @section Why no derived share is reported:
#' No column reports a change from a baseline, because the choice of baseline
#' depends on the user's question. Taking the observed range as the baseline
#' addresses how much of the result is attributable to the box, counting from
#' the smallest box that contains the data; taking the scope on which the model
#' was fitted addresses how sensitive the result is to the chosen box. These
#' are different questions, and the table does not choose between them. The
#' data do not change between rows, so every difference within a column
#' reflects the box.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x A condition name or position. Defaults to the first condition.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param widen Fractions of the observed span to add to each side of each
#'   axis. `0` is the observed range.
#' @param scopes Optional list of explicit `c(xmin, xmax, ymin, ymax)` boxes,
#'   used instead of `widen` when specific theoretical scopes are of interest.
#' @return A data frame with one row per scope.
#' @seealso [nsca_scope()], [nsca_table()], [SCAtools::sca_scope_sensitivity()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_scope_sensitivity(fit)
nsca_scope_sensitivity <- function(model, x = NULL, ceiling = NULL,
                                   widen = c(0, 0.25, 0.5, 1),
                                   scopes = NULL) {
  metadata <- .nsca_metadata(model)
  conditions <- .nsca_select_conditions(model, x)
  if (is.null(x)) {
    conditions <- conditions[[1L]]
  }
  if (length(conditions) != 1L) {
    stop("A scope sensitivity analysis requires exactly one condition.",
         call. = FALSE)
  }
  condition <- conditions[[1L]]
  ceiling <- .nsca_resolve_ceiling(metadata, ceiling)

  item <- metadata$raw$conditions[[condition]]
  if (is.null(item)) {
    stop(sprintf("Unknown condition '%s'.", condition), call. = FALSE)
  }

  if (is.null(scopes)) {
    widen <- suppressWarnings(as.numeric(widen))
    if (length(widen) == 0L || anyNA(widen) || any(widen < 0)) {
      stop("'widen' must be one or more non-negative numbers.", call. = FALSE)
    }
    widen <- sort(unique(widen))
    scopes <- lapply(widen, function(f) .nsca_widened_scope(item$x, item$y, f))
    labels <- ifelse(
      widen == 0,
      "observed range",
      sprintf("observed + %g%% of span per side", 100 * widen)
    )
    declared <- metadata$scope_declared
    if (identical(metadata$scope_source, "declared") && !is.null(declared)) {
      box <- .nsca_scope_for(
        declared, match(condition, metadata$conditions),
        length(metadata$conditions)
      )
      if (length(box) == 4L && all(is.finite(box))) {
        scopes <- c(scopes, list(box))
        labels <- c(labels, "as fitted (declared)")
      }
    }
  } else {
    if (!is.list(scopes) || length(scopes) == 0L) {
      stop("'scopes' must be a non-empty list of c(xmin, xmax, ymin, ymax).",
           call. = FALSE)
    }
    scopes <- lapply(scopes, function(s) suppressWarnings(as.numeric(s)))
    bad <- vapply(scopes, function(s) length(s) != 4L || anyNA(s), logical(1L))
    if (any(bad)) {
      stop("Every entry of 'scopes' must be four finite numbers.",
           call. = FALSE)
    }
    labels <- names(scopes)
    if (is.null(labels) || !all(nzchar(labels))) {
      labels <- vapply(
        scopes,
        function(s) paste(format(s), collapse = ", "),
        character(1L)
      )
    }
  }

  rows <- lapply(seq_along(scopes), function(i) {
    .nsca_sensitivity_row(metadata, condition, ceiling, scopes[[i]],
                          labels[[i]])
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out) || nrow(out) == 0L) {
    stop(
      sprintf("No scope in the series could be fitted for '%s'.", condition),
      call. = FALSE
    )
  }
  row.names(out) <- NULL
  out
}

#' Normalize numeric variables for NSCA
#'
#' Delegated to [SCAtools::sca_normalize()], which wraps the engine's own
#' utility. It is exported from this package so that an analysis does not need
#' to attach the other package to rescale its inputs.
#'
#' Both components of an NSCA model are fitted to the same two columns, so a
#' transformation applied here applies to both. The effect sizes are not
#' invariant to a monotone rescaling of either variable: normalizing changes
#' the values, and it is therefore advisable to report which normalization was
#' used.
#'
#' @param data Input data frame or matrix.
#' @param scale Output minimum-maximum pair, or one pair per column.
#' @param min_max Optional empirical or theoretical minimum-maximum pair, or
#'   one pair per column.
#' @return Data on the requested scales.
#' @seealso [SCAtools::sca_normalize()], [nsca_rescale()]
#' @export
#' @examples
#' dat <- data.frame(X = c(0, 2, 4, 10), Y = c(10, 20, 30, 50))
#' nsca_normalize(dat)
nsca_normalize <- function(data, scale = NULL, min_max = NULL) {
  SCAtools::sca_normalize(data = data, scale = scale, min_max = min_max)
}
