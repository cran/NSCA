# How much the joint result moves when the sample does.
#
# The counterpart of SCAtools' sca_resample(), and it reports four quantities
# rather than one because that is the whole point here. Both frontiers are
# fitted to sample extremes, so both components are overstated in a small
# sample -- and the joint indices are built from the two, so they inherit both
# overstatements rather than canceling them. A table that resampled only the
# components would leave the reader to guess what happened to the number the
# package exists to report.
#
# No confidence interval, for the reasons written out in sca_resample(): the
# n-out-of-n bootstrap is inconsistent for boundary estimators, and the
# subsampling distribution is biased upward by construction, so a percentile
# interval read off it would sit above what it claims to cover. Correcting
# needs a convergence rate this package will not assert. What is reported is
# the distribution at each size; the spread is the sampling variability and the
# trend down the `m` column is the bias.

.nsca_resample_quantities <- c(
  "d_nec", "d_suf", "weakest_effect", "balanced_joint_effect"
)

.nsca_resample_one <- function(frame, condition, outcome, direction, ceiling,
                               scope, cutoff, qr.tau, index) {
  blank <- stats::setNames(
    rep(NA_real_, length(.nsca_resample_quantities)),
    .nsca_resample_quantities
  )
  fit <- tryCatch(
    suppressWarnings(suppressMessages(nsca_analysis(
      frame[index, , drop = FALSE], condition, outcome,
      direction = direction, ceilings = ceiling,
      # Pinned to the full sample's box, so only the frontiers move. Letting
      # the scope shrink with the data would put the numerator and the
      # denominator in motion at once.
      scope = scope,
      cutoff = cutoff, qr.tau = qr.tau, test.rep = 0
    ))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(blank)
  }
  row <- tryCatch(suppressWarnings(nsca_table(fit)), error = function(e) NULL)
  if (is.null(row)) {
    return(blank)
  }
  row <- row[row$frontier == ceiling, , drop = FALSE]
  if (nrow(row) == 0L) {
    return(blank)
  }
  vapply(.nsca_resample_quantities, function(q) {
    value <- suppressWarnings(as.numeric(row[[q]][[1L]]))
    if (length(value) == 0L || !is.finite(value)) NA_real_ else value
  }, numeric(1L))
}

#' Variability of the joint result across subsamples
#'
#' Refits one condition on random subsamples of several sizes and reports the
#' distribution of both component effect sizes and both joint indices at each
#' size.
#'
#' @section Why four quantities:
#' Both frontiers are fitted to the sample's own extremes, so both components
#' tend to be overstated when the sample is small. The joint indices are built
#' from the two components, so they inherit both overstatements rather than
#' canceling them: `weakest_effect` takes the smaller component and rises with
#' it, and `balanced_joint_effect` is a geometric mean of two upwardly biased
#' values. Resampling only the components would leave this effect unreported.
#'
#' @section Why there is no confidence interval:
#' The reasons are those set out in [SCAtools::sca_resample()]. In brief, the
#' ordinary n-out-of-n bootstrap is inconsistent for boundary estimators,
#' because a duplicated extreme cannot extend an envelope further outward;
#' subsampling restores consistency but is biased upward, since a smaller
#' sample tends to leave more of the corner empty by chance; and correcting for
#' that bias requires a convergence rate that depends on the technique and the
#' dimension, which this package does not assume.
#'
#' @section Reading it:
#' `observed` is the full-sample value of the quantity, repeated on every row
#' as a reference point. The spread at a given size indicates how much the
#' value depends on the particular cases that were sampled. The trend across
#' the `m` column indicates how much of the value may be attributable to the
#' small size of the sample.
#'
#' These results may be compared with `d_nec_null` and `d_suf_null` from
#' [nsca_table()], which measure the same overstatement differently: by random
#' pairing at the full sample size rather than by refitting at a smaller one.
#'
#' @section Cost:
#' Each replication at each size requires a full two-engine refit, including
#' the band reconstruction. The function is therefore noticeably slower than
#' its SCAtools counterpart, which is why the default for `reps` is lower.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x A condition name or position. Defaults to the first condition.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param sizes Subsample sizes. A value at or below 1 is interpreted as a
#'   fraction of `n`, and a value above 1 as a count.
#' @param reps Subsamples drawn per size.
#' @param seed Optional seed. The caller's random stream is restored
#'   afterward.
#' @return A data frame with one row per subsample size and quantity.
#' @seealso [nsca_table()], [nsca_outliers()], [SCAtools::sca_resample()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' # A report would use the default reps = 100
#' nsca_resample(fit, reps = 10, seed = 1)
nsca_resample <- function(model, x = NULL, ceiling = NULL,
                          sizes = c(0.4, 0.6, 0.8), reps = 100, seed = NULL) {
  metadata <- .nsca_metadata(model)
  conditions <- .nsca_select_conditions(model, x)
  if (is.null(x)) {
    conditions <- conditions[[1L]]
  }
  if (length(conditions) != 1L) {
    stop("A resampling analysis requires exactly one condition.",
         call. = FALSE)
  }
  condition <- conditions[[1L]]
  ceiling <- .nsca_resolve_ceiling(metadata, ceiling)

  reps <- suppressWarnings(as.integer(reps)[[1L]])
  if (!isTRUE(is.finite(reps) && reps >= 1L)) {
    stop("'reps' must be a positive whole number.", call. = FALSE)
  }

  item <- metadata$raw$conditions[[condition]]
  xv <- as.numeric(item$x)
  yv <- as.numeric(item$y)
  keep <- is.finite(xv) & is.finite(yv)
  xv <- xv[keep]
  yv <- yv[keep]
  n <- length(xv)
  if (n < 4L) {
    stop("Too few complete observations to resample.", call. = FALSE)
  }

  frame <- data.frame(xv, yv, stringsAsFactors = FALSE)
  names(frame) <- c(condition, metadata$outcome)
  outcome <- metadata$outcome
  direction <- unname(metadata$direction[[condition]])
  scope <- c(unname(item$x_bounds), unname(item$y_bounds))
  qr.tau <- if (is.null(metadata$qr.tau)) 0.95 else metadata$qr.tau

  full <- suppressWarnings(nsca_table(model))
  full <- full[full$condition == condition &
                 full$frontier == ceiling, , drop = FALSE]
  observed <- vapply(.nsca_resample_quantities, function(q) {
    if (nrow(full) == 0L) NA_real_ else suppressWarnings(as.numeric(full[[q]][[1L]]))
  }, numeric(1L))

  m_values <- .nsca_resample_sizes(sizes, n)

  rows <- .nsca_with_seed(seed, lapply(m_values, function(m) {
    draws <- if (m >= n) 1L else reps
    drawn <- vapply(seq_len(draws), function(i) {
      .nsca_resample_one(
        frame, condition, outcome, direction, ceiling, scope,
        metadata$cutoff, qr.tau, sample.int(n, m)
      )
    }, numeric(length(.nsca_resample_quantities)))
    drawn <- matrix(drawn, nrow = length(.nsca_resample_quantities))
    rownames(drawn) <- .nsca_resample_quantities

    do.call(rbind, lapply(.nsca_resample_quantities, function(q) {
      values <- drawn[q, ]
      usable <- values[is.finite(values)]
      data.frame(
        condition = condition,
        outcome = outcome,
        frontier = ceiling,
        quantity = q,
        n = n,
        m = m,
        fraction = m / n,
        reps = draws,
        converged = length(usable),
        observed = unname(observed[[q]]),
        subsample_median = if (length(usable)) {
          stats::median(usable)
        } else NA_real_,
        subsample_mean = if (length(usable)) mean(usable) else NA_real_,
        subsample_sd = if (length(usable) > 1L) stats::sd(usable) else NA_real_,
        subsample_p05 = if (length(usable)) {
          unname(stats::quantile(usable, 0.05, names = FALSE))
        } else NA_real_,
        subsample_p95 = if (length(usable)) {
          unname(stats::quantile(usable, 0.95, names = FALSE))
        } else NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
  }))

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Duplicated from SCAtools rather than exported from it: the counterpart is
# private there, and a size rule is four lines. Kept identical on purpose --
# a value at or below 1 is a fraction of n, above 1 a count, and both are
# clamped into [3, n] so that a fraction small enough to leave two points
# cannot reach the engine.
.nsca_resample_sizes <- function(sizes, n) {
  sizes <- suppressWarnings(as.numeric(sizes))
  if (length(sizes) == 0L || anyNA(sizes) || any(sizes <= 0)) {
    stop("'sizes' must be one or more positive numbers.", call. = FALSE)
  }
  m <- ifelse(sizes <= 1, round(sizes * n), round(sizes))
  m <- sort(unique(pmin(pmax(m, 3L), n)))
  if (length(m) == 0L) {
    stop("No usable subsample size.", call. = FALSE)
  }
  as.integer(m)
}
