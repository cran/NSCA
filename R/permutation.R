# The shared permutation sequence.
#
# The two component p-values that the engines produce come from two separate
# permutation runs with two separate random draws. That is fine while the two
# tests are only ever combined by an intersection-union rule, which needs no
# assumption about their joint distribution. It stops being fine the moment a
# statistic of both components is tested, because
#
#   W = min(d_nec, d_suf)
#
# has a null distribution that depends on how the two components co-vary under
# random pairing, and that dependence is not free: the two empty zones are
# disjoint subsets of one scope box, so d_nec + d_suf is bounded and the two
# cannot both be large by accident. Drawing them from independent permutations
# throws that constraint away and produces a null distribution for the minimum
# that no dataset could have generated.
#
# So this file drives the permutation itself. One shuffle of the outcome breaks
# the X-Y pairing for every condition at once; both engines then estimate on
# that same shuffled frame; and d_nec, d_suf and their minimum are read off
# together. All three p-values come from one sequence and are mutually
# consistent, which is what a joint statistic requires.
#
# The estimators are the engines' own. Nothing here re-implements a frontier,
# so the null draws are on exactly the same scale as the observed values they
# are compared with. The cost is two full engine fits per replication, which is
# why this is off by default.

# Phipson and Smyth (2010): a permutation p-value must count the observed
# arrangement itself, or it can report zero for something that happened once.
#
# The denominator is the number of draws actually supplied. Until 0.5.0 this
# function filtered out non-finite draws first and then divided by what was
# left, so replications whose engine fits had failed vanished from the
# denominator as well as the numerator: a run asked for 1000 could report a
# p-value computed on 600 without saying so. Failed replications are now
# redrawn by .nsca_shared_null(), and any shortfall is reported there, so the
# array reaching this function is complete by construction. The guard below is
# kept as a floor, not as a silent repair.
.nsca_permutation_p <- function(draws, observed) {
  if (!is.finite(observed) || length(draws) == 0L || anyNA(draws)) {
    return(NA_real_)
  }
  (1 + sum(draws >= observed)) / (length(draws) + 1)
}

# Preserve the caller's random stream. A seed given here is for reproducing
# this test, not for silently reseeding the session.
#
# The restore is unconditional from 0.5.0. Previously it was attached only when
# .Random.seed already existed, so in a fresh session -- where it does not
# exist until the first draw -- set.seed() left its state behind in the global
# environment, which is precisely what this helper documents itself as
# avoiding. Removing the variable is the correct inverse of "there was no seed
# before".
.nsca_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    previous <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", previous, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(".Random.seed", envir = globalenv())),
      add = TRUE
    )
  }
  set.seed(seed)
  force(expr)
}

# How many times a replication may be redrawn before the sequence gives up.
# A failure that is not random -- a scope the engine cannot fit, say -- would
# otherwise loop forever.
.nsca_retry_factor <- 10L

# The sampling rule of NCA's own permutation test, so that a shared-sequence
# p-value and an engine p-value are drawn the same way. Up to this many
# permutations, every ordering of the outcome is distinct; and no more
# permutations are drawn than the n! orderings a sample of n has.
.nsca_distinct_limit <- 720L

# n! for small n, Inf beyond the point where it could matter for any count of
# permutations a user would ask for.
.nsca_orderings <- function(n) {
  if (n <= 12L) factorial(n) else Inf
}

#' Null draws for both components from one permutation sequence
#'
#' @param frame The complete data frame, as supplied to [nsca_analysis()].
#' @param x_names,y_name Condition and outcome column names.
#' @param direction Validated direction vector, one per condition.
#' @param ceilings Frontier techniques.
#' @param scope Theoretical scope, passed through unchanged.
#' @param necessity_corner Unnamed corner vector, one per condition.
#' @param cutoff,qr.tau Engine settings, kept identical to the observed fit.
#' @param reps Number of permutations.
#' @param seed Optional seed for reproducibility.
#' @return A list of four arrays with dimensions
#'   `reps x conditions x ceilings`, holding the necessity, sufficiency,
#'   weakest-effect and component-sum draws. Attributes `requested`, `target`
#'   and `exhaustive` record the requested count, the target count after
#'   capping at `n!`, and whether every ordering of the outcome was
#'   evaluated.
#' @noRd
.nsca_shared_null <- function(frame, x_names, y_name, direction, ceilings,
                              scope, necessity_corner, cutoff, qr.tau,
                              reps, seed = NULL) {
  observations <- nrow(frame)
  requested <- reps
  orderings <- .nsca_orderings(observations)
  if (reps > orderings) {
    reps <- as.integer(orderings)
    message(sprintf(
      paste0(
        "Lowered shared.test.rep to %d: a sample of %d has only %d distinct ",
        "orderings of the outcome, all of which are used."
      ),
      reps, observations, reps
    ))
  }
  # Distinct orderings, as NCA draws them, when the count is small enough for
  # repeats to be likely. Larger counts draw freely, and the stream of draws is
  # then exactly what it was before this rule existed.
  distinct <- reps <= .nsca_distinct_limit
  seen <- character(0L)

  shape <- c(reps, length(x_names), length(ceilings))
  labels <- list(NULL, x_names, ceilings)
  necessity <- array(NA_real_, dim = shape, dimnames = labels)
  sufficiency <- array(NA_real_, dim = shape, dimnames = labels)

  obtained <- 0L
  attempts <- 0L
  max_attempts <- reps * .nsca_retry_factor

  .nsca_with_seed(seed, {
    # Redraw rather than skip. A replication whose engine fits fail used to be
    # left as NA and then dropped from the denominator, so the reported number
    # of replications was an upper bound rather than a count. Here a failure
    # costs another draw, the sequence still delivers `reps` usable
    # replications, and a shortfall against the cap is reported rather than
    # absorbed.
    while (obtained < reps && attempts < max_attempts) {
      ordering <- sample.int(observations)
      if (distinct) {
        # A repeated ordering is drawn again rather than evaluated twice. It
        # does not count as an attempt: attempts measure engine failures.
        # Once every ordering has been seen, there is nothing left to draw.
        if (length(seen) >= orderings) {
          break
        }
        key <- paste(ordering, collapse = ",")
        if (key %in% seen) {
          next
        }
        seen <- c(seen, key)
      }
      attempts <- attempts + 1L
      shuffled <- frame
      shuffled[[y_name]] <- frame[[y_name]][ordering]

      # A permutation of the outcome leaves both marginals unchanged, so the
      # empirical scope of the shuffled frame is identical to the observed
      # one. The areas are therefore measured against the same denominator
      # without having to pin the scope down by hand.
      nec_fit <- tryCatch(
        suppressWarnings(NCA::nca_analysis(
          data = shuffled, x = x_names, y = y_name, ceilings = ceilings,
          corner = necessity_corner, scope = scope,
          bottleneck.x = "actual", bottleneck.y = "actual",
          cutoff = cutoff, qr.tau = qr.tau, effect_aggregation = 1,
          test.rep = 0
        )),
        error = function(e) NULL
      )
      suf_fit <- tryCatch(
        suppressWarnings(SCAtools::sca_analysis(
          data = shuffled, x = x_names, y = y_name, direction = direction,
          ceilings = ceilings, scope = scope, cutoff = cutoff,
          qr.tau = qr.tau, test.rep = 0
        )),
        error = function(e) NULL
      )
      if (is.null(nec_fit) || is.null(suf_fit)) {
        next
      }

      # Collect the whole replication before committing it. A draw counts only
      # when every requested cell came back finite; a partial one is a failure
      # and is redrawn, so "reps" means reps complete replications.
      nec_draw <- matrix(NA_real_, nrow = length(x_names), ncol = length(ceilings))
      suf_draw <- nec_draw
      for (i in seq_along(x_names)) {
        for (k in seq_along(ceilings)) {
          nec_draw[i, k] <- .nsca_nca_param(
            nec_fit, x_names[[i]], ceilings[[k]], "Effect size"
          )
          suf_draw[i, k] <- .nsca_sca_param(
            suf_fit, x_names[[i]], ceilings[[k]], "Effect size"
          )
        }
      }
      if (anyNA(nec_draw) || anyNA(suf_draw)) {
        next
      }

      obtained <- obtained + 1L
      necessity[obtained, , ] <- nec_draw
      sufficiency[obtained, , ] <- suf_draw
    }
  })

  if (obtained < reps) {
    warning(
      sprintf(
        paste0(
          "The shared permutation sequence obtained %d of the %d replications ",
          "requested, after %d attempts. The reported p-values are computed ",
          "on the %d that succeeded. Repeated failure is often not random, so ",
          "it is advisable to check the scope and the frontier technique; a ",
          "null distribution built from a non-random subset of shuffles may ",
          "differ from the intended null."
        ),
        obtained, reps, attempts, obtained
      ),
      call. = FALSE
    )
    if (obtained == 0L) {
      return(NULL)
    }
    keep <- seq_len(obtained)
    necessity <- necessity[keep, , , drop = FALSE]
    sufficiency <- sufficiency[keep, , , drop = FALSE]
    shape <- c(obtained, length(x_names), length(ceilings))
    dimnames(necessity) <- labels
    dimnames(sufficiency) <- labels
  }

  weakest <- array(
    suppressWarnings(pmin(necessity, sufficiency)),
    dim = shape, dimnames = labels
  )
  # The two areas together, from the same shuffle: the null distribution of
  # component_sum, which only a shared sequence can supply.
  total <- array(necessity + sufficiency, dim = shape, dimnames = labels)
  out <- list(necessity = necessity, sufficiency = sufficiency,
              weakest = weakest, sum = total)
  attr(out, "requested") <- requested
  attr(out, "target") <- reps
  attr(out, "exhaustive") <- distinct && obtained == orderings
  out
}

# The Monte Carlo margin of a permutation p-value, as NCA reports it as
# `p-accuracy`: the half-width of a normal-approximation interval at
# `confidence`. Zero when every ordering of the outcome was evaluated, since
# the p-value is then exact rather than estimated.
.nsca_p_accuracy <- function(p, reps, confidence, exhaustive = FALSE) {
  confidence <- suppressWarnings(as.numeric(confidence))
  if (!is.finite(p) || !is.finite(reps) || reps < 1 ||
      length(confidence) != 1L || !is.finite(confidence)) {
    return(NA_real_)
  }
  if (isTRUE(exhaustive)) {
    return(0)
  }
  z <- stats::qnorm(confidence + 0.5 * (1 - confidence))
  z * sqrt(p * (1 - p) / reps)
}

# The 1 - alpha quantile of a null distribution: the value random pairing
# exceeds only alpha of the time, which an observed statistic has to clear for
# its permutation test to reject.
.nsca_critical <- function(draws, alpha) {
  draws <- suppressWarnings(as.numeric(draws))
  draws <- draws[is.finite(draws)]
  alpha <- suppressWarnings(as.numeric(alpha))
  if (length(draws) == 0L || length(alpha) != 1L || !is.finite(alpha)) {
    return(NA_real_)
  }
  unname(stats::quantile(draws, 1 - alpha, names = FALSE))
}

# The null draws behind a joint-table p-value, from the same source as the
# p-value itself: the shared sequence when there is one, otherwise the
# engine's own test of that component. The weakest effect and the component
# sum have a null distribution only under a shared sequence.
.nsca_null_draws <- function(model, shared, side, condition, ceiling) {
  if (!is.null(shared)) {
    return(as.numeric(shared[[side]][, condition, ceiling]))
  }
  if (!(side %in% c("necessity", "sufficiency"))) {
    return(numeric(0L))
  }
  item <- model[[side]]$tests[[condition]]
  if (is.null(item) || is.null(item[[ceiling]])) {
    return(numeric(0L))
  }
  suppressWarnings(as.numeric(item[[ceiling]]$data))
}

.nsca_draw_mean <- function(draws) {
  draws <- draws[is.finite(draws)]
  if (length(draws) == 0L) NA_real_ else mean(draws)
}

# One seed in, two out.
#
# The two components are tested on two independent shuffles. Each engine
# derives its whole permutation sequence from a single draw off the caller's
# random stream, so seeding both with the same number gives both the same
# sequence -- the shared-sequence design, arrived at by accident, with
# p_source still saying two independent engine runs. Drawing two seeds from
# the one the caller gave keeps the sequences independent and keeps both
# reproducible from a single number.
#
# NULL in, NULL out: no seed means no reseeding, which is what the engines'
# own behavior already is.
.nsca_component_seeds <- function(seed) {
  if (is.null(seed)) {
    return(list(sufficiency = NULL, necessity = NULL))
  }
  seed <- suppressWarnings(as.integer(seed))
  if (length(seed) != 1L || is.na(seed)) {
    stop("'test.seed' must be NULL or a single integer.", call. = FALSE)
  }
  drawn <- .nsca_with_seed(seed, sample.int(.Machine$integer.max, 2L))
  list(sufficiency = drawn[[1L]], necessity = drawn[[2L]])
}

# Data with both corners empty.
#
# `NCA::nca_random()` fills one corner: it draws points that leave a single
# frontier clear. A necessary-and-sufficient pattern needs two, on diagonally
# opposite corners, so the cloud has to sit inside a band rather than above or
# below one line. Neither engine offers that, so it is generated here.
#
# Generated for "HH" and then mirrored, which is the same trick the corner
# mapping uses: every direction is "HH" with one or both axes reflected, so
# writing the generator once and reflecting is the only way to guarantee the
# four directions stay each other's mirror images rather than four separate
# implementations that drift.

# The exact area under min(max(a + b * x, 0), 1) over [0, 1], for b >= 0.
.nsca_clipped_integral <- function(intercept, slope) {
  if (slope == 0) {
    return(min(max(intercept, 0), 1))
  }
  x0 <- min(max(-intercept / slope, 0), 1)
  x1 <- min(max((1 - intercept) / slope, 0), 1)
  intercept * (x1 - x0) + slope * (x1^2 - x0^2) / 2 + (1 - x1)
}

# The empty areas a pair of lines leaves on the unit square: above the upper
# (necessity) line and below the lower (sufficiency) line. These are the
# population effect sizes a power analysis aims at. Computed for "HH" and true
# for all four directions: the other three are the same band with one or both
# axes reflected, and a reflection moves an area without changing it.
.nsca_band_truth <- function(lines) {
  list(
    d_suf = .nsca_clipped_integral(lines$lower[[1L]], lines$lower[[2L]]),
    d_nec = 1 - .nsca_clipped_integral(lines$upper[[1L]], lines$upper[[2L]])
  )
}

# The intercept of a line of slope `slope` that leaves `d` empty above it.
#
# One effect size names a family of lines, one for each slope; the slope picks
# the member. In the usual case the line runs from the left edge to the top
# edge of the square, the empty region is a triangle of area
# (1 - a)^2 / (2 * slope), and the intercept follows in closed form. Otherwise
# it is found by root-finding on the exact area, which falls monotonically as
# the intercept rises.
.nsca_intercept_above <- function(d, slope) {
  if (slope == 0) {
    return(1 - d)
  }
  if (d <= min(slope, 1 / slope) / 2) {
    return(1 - sqrt(2 * slope * d))
  }
  stats::uniroot(
    function(a) 1 - .nsca_clipped_integral(a, slope) - d,
    lower = -slope, upper = 1, tol = 1e-12
  )$root
}

# The two boundary lines, each c(intercept, slope), for "HH".
#
# The sufficiency line leaves its area below it, in the lower-right corner. A
# half-turn of the square carries that corner onto the upper-left one and the
# line y = c + b x onto y = (1 - c - b) + b x, so the necessity solution
# serves both sides.
.nsca_simulation_band <- function(d_nec, d_suf, slope_nec, slope_suf) {
  check_d <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value >= 1) {
      stop(sprintf("'%s' must be a single number in [0, 1).", label),
           call. = FALSE)
    }
    as.numeric(value)
  }
  check_slope <- function(value, label) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0) {
      stop(
        sprintf(
          paste0(
            "'%s' must be a single number of at least 0. The steepness of ",
            "the boundary is given as for \"HH\"; the other directions are ",
            "reflections, so their sign follows from `direction`."
          ),
          label
        ),
        call. = FALSE
      )
    }
    as.numeric(value)
  }
  d_nec <- check_d(d_nec, "d_nec")
  d_suf <- check_d(d_suf, "d_suf")
  slope_nec <- check_slope(slope_nec, "slope_nec")
  slope_suf <- check_slope(slope_suf, "slope_suf")

  upper <- c(.nsca_intercept_above(d_nec, slope_nec), slope_nec)
  lower <- c(1 - slope_suf - .nsca_intercept_above(d_suf, slope_suf), slope_suf)

  # Both boundaries are straight lines clipped to the square, so the gap
  # between them is piecewise linear and is smallest at an end of [0, 1] or
  # where one of the lines enters or leaves the square.
  clipped <- function(line, x) pmin(pmax(line[[1L]] + line[[2L]] * x, 0), 1)
  kinks <- function(line) {
    if (line[[2L]] == 0) return(numeric(0))
    (c(0, 1) - line[[1L]]) / line[[2L]]
  }
  x <- c(0, 1, kinks(lower), kinks(upper))
  x <- x[x >= 0 & x <= 1]
  if (any(clipped(lower, x) > clipped(upper, x) + 1e-12)) {
    stop(
      paste0(
        "These effect sizes and slopes leave no band between the two ",
        "boundaries: together the two empty corners would cover more than the ",
        "square allows. A smaller `d_nec` or `d_suf`, or a different slope, ",
        "is required."
      ),
      call. = FALSE
    )
  }
  list(lower = lower, upper = upper)
}

#' Simulate a necessary-and-sufficient pattern
#'
#' Draws observations inside a band, so that both corners that a
#' necessary-and-sufficient claim requires to be empty are empty. This is the
#' counterpart of [SCAtools::sca_random()], which produces a single empty
#' corner because a sufficiency claim requires one.
#'
#' The band is described by what it leaves empty. `d_nec` and `d_suf` are the
#' two population effect sizes: the shares of the unit square left empty in the
#' necessity corner and in the sufficiency corner. Each boundary is a straight
#' line, and since the same area can be left empty by lines of many slopes,
#' `slope_nec` and `slope_suf` specify which line is meant; the intercepts then
#' follow. A slope of 1, the default, is parallel to the diagonal. `d_suf = 0`
#' gives a pattern with no empty sufficiency corner, and `d_nec = 0` one with
#' no empty necessity corner.
#'
#' The condition is drawn uniformly on `[0, 1]`, and at each value of the
#' condition the outcome is drawn uniformly between the two boundaries. Where a
#' boundary leaves the square, the band is truncated at its edge, so the points
#' are uniform across the band at each value of the condition rather than over
#' its whole area.
#'
#' All data are generated for `"HH"` and then reflected. Every direction is
#' `"HH"` with one or both axes mirrored, so a single generator combined with a
#' reflection ensures that the four directions are exact mirror images of each
#' other. Slopes are
#' therefore given as for `"HH"`, as non-negative numbers; for the other
#' directions the reflection supplies the sign. A slope estimated in the units
#' of a study corresponds to `slope * (range of X) / (range of Y)` on the unit
#' square.
#'
#' @param n Number of observations.
#' @param direction One of `"HH"`, `"LH"`, `"HL"`, `"LL"`.
#' @param d_nec,d_suf Population effect sizes of the necessity and the
#'   sufficiency corner, each in `[0, 1)`.
#' @param slope_nec,slope_suf Slopes of the two boundary lines on the unit
#'   square, as for `"HH"`; non-negative.
#' @param names Column names for the condition and the outcome.
#' @return A data frame with two columns, both on `[0, 1]`.
#' @seealso [nsca_power()], [SCAtools::sca_random()]
#' @export
#' @examples
#' dat <- nsca_random(n = 60, d_nec = 0.3, d_suf = 0.3)
#' nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#'
#' # A pattern in which only necessity holds
#' only_nec <- nsca_random(n = 60, d_nec = 0.3, d_suf = 0)
nsca_random <- function(n = 100, direction = "HH",
                        d_nec = 0.10, d_suf = 0.10,
                        slope_nec = 1, slope_suf = 1,
                        names = c("X", "Y")) {
  direction <- .nsca_validate_direction(direction, 1L)
  n <- suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n < 3L) {
    stop("'n' must be a single integer of at least 3.", call. = FALSE)
  }
  lines <- .nsca_simulation_band(d_nec, d_suf, slope_nec, slope_suf)

  x <- stats::runif(n)
  low <- lines$lower[[1L]] + lines$lower[[2L]] * x
  high <- lines$upper[[1L]] + lines$upper[[2L]] * x
  low <- pmin(pmax(low, 0), 1)
  high <- pmin(pmax(high, 0), 1)
  y <- stats::runif(n, min = low, max = high)

  if (substr(direction, 1L, 1L) == "L") {
    x <- 1 - x
  }
  if (substr(direction, 2L, 2L) == "L") {
    y <- 1 - y
  }

  out <- data.frame(x, y)
  names(out) <- names
  attr(out, "nsca_direction") <- direction
  attr(out, "band") <- lines
  out
}

# Dependence between the two component tests, measured rather than assumed.
#
# The obvious summary is `power_joint - power_nec * power_suf`, and it is what
# this function reported first. The product is the joint power the two tests
# would have if they rejected independently, which they cannot: the frontiers
# bound the same cloud and both tests run on the same observations. A
# difference of two rates from a benchmark known to be false is not a
# dependence measure, though its name said it was.
#
# The 2x2 table of who rejected has a correlation, and that is one. phi is the
# Pearson correlation of the two binary indicators: zero under independence,
# positive when the components tend to succeed together, negative when they
# trade off. Undefined when either test rejected always or never, which is
# exactly when there is nothing to correlate.
.nsca_rejection_phi <- function(nec, suf) {
  if (length(nec) < 2L) {
    return(NA_real_)
  }
  n11 <- sum(nec & suf)
  n10 <- sum(nec & !suf)
  n01 <- sum(!nec & suf)
  n00 <- sum(!nec & !suf)
  denominator <- sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  (n11 * n00 - n10 * n01) / denominator
}

#' Power for a necessary-and-sufficient claim
#'
#' Simulates from [nsca_random()] and reports how often each component rejects,
#' how often both reject together, and how far the realized effect sizes
#' differ from the population values.
#'
#' @section Describing the design:
#' `d_nec` and `d_suf` are the population effect sizes of the two corners, and
#' `slope_nec` and `slope_suf` the slopes of the two boundary lines on the unit
#' square; see [nsca_random()]. The effect sizes may come from theory, from a
#' pilot study or from comparable studies. A slope estimated in the units of a
#' study corresponds to `slope * (range of X) / (range of Y)` on the unit
#' square. A power analysis is a planning tool: values estimated from the data
#' of the study itself provide little information about that study's power.
#'
#' @section Three power estimates:
#' `power_necessity` and `power_sufficiency` are the rejection rates of the two
#' component tests; `power_joint` is the rate at which the intersection-union
#' test rejects, which requires both. The joint rate is the one for which a
#' necessary-and-sufficient study is powered, and it cannot exceed either of
#' the two component rates.
#'
#' `power_if_independent` is the product of the two component rates: the joint
#' power the tests would have if they rejected independently. The tests are not
#' independent, since the two frontiers bound the same cloud and both tests use
#' the same observations, so the product serves only as a reference point
#' whose independence assumption does not hold. `rejection_phi` measures the
#' departure from it directly.
#' It is the correlation of the two rejection indicators across replications:
#' zero under independence, positive when the components tend to reject
#' together, negative when they trade off, and `NA` when one of them rejected
#' in every replication or in none.
#'
#' @section The two supplementary tests:
#' With `shared = TRUE` every replication draws one shared permutation
#' sequence (see [nsca_analysis()]) instead of the engines' two separate runs,
#' all five p-values come from it, and two more rates are reported:
#' `power_weakest_perm` and `power_sum_perm`, how often `p_weakest_perm` and
#' `p_sum_perm` reject. With `shared = FALSE` both columns are `NA`.
#'
#' They address the questions that those p-values address (see
#' [nsca_table()]), not the joint claim, and their rates are best interpreted
#' with that in mind. In our simulations with both corners equally empty, both
#' supplementary tests rejected more often than the joint test. With one corner
#' substantially empty and the other with an effect size close to zero,
#' `power_sum_perm` remained high, and `power_weakest_perm` exceeded the rate
#' at which the weaker corner's own test rejected. One largely empty corner can
#' be enough to make the sum large, and the smaller of the two areas is
#' compared with a null distribution that lies below the null of either area
#' alone, since random pairing rarely leaves both corners empty at the same
#' time. These are among the reasons why `joint_support` uses `p_nsca_iut`.
#'
#' @section Checking a design in which one side is absent:
#' Every p-value in this package compares the data with random pairing of the
#' same values (see [nsca_table()]). When one corner is substantially empty,
#' random pairing can fill the opposite corner more than the data do, so the
#' test of that corner may reject more often than its nominal level even where
#' the corner is not empty. In our simulations this tendency was noticeable,
#' and its size depended on the strength of the other side, the sample size and
#' the distributions involved. The scripts behind these statements, with their
#' results, are installed in `system.file("simulations", package = "NSCA")`.
#'
#' Researchers who wish to examine this in their own design can simulate it
#' with one side removed, setting the side under examination to zero and the
#' other to the expected value, for example
#' `nsca_power(n = 60, d_nec = 0.25, d_suf = 0, rep = 500, test.rep = 999)`.
#' In such a design the claim that both sides hold is false, so
#' `power_sufficiency` and `power_joint` (and, with `shared = TRUE`,
#' `power_weakest_perm` and `power_sum_perm`) estimate false-positive rates
#' rather than power, while `power_necessity` remains the power of the side
#' that holds. `d_nec = 0` examines the other side in the same way.
#'
#' Rates noticeably above `alpha` suggest that the corresponding p-value is
#' best interpreted together with the effect sizes, the scatter plot and theory
#' rather than in isolation. A possible further step is to find the smaller
#' `alpha` at which the rate in the side-removed design falls to the intended
#' level; how well that calibration transfers depends on how closely the
#' simulated design resembles the data. The simulated design is a
#' simplification (a uniform condition and an outcome uniform between straight
#' boundaries), so its rates are best read as indications for similar
#' designs.
#'
#' @section Precision of the simulated rates:
#' Each rate is a proportion over the replications, so it carries Monte Carlo
#' error of about `sqrt(rate * (1 - rate) / rep)`: near 0.05 that is about 0.01
#' with `rep = 500` and 0.007 with `rep = 1000`. `test.rep` sets how finely
#' each p-value is resolved; 999 or more is advisable for a table that will be
#' reported.
#'
#' @section The realized effect size exceeds the population one:
#' The two boundary lines are the population frontiers, so `d_nec_declared` and
#' `d_suf_declared` are the true effect sizes. A finite sample does not
#' reproduce them, and the deviation has a systematic direction: the CE-FDH
#' ceiling is the running maximum of the sample, which lies below the
#' population boundary, so the empty region above it is larger than the true
#' one. The same argument applies, in reverse, to the floor. Both empty areas
#' are therefore overstated, and the overstatement decreases as n grows: in one
#' simulation with both effect sizes at 0.32, the realized necessity effect
#' size averaged about 0.43 at n = 40 and 0.37 at n = 200.
#'
#' `d_nec_realized` and `d_suf_realized` report this quantity, so that the size
#' of the bias appears in the same table as the corresponding power.
#'
#' This is the same overstatement that [nsca_table()] reports as `d_nec_null`
#' and `d_suf_null`, measured differently. There it is obtained from the
#' permutation null on one observed sample, that is, from what random pairing
#' alone produces at that n, scope and technique. Here it is the gap between a
#' population of known effect size and what fitting recovers from a finite
#' sample drawn from it. A simulation can compare against a known true value;
#' an analysis cannot, which is why the analysis uses the null instead. In our
#' simulations the two agreed in direction and in order of magnitude, and both
#' decreased as n grew.
#'
#' No claim is made about what this implies for power. The observed statistic
#' is inflated, and so is the permutation null, which the same estimator
#' computes on shuffled data; the extent to which the inflation cancels is an
#' empirical question that this function does not address. The gap is best
#' read as an indication of the estimator's small-sample behavior under this
#' design, and a power computed where the gap is large rests on a statistic
#' that may not yet measure what it measures in larger samples.
#'
#' @section Cost:
#' Every replication fits two engines and runs two permutation tests, so the
#' cost is `rep * length(n) * test.rep * 2`. We suggest starting with
#' `rep = 50` and `test.rep = 100` to see the general shape, and increasing
#' both only for the values that will be reported. The shared sequence refits
#' both engines on every shuffle, which takes two to three times as long as
#' the engines' own tests.
#'
#' The permutation test cannot return a p-value below `1 / test.rep`
#' (`1 / (test.rep + 1)` with `shared = TRUE`), so with fewer than `1 / alpha`
#' permutations (20 at `alpha = 0.05`) no rejection is possible, and the
#' function warns that the estimated power would then be 0 regardless of the
#' design.
#'
#' @param n Sample sizes to evaluate.
#' @param direction One of `"HH"`, `"LH"`, `"HL"`, `"LL"`.
#' @param d_nec,d_suf Population effect sizes of the necessity and the
#'   sufficiency corner, each in `[0, 1)`, passed to [nsca_random()]. The
#'   default of 0.10 on each side matches the default effect size of
#'   [SCAtools::sca_power()].
#' @param slope_nec,slope_suf Slopes of the two boundary lines on the unit
#'   square, as for `"HH"`; non-negative.
#' @param ceiling Frontier technique.
#' @param alpha Significance level for both component tests.
#' @param rep Replications per sample size.
#' @param test.rep Permutation replications within each fit.
#' @param seed Optional seed. The caller's random stream is restored.
#' @param shared Whether to use one shared permutation sequence of `test.rep`
#'   shuffles per replication, and also to report the power of
#'   `p_weakest_perm` and `p_sum_perm`.
#' @return A data frame with one row per sample size. `power_weakest_perm`
#'   and `power_sum_perm` are `NA` unless `shared = TRUE`.
#' @seealso [nsca_random()], [SCAtools::sca_power()]
#' @export
#' @examples
#' # Kept small for speed: a run intended for reporting uses the defaults and
#' # takes minutes rather than seconds.
#' power <- nsca_power(n = c(30, 60), d_nec = 0.25, d_suf = 0.25, rep = 2,
#'                     test.rep = 39, seed = 1)
#' power
nsca_power <- function(n = c(50, 100, 200), direction = "HH",
                       d_nec = 0.10, d_suf = 0.10,
                       slope_nec = 1, slope_suf = 1,
                       ceiling = "ce_fdh", alpha = 0.05,
                       rep = 50, test.rep = 100, seed = NULL,
                       shared = FALSE) {
  direction <- .nsca_validate_direction(direction, 1L)
  n <- suppressWarnings(as.integer(n))
  if (length(n) == 0L || anyNA(n) || any(n < 3L)) {
    stop("'n' must be one or more integers of at least 3.", call. = FALSE)
  }
  if (length(ceiling) != 1L) {
    stop("A power analysis requires exactly one frontier technique.",
         call. = FALSE)
  }
  rep <- suppressWarnings(as.integer(rep))
  if (length(rep) != 1L || is.na(rep) || rep < 1L) {
    stop("'rep' must be a single positive integer.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  test.rep <- .nsca_validate_test_rep(test.rep, "test.rep")
  if (test.rep < 1L) {
    stop("'test.rep' must be at least 1: power is the rate at which a test ",
         "rejects, and without a permutation test no rejection is possible.",
         call. = FALSE)
  }
  if (!is.logical(shared) || length(shared) != 1L || is.na(shared)) {
    stop("'shared' must be TRUE or FALSE.", call. = FALSE)
  }
  # Validated once, here, so a malformed band fails before the first fit
  # rather than on replication one of sample size one.
  lines <- .nsca_simulation_band(d_nec, d_suf, slope_nec, slope_suf)
  truth <- .nsca_band_truth(lines)
  # A test that cannot return a p-value at or below alpha makes every power
  # 0, which the table would report as a property of the design. Said before
  # the run, which takes minutes, rather than after it.
  .nsca_warn_power_floor(test.rep, alpha, n, shared)

  .nsca_with_seed(seed, {
    rows <- lapply(n, function(size) {
      nec <- logical(0)
      suf <- logical(0)
      weakest <- logical(0)
      total <- logical(0)
      nec_d <- numeric(0)
      suf_d <- numeric(0)
      attempts <- 0L
      # Redrawn, capped, and the shortfall reported: the same three parts as
      # the shared permutation sequence, for the same reason. Redrawing alone
      # would be wrong, because if the replications that fail are the ones
      # where the band is tightest then quietly replacing them biases every
      # rate. The cap stops a design the engines cannot fit from looping
      # forever, and the shortfall is what tells the reader which case they
      # are in.
      max_attempts <- rep * .nsca_retry_factor

      while (length(nec) < rep && attempts < max_attempts) {
        attempts <- attempts + 1L
        # Two seeds per replication, drawn before anything else touches the
        # stream: one for the data, one for the tests. The engines reseed the
        # session from the clock when a permutation test finishes, so until
        # 0.5.2 every replication after the first drew its data from an
        # unseeded stream and `seed` reproduced nothing past it.
        seeds <- sample.int(.Machine$integer.max, 2L)
        dat <- .nsca_with_seed(
          seeds[[1L]],
          nsca_random(size, direction = direction, d_nec = d_nec,
                      d_suf = d_suf, slope_nec = slope_nec,
                      slope_suf = slope_suf)
        )
        # The engines print progress for every test; a simulation running
        # hundreds of them keeps the console for its own result.
        utils::capture.output(fit <- tryCatch(
          suppressWarnings(suppressMessages(nsca_analysis(
            dat, names(dat)[[1L]], names(dat)[[2L]],
            direction = direction, ceilings = ceiling,
            # The unit square, so the realized effect sizes are areas over the
            # same box as the declared ones. Without it the denominator would
            # be the observed range and the two would not be comparable.
            scope = c(0, 1, 0, 1),
            test.rep = if (shared) 0L else test.rep,
            shared.test.rep = if (shared) test.rep else 0L,
            test.p_threshold = alpha,
            test.seed = seeds[[2L]], shared.seed = seeds[[2L]]
          ))),
          error = function(e) NULL
        ))
        if (is.null(fit)) {
          next
        }
        row <- suppressWarnings(nsca_table(fit))
        row <- row[row$frontier == ceiling, , drop = FALSE]
        if (nrow(row) == 0L ||
            !is.finite(row$p_nec[[1L]]) || !is.finite(row$p_suf[[1L]])) {
          next
        }
        if (shared && (!is.finite(row$p_weakest_perm[[1L]]) ||
                       !is.finite(row$p_sum_perm[[1L]]))) {
          next
        }
        nec <- c(nec, row$p_nec[[1L]] <= alpha)
        suf <- c(suf, row$p_suf[[1L]] <= alpha)
        if (shared) {
          weakest <- c(weakest, row$p_weakest_perm[[1L]] <= alpha)
          total <- c(total, row$p_sum_perm[[1L]] <= alpha)
        }
        nec_d <- c(nec_d, row$d_nec[[1L]])
        suf_d <- c(suf_d, row$d_suf[[1L]])
      }

      obtained <- length(nec)
      rate <- function(v) if (obtained == 0L) NA_real_ else mean(v)
      power_nec <- rate(nec)
      power_suf <- rate(suf)

      data.frame(
        n = size,
        direction = direction,
        frontier = ceiling,
        alpha = alpha,
        replications = obtained,
        attempts = attempts,
        shortfall = rep - obtained,
        d_nec_declared = truth$d_nec,
        d_suf_declared = truth$d_suf,
        d_nec_realized = if (obtained == 0L) NA_real_ else mean(nec_d),
        d_suf_realized = if (obtained == 0L) NA_real_ else mean(suf_d),
        power_necessity = power_nec,
        power_sufficiency = power_suf,
        power_joint = rate(nec & suf),
        power_if_independent = power_nec * power_suf,
        rejection_phi = .nsca_rejection_phi(nec, suf),
        power_weakest_perm = if (shared) rate(weakest) else NA_real_,
        power_sum_perm = if (shared) rate(total) else NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    out <- do.call(rbind, rows)
    row.names(out) <- NULL

    if (any(out$shortfall > 0L)) {
      warning(
        sprintf(
          paste0(
            "%d of %d replications could not be obtained after redrawing up ",
            "to %d times per sample size. The rates are computed over the ",
            "replications that were obtained, and `replications` reports ",
            "their number. A large shortfall indicates that the design is at ",
            "the limit of what the engines can fit, not that power is low; if ",
            "the draws that fail are those where the band is narrowest, the ",
            "rates may be biased."
          ),
          sum(out$shortfall), length(n) * rep, .nsca_retry_factor
        ),
        call. = FALSE
      )
    }
    out
  })
}
