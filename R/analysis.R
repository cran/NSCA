.nsca_metadata <- function(model) {
  metadata <- attr(model, "nsca")
  if (is.null(metadata)) {
    stop("The supplied object is not a valid 'nsca_result'.", call. = FALSE)
  }
  metadata
}

# Lines that produce no empty-zone estimate, and therefore no component of a
# necessary-and-sufficient claim.
#
# "ols" is an ordinary least-squares regression of the outcome on the
# condition, fitted to all the data: the line lm(y ~ x) returns. It summarizes
# central tendency, not an empty space, so it can be neither the necessity nor
# the sufficiency component. Since 0.4.1 it is still available, through
# `reference`, as a comparison line: section 6.1 of the condition analysis in
# degree framework treats average-effect and condition analysis as different
# targets that are most informative shown together rather than merged into one
# number. Drawing it changes no NSCA statistic.
.nsca_unusable_ceilings <- c("ols")
.nsca_reference_lines <- c("ols", "qr")

.nsca_validate_reference <- function(reference) {
  if (is.null(reference)) {
    return(character(0L))
  }
  reference <- as.character(reference)
  reference <- unique(tolower(trimws(reference[!is.na(reference)])))
  reference <- setdiff(reference, c("none", ""))
  if (length(reference) == 0L) {
    return(character(0L))
  }
  unknown <- setdiff(reference, .nsca_reference_lines)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown reference line(s): %s. Available: %s, or NULL for none.",
        paste(unknown, collapse = ", "),
        paste(.nsca_reference_lines, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  reference
}

.nsca_validate_ceilings <- function(ceilings) {
  if (!is.character(ceilings) || length(ceilings) == 0L) {
    stop("'ceilings' must name at least one frontier technique.", call. = FALSE)
  }
  ceilings <- unique(tolower(trimws(ceilings)))
  unknown <- setdiff(ceilings, .nsca_known_ceilings)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        paste0(
          "Unknown frontier technique(s) in 'ceilings': %s. Available: %s. ",
          "In earlier versions an unrecognized name was passed to the ",
          "engines, which then used their own defaults while the recorded ",
          "settings still listed the requested name."
        ),
        paste(sprintf("'%s'", unknown), collapse = ", "),
        paste(.nsca_known_ceilings, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  dropped <- intersect(ceilings, .nsca_unusable_ceilings)
  if (length(dropped) > 0L) {
    warning(
      sprintf(
        paste0(
          "Ignoring %s in 'ceilings': it is a central-tendency line, not an ",
          "empty-space frontier, so it yields no necessity or sufficiency ",
          "component. It can be drawn for comparison with reference = \"ols\"."
        ),
        paste(dropped, collapse = ", ")
      ),
      call. = FALSE
    )
    ceilings <- setdiff(ceilings, dropped)
  }
  if (length(ceilings) == 0L) {
    stop("No usable frontier technique remains.", call. = FALSE)
  }
  ceilings
}

# Practical-relevance thresholds for the two components. Significance and
# relevance are different questions, and a support decision that silently
# collapses them is the classic way to turn a large sample into a finding.
.nsca_components <- c("necessity", "sufficiency")

.nsca_validate_relevance <- function(relevance) {
  if (is.null(relevance)) {
    return(stats::setNames(c(NA_real_, NA_real_), .nsca_components))
  }

  # Names are a statement about which component gets which number, so a name
  # that is not one of the two is an error rather than a hint to fall back on
  # position.
  #
  # Before 0.5.0 the test was all(c("necessity", "sufficiency") %in% names(.)),
  # which required *both* names to be right and silently fell back to positional
  # order otherwise. c(sufficiency = 0.20, nec = 0.10) therefore gave necessity
  # 0.20 and sufficiency 0.10 -- exactly reversed, with no warning -- and
  # c(necessity = 0.10) had its single value copied to both sides. The argument
  # decides necessity_supported, sufficiency_supported and joint_support, so a
  # typo could flip a reported conclusion. The check now happens *before*
  # as.numeric() strips the names, which is what let the length-one case
  # through.
  labels <- names(relevance)
  if (!is.null(labels) && any(nzchar(labels))) {
    if (!setequal(labels, .nsca_components)) {
      stop(
        sprintf(
          paste0(
            "'relevance' was given the name(s) %s. When names are used, they ",
            "must be exactly 'necessity' and 'sufficiency' (in any order). ",
            "Without names, the two numbers are read by position."
          ),
          paste(sprintf("'%s'", labels), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    values <- suppressWarnings(as.numeric(relevance[.nsca_components]))
  } else {
    values <- suppressWarnings(as.numeric(relevance))
    if (length(values) == 1L) {
      values <- c(values, values)
    }
  }

  if (length(values) != 2L || anyNA(values) || any(values < 0)) {
    stop(
      paste0(
        "'relevance' must be NULL, one non-negative number applied to both ",
        "components, or two such numbers as c(necessity, sufficiency)."
      ),
      call. = FALSE
    )
  }
  stats::setNames(values, .nsca_components)
}

# Read one engine parameter without failing when it is absent.
#
# "Absent" means this condition or this frontier technique produced no value --
# a step frontier has no slope, `ols` never reaches the parameter matrix. That
# is data, and NA is the right answer.
#
# It does not mean the engine stopped reporting a quantity SCAtools asks for by
# name. That is an interface mismatch, and until 0.5.0 this tryCatch turned it
# into NA along with everything else, so a renamed engine parameter would have
# emptied a column in silence. SCAtools now raises a classed condition for that
# case and it is re-thrown here rather than absorbed.
.nsca_param <- function(fit, condition, ceiling, param, extractor) {
  value <- tryCatch(
    extractor(fit, x = condition, ceiling = ceiling, param = param),
    sca_engine_contract_error = function(e) stop(e),
    error = function(e) NA_real_
  )
  suppressWarnings(as.numeric(value)[[1L]])
}

# The necessity side reads NCA directly, so it needs the same contract check
# SCAtools applies to the sufficiency side.
.nsca_engine_params <- c(
  "Ceiling zone", "Effect size", "Slope", "Intercept",
  "p-value", "p-accuracy", "Fit", "Ceiling accuracy"
)

.nsca_require_engine_contract <- function(fit) {
  item <- fit$summaries[[1L]]
  if (is.null(item) || is.null(item$params) || ncol(item$params) == 0L) {
    return(invisible(fit))
  }
  missing <- c(
    setdiff(.nsca_engine_params, rownames(item$params)),
    setdiff("Number of observations", rownames(item$global))
  )
  if (length(missing) > 0L) {
    stop(
      sprintf(
        paste0(
          "The NCA engine (version %s) does not report '%s'. NSCA reads the ",
          "necessity results from the engine's output by name, so this ",
          "indicates an interface mismatch rather than a property of the data."
        ),
        as.character(utils::packageVersion("NCA")),
        paste(missing, collapse = "', '")
      ),
      call. = FALSE
    )
  }
  invisible(fit)
}

.nsca_nca_param <- function(fit, condition, ceiling, param) {
  .nsca_param(fit, condition, ceiling, param, NCA::nca_extract)
}

.nsca_sca_param <- function(fit, condition, ceiling, param) {
  .nsca_param(fit, condition, ceiling, param, SCAtools::sca_extract)
}

#' Necessary and sufficient condition analysis
#'
#' Estimates the two empty-space frontiers implied by a
#' necessary-and-sufficient statement and assesses them jointly. The
#' sufficiency side is estimated by [SCAtools::sca_analysis()] and the
#' necessity side by [NCA::nca_analysis()], on diagonally opposite corners of
#' the same scatter plot.
#'
#' @section Matched estimation:
#' Both sides are estimated with the same frontier technique, the same
#' theoretical scope, the same observations and the same outcome levels. Each
#' joint index compares the two empty zones, and such a comparison is
#' meaningful only when both sides are measured in the same way: envelopment
#' frontiers follow the data closely and yield the largest empty areas,
#' whereas regression frontiers cut into them, so mixing techniques across
#' sides would let the choice of technique, rather than the evidence,
#' determine which side appears weaker. Supplying several techniques produces
#' one internally matched row per technique, which offers a way to check
#' whether a conclusion is robust to that choice.
#'
#' @section Degenerate scopes:
#' A constant outcome is refused when no theoretical scope is given, because
#' it has no scope of its own and every empty area would be measured against
#' a scope of zero height. With a scope imposed, the analysis runs but issues
#' a warning and marks the affected rows `degenerate` in [nsca_table()]. An
#' axis that does not vary inside the declared scope makes every reported
#' area a property of the scope rather than of the data, and the resulting
#' areas overstate the evidence: a constant outcome at the center of the
#' scope leaves the upper and lower halves both empty, so both components
#' reach 0.5, the admissible region closes, and every joint index reaches its
#' maximum for data that carry no information. A milder warning is issued
#' when the observations span less than five percent of the declared scope on
#' either axis.
#'
#' @section The shared permutation sequence:
#' `shared.test.rep` replaces the engines' two separate permutation runs with
#' a single sequence used for both. Each replication shuffles the outcome
#' once, refits both engines on that same shuffled frame, and records
#' `d_nec`, `d_suf` and their minimum together. This is what makes
#' `p_weakest_perm` interpretable: the null distribution of a statistic
#' involving both components depends on how they co-vary under random
#' pairing, and independent permutations would remove that dependence. When
#' it is used, `p_nec` and `p_suf` are taken from the same sequence, the
#' `p_source` column reports `"shared"`, and all five p-values -- `p_nec`,
#' `p_suf`, `p_nsca_iut`, `p_weakest_perm` and `p_sum_perm` -- are mutually
#' consistent. It requires two engine fits per replication and is therefore
#' off by default.
#'
#' The orderings of the outcome are drawn following the rule NCA uses for its
#' own permutation test: when `shared.test.rep` is 720 or fewer, all orderings
#' are distinct, and a count larger than the `n!` orderings a sample of `n`
#' has is lowered to `n!`, so that each ordering is used once and the p-values
#' are exact.
#'
#' @param data A data frame or object coercible to a data frame.
#' @param x Columns containing one or more conditions.
#' @param y A single outcome column.
#' @param direction Necessary-and-sufficient direction(s): `"HH"`, `"LH"`,
#'   `"HL"`, or `"LL"`. The first letter is the condition level and the second
#'   the outcome level.
#' @param ceilings One or more empty-space frontier techniques, applied
#'   identically to both sides. `"ols"` is not accepted here because it
#'   estimates central tendency rather than an empty space; it can be supplied
#'   as `reference` instead.
#' @param reference Central-tendency lines drawn beside the two frontiers for
#'   comparison, or `NULL` for none. `"ols"` is the ordinary least-squares
#'   regression of the outcome on the condition, fitted to all the data. It is
#'   an average-effect summary; it does not enter any joint index, threshold
#'   or support decision, and is reported by [nsca_reference()].
#' @param scope Optional theoretical scope `c(xmin, xmax, ymin, ymax)`, shared
#'   by both sides.
#' @param threshold.x,threshold.y Reporting scales for [nsca_thresholds()], as
#'   in [SCAtools::sca_analysis()].
#' @param convention Reporting convention, as in [SCAtools::sca_analysis()].
#' @param steps Number of outcome levels, or an explicit vector of levels on
#'   the `threshold.y` scale.
#' @param step.size Optional spacing between outcome levels.
#' @param cutoff How out-of-range threshold values are represented.
#' @param reference.tau Quantile for a `"qr"` reference line, defaulting to
#'   `0.5` -- the median line. This is separate from `qr.tau`: a `qr`
#'   frontier lies at the edge of the cloud and uses `qr.tau`, whereas a `qr`
#'   reference line lies in its middle and uses `reference.tau`. Both can be
#'   requested at once, which the engine cannot express in a single call, so
#'   reference lines are computed by SCAtools rather than by the engine.
#' @param qr.tau Quantile used by the quantile-regression frontier, stated as
#'   the share of observations the line is intended to leave on the non-empty
#'   side. It is read that way in both directions: the engine fits `tau` when
#'   the empty space is an upper corner and `1 - tau` when it is a lower one
#'   (corners 3 and 4), so the default 0.95 is the 95th-percentile line where
#'   the empty space lies above and the 5th-percentile line where it lies
#'   below. Sufficiency in the `"HH"` direction has its empty space at high X
#'   with low Y, which is corner 4, and therefore fits the 5th percentile: a
#'   floor, not a ceiling. No manual inversion of `qr.tau` is needed to
#'   compensate.
#'
#'   Two reported quantities behave differently for `qr` than for the other
#'   frontiers, and both follow from the line being placed by a quantile
#'   rather than by the observations at the edge of the cloud.
#'   `frontier_accuracy` counts the share of observations on the non-empty
#'   side, which quantile regression sets to `qr.tau` by construction, so a
#'   value of 0.95 reflects the input rather than evidence of a good fit.
#'   `fit` compares this frontier's empty area with CE-FDH's, and the engine
#'   returns `NA` whenever this area is the larger of the two. This often
#'   happens for a quantile line, because it may pass through observations
#'   where an envelope cannot. An `NA` fit for `qr` is therefore an ordinary
#'   result rather than a failure.

#' @param test.rep Number of permutation resamples for the engines' own
#'   component tests. Zero skips testing, in which case no joint-support
#'   decision can be reached. The p-value cannot fall below `1 / test.rep`,
#'   nor below `1 / n!` for a sample of `n` observations, so a count or a
#'   sample too small to reach `test.p_threshold` is reported with a warning.
#' @param test.p_confidence Confidence level for permutation p-value accuracy.
#' @param test.p_threshold Significance level used by the intersection-union
#'   combination of the two directional component tests.
#' @param relevance Optional pre-specified practical-relevance thresholds for
#'   the component effect sizes, as one number applied to both or
#'   `c(necessity, sufficiency)`. Support requires an effect at least this
#'   large as well as a significant test. `NULL` means no relevance criterion
#'   was pre-specified, and support then rests on significance alone.
#' @param shared.test.rep Number of replications of the shared permutation
#'   sequence. Zero, the default, skips it; `p_weakest_perm` and `p_sum_perm`
#'   are then `NA` and the component p-values come from the engines. Its
#'   p-values cannot fall below `1 / (shared.test.rep + 1)`, and a count too
#'   small to reach `test.p_threshold` is reported with a warning.
#' @param test.seed Optional seed for the two component permutation tests, so
#'   that `p_nec`, `p_suf` and therefore `p_nsca_iut` and `joint_support` are
#'   reproducible. The caller's random stream is restored afterward.
#'
#'   A single seed yields two sequences. The two components are tested on two
#'   independent shuffles, and seeding both engines with the same number
#'   would give them the same sequence, that is, the shared-sequence design
#'   obtained unintentionally and without the bookkeeping that makes it
#'   interpretable. The seed given here is therefore used once, to draw two
#'   component seeds; `shared.test.rep` and `shared.seed` remain the way to
#'   request a single sequence explicitly.
#' @param shared.seed Optional seed for the shared sequence. The caller's
#'   random stream is restored afterward.
#' @param correspondence.tolerance How wide a necessity-sufficiency interval
#'   still counts as threshold correspondence, in the units of the condition.
#'   In a perfectly necessary-and-sufficient relation the two thresholds
#'   coincide at every outcome level, so the interval between them is zero;
#'   this argument sets how far from zero still counts as close enough.
#'
#'   The default of `0` is not a recommendation. The width at which a
#'   relation ceases to be necessary-and-sufficient is a substantive claim
#'   about the subject matter, and a package-supplied default would be the
#'   kind of unjustified tolerance the framework cautions against. We suggest
#'   setting it only when a particular width can be justified, and otherwise
#'   leaving it at zero and reading the interval itself.
#'
#'   It is distinct from the floating-point tolerance that [nsca_thresholds()]
#'   uses to decide whether two computed numbers are equal. That tolerance is
#'   a numerical matter; this one is substantive.
#' @param geometry.max_overlap Largest excess frontier overlap, beyond what a
#'   step frontier produces by construction, that still counts as acceptable
#'   geometry.
#' @param geometry.max_reconstruction Largest residual of the area identity
#'   that still counts as acceptable geometry.
#' @param na.action How missing values are dropped, as in
#'   [SCAtools::sca_analysis()]. Both sides use the same rule, because a
#'   necessity effect size and a sufficiency effect size are comparable only
#'   when both are areas over the same box.
#' @param purity Whether to compute the engine's additional purity metrics
#'   where it offers them. Off by default: the engine computes them only for a
#'   corner with neither axis flipped, so at most one side of an NSCA model
#'   can have them, and which side, if any, depends on the direction. No NSCA
#'   statistic uses them, and they are computationally expensive on data with
#'   many frontier points.
#' @return An object of class `nsca_result`.
#' @references
#' Chan, Y. (in press). *From Necessity to Sufficiency and
#' Necessity-and-Sufficiency: Extending Condition Analysis in Degree*.
#' University of Leeds.
#'
#' Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology
#' of "necessary but not sufficient" causality. *Organizational Research
#' Methods*, 19(1), 10-52. \doi{10.1177/1094428115584005}
#' @seealso [nsca_table()], [nsca_results()], [nsca_joint()],
#'   [nsca_thresholds()], [nsca_corners()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", direction = "HH", ceilings = "ce_fdh")
#' fit
nsca_analysis <- function(
    data,
    x,
    y,
    direction = "HH",
    ceilings = c("ce_fdh", "cr_fdh"),
    reference = NULL,
    scope = NULL,
    threshold.x = "percentage.range",
    threshold.y = "percentage.range",
    convention = c("absolute", "directional"),
    steps = 10,
    step.size = NULL,
    cutoff = 0,
    qr.tau = 0.95,
    reference.tau = 0.5,
    test.rep = 0,
    test.p_confidence = 0.95,
    test.p_threshold = 0.05,
    relevance = NULL,
    shared.test.rep = 0,
    test.seed = NULL,
    shared.seed = NULL,
    correspondence.tolerance = 0,
    geometry.max_overlap = 0.01,
    geometry.max_reconstruction = 0.02,
    na.action = c("pairwise", "listwise"),
    purity = FALSE) {

  test.rep <- .nsca_validate_test_rep(test.rep, "test.rep")
  shared.test.rep <- .nsca_validate_test_rep(shared.test.rep,
                                             "shared.test.rep")
  test.p_threshold <- .nsca_validate_alpha(test.p_threshold,
                                           "test.p_threshold")
  test.p_confidence <- .nsca_validate_alpha(test.p_confidence,
                                            "test.p_confidence")

  frame <- .nsca_validate_frame(data)
  x_names <- .nsca_select_columns(frame, x, "x")
  y_name <- .nsca_select_columns(frame, y, "y")
  if (length(y_name) != 1L) {
    stop("'y' must select exactly one outcome column.", call. = FALSE)
  }
  if (length(x_names) == 0L) {
    stop("At least one condition must be selected in 'x'.", call. = FALSE)
  }
  # A variable compared with itself puts both frontiers on the same
  # diagonal, so d_nec and d_suf are fixed by the geometry rather than by
  # the data and the admissible region collapses to nothing.
  if (y_name %in% x_names) {
    stop(
      sprintf(
        paste0(
          "'%s' is named as both a condition and the outcome. Both ",
          "frontiers would lie on the same diagonal, so the two effect ",
          "sizes would reflect the geometry rather than the data."
        ),
        y_name
      ),
      call. = FALSE
    )
  }

  direction <- .nsca_validate_direction(direction, length(x_names))
  y_letters <- unique(substr(direction, 2L, 2L))
  if (length(y_letters) > 1L) {
    stop(
      paste0(
        "All conditions in one analysis must use the same outcome direction. ",
        "High-Y and low-Y statements require separate models."
      ),
      call. = FALSE
    )
  }
  ceilings <- .nsca_validate_ceilings(ceilings)
  reference <- .nsca_validate_reference(reference)
  convention <- match.arg(convention)
  relevance <- .nsca_validate_relevance(relevance)
  correspondence.tolerance <- .nsca_validate_correspondence(
    correspondence.tolerance, x_names
  )
  geometry.max_overlap <- .nsca_validate_tolerance(
    geometry.max_overlap, "geometry.max_overlap")
  geometry.max_reconstruction <- .nsca_validate_tolerance(
    geometry.max_reconstruction, "geometry.max_reconstruction")
  if (!is.numeric(shared.test.rep) || length(shared.test.rep) != 1L ||
      is.na(shared.test.rep) || shared.test.rep < 0) {
    stop("'shared.test.rep' must be a single non-negative number.",
         call. = FALSE)
  }
  shared.test.rep <- as.integer(shared.test.rep)

  necessity_corner <- unname(.nsca_necessity_corner[direction])
  sufficiency_corner <- unname(.nsca_sufficiency_corner[direction])

  na.action <- .nsca_validate_na_action(na.action)
  # Before anything reads the scope: one outcome column may not be given two
  # theoretical ranges. See .nsca_validate_scope_y().
  .nsca_validate_scope_length(scope, length(x_names))
  .nsca_validate_scope_bounds(scope, length(x_names))
  .nsca_validate_scope_y(scope, x_names, y_name)

  raw <- .nsca_raw_data(frame, x_names, y_name, scope, na.action)

  # Both engines do their own per-condition complete-case filtering, so
  # `listwise` has to be applied here or it changes this package's bookkeeping
  # and nothing either engine reports -- including putting d_nec and d_suf on
  # different samples.
  frame <- .nsca_apply_na_action(frame, x_names, y_name, na.action)
  scope_adjustment <- .nsca_scope_adjustment(scope, raw, y_name)

  # Before either engine runs: a replication count, or a sample, too small for
  # any p-value to reach the threshold makes every test a non-rejection by
  # arithmetic. The engines cap their counts at n! of the rows they are handed;
  # the shared sequence at n! of the complete rows it shuffles.
  .nsca_warn_test_floor(
    test.rep, shared.test.rep, test.p_threshold,
    observations = nrow(frame),
    shared_observations = sum(stats::complete.cases(frame[c(x_names, y_name)]))
  )

  # SCAtools clears the flag its own engine call leaves behind, but the
  # NCA::nca_analysis() call below creates it again in the caller's global
  # environment. Clear it here unless it was already there.
  flag_existed <- exists("scope_warnings_nca", envir = globalenv(),
                         inherits = FALSE)
  on.exit(.nsca_clear_engine_scope_flag(flag_existed), add = TRUE)
  .nsca_require_varying_condition(raw)

  if (diff(raw$y_bounds) <= 0) {
    stop(
      paste0(
        "The outcome takes a single value, so it has no scope of its own and ",
        "neither component is defined: every empty area would be measured ",
        "against a scope of zero height. If the outcome is bounded by design, ",
        "a theoretical 'scope' with distinct ymin and ymax can be supplied; ",
        "the result then depends on where the constant lies within that scope."
      ),
      call. = FALSE
    )
  }

  # The engine computes its purity metrics only when neither axis is flipped.
  # Of the two corners an NSCA model uses, at most one is unflipped, so these
  # metrics can never be available symmetrically, and which side gets them is
  # decided by the direction rather than by anything substantive. They also
  # cost a great deal on data with many frontier points. Off unless asked for.
  if (!isTRUE(purity)) {
    previous <- Sys.getenv("NCA_SKIP_PURITY", unset = NA)
    Sys.setenv(NCA_SKIP_PURITY = "TRUE")
    on.exit(
      if (is.na(previous)) {
        Sys.unsetenv("NCA_SKIP_PURITY")
      } else {
        Sys.setenv(NCA_SKIP_PURITY = previous)
      },
      add = TRUE
    )
  }

  # One seed in, two out. Seeding both engines with the same number would hand
  # them the same permutation sequence: each derives its whole sequence from a
  # single draw off the caller's stream, so identical starting states give
  # identical shuffles. The two component p-values would then be computed on
  # one sequence while p_source still reported two independent engine runs --
  # the shared-sequence design with none of its guarantees and none of its
  # labeling. Drawing two seeds keeps the components independent and keeps
  # both reproducible.
  component_seeds <- .nsca_component_seeds(test.seed)

  # Sufficiency first: it settles the outcome levels, which the necessity side
  # then reuses so that the two threshold tables line up row for row.
  # SCAtools warns on its own when its test cannot reach the threshold; the
  # warning above already covers both components, so that one is muffled.
  sufficiency <- withCallingHandlers(SCAtools::sca_analysis(
    data = frame,
    x = x,
    y = y,
    direction = direction,
    ceilings = ceilings,
    reference = reference,
    scope = scope,
    threshold.x = threshold.x,
    threshold.y = threshold.y,
    convention = convention,
    na.action = na.action,
    steps = steps,
    step.size = step.size,
    cutoff = cutoff,
    qr.tau = qr.tau,
    reference.tau = reference.tau,
    test.rep = test.rep,
    test.p_confidence = test.p_confidence,
    test.p_threshold = test.p_threshold,
    test.seed = component_seeds$sufficiency
  ), sca_test_floor = function(w) invokeRestart("muffleWarning"))
  sca_meta <- attr(sufficiency, "sca")
  outcome_levels <- sca_meta$outcome_levels

  # The engine calls set.seed(NULL) when it finishes, which reseeds from the
  # clock rather than restoring anything, so the seed has to be applied and
  # the caller's stream put back around the call. SCAtools does the same for
  # its own engine call through .sca_with_seed().
  necessity <- .nsca_with_seed(component_seeds$necessity, NCA::nca_analysis(
    data = frame,
    x = x,
    y = y,
    ceilings = ceilings,
    corner = necessity_corner,
    scope = scope,
    bottleneck.x = "actual",
    bottleneck.y = "actual",
    steps = outcome_levels,
    step.size = NULL,
    cutoff = cutoff,
    qr.tau = qr.tau,
    effect_aggregation = 1,
    test.rep = test.rep,
    test.p_confidence = test.p_confidence,
    test.p_threshold = test.p_threshold
  ))

  .nsca_require_engine_contract(necessity)

  condition_names <- names(sufficiency$summaries)
  outcome_name <- sca_meta$outcome
  names(direction) <- condition_names
  names(necessity_corner) <- condition_names
  names(sufficiency_corner) <- condition_names

  raw$conditions <- raw$conditions[match(condition_names, names(raw$conditions))]
  names(raw$conditions) <- condition_names

  # Every reported quantity is an area over the declared scope. When one axis
  # barely moves inside that scope, the areas are decided by the scope and not
  # by the data, and they are decided in the flattering direction. Flag it
  # here, once, rather than leaving it to be read off a table.
  degeneracy <- lapply(condition_names, function(condition) {
    item <- raw$conditions[[condition]]
    x_coverage <- .nsca_axis_coverage(item$x, item$x_bounds)
    y_coverage <- .nsca_axis_coverage(item$y, item$y_bounds)
    worst <- min(c(x_coverage, y_coverage))
    list(
      x_coverage = x_coverage,
      y_coverage = y_coverage,
      degenerate = !is.finite(worst) || worst <= .nsca_degenerate_tolerance,
      thin = is.finite(worst) && worst > .nsca_degenerate_tolerance &&
        worst < .nsca_thin_coverage
    )
  })
  names(degeneracy) <- condition_names

  flagged <- function(field) {
    condition_names[vapply(degeneracy, function(item) isTRUE(item[[field]]),
                           logical(1L))]
  }
  flat <- flagged("degenerate")
  if (length(flat) > 0L) {
    warning(
      sprintf(
        paste0(
          "%s: the condition or the outcome does not vary inside the declared ",
          "scope. Every empty area is then determined by where that constant ",
          "lies in the scope rather than by the data, and the resulting areas ",
          "overstate the evidence: a constant outcome at the center of an ",
          "imposed scope gives d_nec = d_suf = 0.5, closes the admissible ",
          "region, and scores 1 on every joint index although it carries no ",
          "information. The permutation p-values are more informative here; ",
          "nsca_table() marks these rows as degenerate, with ",
          "geometry_acceptable set to FALSE."
        ),
        paste(flat, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  thin <- flagged("thin")
  if (length(thin) > 0L) {
    warning(
      sprintf(
        paste0(
          "%s: the observations span less than %g%% of the declared scope on ",
          "at least one axis, so most of each reported area consists of ",
          "unobserved scope rather than evidence. It may be worth confirming ",
          "that this is the intended scope."
        ),
        paste(thin, collapse = ", "), 100 * .nsca_thin_coverage
      ),
      call. = FALSE
    )
  }

  statements <- vapply(
    condition_names,
    function(condition) {
      .nsca_statement(direction[[condition]], condition, outcome_name)
    },
    character(1L)
  )

  # One permutation sequence for both engines, when asked for.
  shared <- NULL
  if (shared.test.rep > 0L) {
    complete <- frame[stats::complete.cases(frame[c(x_names, y_name)]), ,
                      drop = FALSE]
    fitted_counts <- vapply(raw$conditions, function(item) length(item$x),
                            integer(1L))
    if (nrow(complete) < max(fitted_counts)) {
      warning(
        paste0(
          "The shared permutation sequence uses listwise-complete rows (",
          nrow(complete), "), while at least one component was fitted on more ",
          "(", max(fitted_counts), "). Its p-values therefore describe a ",
          "smaller sample than the observed effect sizes."
        ),
        call. = FALSE
      )
    }
    # Address the conditions by the engines' own names throughout, so that the
    # direction and the corner cannot drift out of step with the column order
    # if an engine ever returns its summaries in a different order than it was
    # given them.
    shared <- .nsca_shared_null(
      frame = complete, x_names = condition_names, y_name = y_name,
      direction = unname(direction[condition_names]), ceilings = ceilings,
      scope = scope,
      necessity_corner = unname(necessity_corner[condition_names]),
      cutoff = cutoff, qr.tau = qr.tau,
      reps = shared.test.rep, seed = shared.seed
    )
  }

  model <- list(necessity = necessity, sufficiency = sufficiency)
  attr(model, "nsca") <- list(
    call = match.call(),
    direction = direction,
    necessity_corner = necessity_corner,
    sufficiency_corner = sufficiency_corner,
    statement = statements,
    conditions = condition_names,
    outcome = outcome_name,
    outcome_letter = y_letters[[1L]],
    ceilings = ceilings,
    reference = reference,
    reference.tau = reference.tau,
    outcome_levels = outcome_levels,
    threshold.x = threshold.x,
    threshold.y = threshold.y,
    convention = convention,
    cutoff = cutoff,
    purity = isTRUE(purity),
    degeneracy = degeneracy,
    relevance = relevance,
    na.action = na.action,
    # Both halves of the scope decision, so that print() can say which one was
    # in force. Every d in this package is an area over the scope area.
    scope_declared = scope,
    scope_source = if (is.null(scope)) "observed" else "declared",
    # One row per bound the engines had to widen to hold the data; empty
    # when the declared scope already contained every observation.
    scope_adjustment = scope_adjustment,
    shared.seed = shared.seed,
    correspondence.tolerance = correspondence.tolerance,
    geometry = c(max_overlap = geometry.max_overlap,
                 max_reconstruction = geometry.max_reconstruction),
    test.rep = test.rep,
    test.seed = test.seed,
    # The rows handed to the engines, which cap their permutation counts at
    # n! of them and so set the smallest p-value a test can return.
    observations = nrow(frame),
    # The count aimed at: the request, or n! when that is smaller.
    shared.test.rep = if (is.null(shared)) {
      shared.test.rep
    } else {
      attr(shared, "target")
    },
    shared = shared,
    alpha = test.p_threshold,
    p_confidence = test.p_confidence,
    raw = raw,
    engines = c(
      NCA = as.character(utils::packageVersion("NCA")),
      SCAtools = as.character(utils::packageVersion("SCAtools"))
    )
  )
  class(model) <- "nsca_result"
  model
}

#' Quick necessary and sufficient condition analysis
#'
#' A convenience wrapper around [nsca_analysis()] that also runs the
#' permutation tests.
#'
#' Because the tests run by default, we suggest supplying `test.seed` so that
#' the p-values are reproducible. An unseeded permutation test draws its seed
#' from the clock and reports a different p-value on each run, which makes
#' the value unsuitable for reporting; the `Settings` block states whether a
#' seed was fixed and, if it was not, adds a note that the p-values may
#' differ on the next run.
#'
#' @inheritParams nsca_analysis
#' @return An object of class `nsca_result`.
#' @seealso [nsca_analysis()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' \donttest{
#' # The default 1000 permutations take longer than an example should;
#' # nsca_analysis() runs none unless asked.
#' nsca(dat, "X", "Y")
#' }
nsca <- function(data, x, y, direction = "HH",
                 ceilings = c("ce_fdh", "cr_fdh"),
                 reference = NULL, test.rep = 1000, test.seed = NULL) {
  model <- nsca_analysis(
    data = data, x = x, y = y, direction = direction,
    ceilings = ceilings, reference = reference, test.rep = test.rep,
    test.seed = test.seed
  )
  # nsca() forwards all of its arguments to nsca_analysis() by name, so that
  # function's match.call() cannot tell a caller's choice from a default.
  # Record the call the researcher actually wrote; print() reads it to mark
  # the settings block.
  attr(model, "nsca")$call <- match.call()
  attr(model, "show.plots") <- TRUE
  model
}

# Same two rules SCAtools offers, validated here rather than reached for with
# ::: so that NSCA does not depend on another package's internals.
.nsca_validate_na_action <- function(na.action) {
  if (is.null(na.action) || !is.character(na.action) || length(na.action) == 0L) {
    stop("'na.action' must be \"pairwise\" or \"listwise\".", call. = FALSE)
  }
  na.action <- tolower(trimws(na.action[[1L]]))
  if (!(na.action %in% c("pairwise", "listwise"))) {
    stop(
      sprintf("Invalid 'na.action': '%s'. Expected \"pairwise\" or \"listwise\".",
              na.action),
      call. = FALSE
    )
  }
  na.action
}

# The geometry tolerances, validated like their neighbors.
#
# `relevance` and `shared.test.rep` have had validators since 0.4.0; these two
# had none, and they are stored with c(max_overlap = , max_reconstruction = ).
# A character value therefore coerced the whole vector to character, and the
# later comparison `excess_overlap <= "0.01"` became a string comparison that
# silently produced the wrong geometry_acceptable. NA propagated instead,
# turning joint_support into NA for reasons nothing reported.
.nsca_validate_tolerance <- function(value, name) {
  values <- suppressWarnings(as.numeric(value))
  if (length(values) != 1L || anyNA(values) || values < 0) {
    stop(
      sprintf("'%s' must be a single non-negative number.", name),
      call. = FALSE
    )
  }
  values
}

# A correspondence tolerance is a width in the units of one condition.
#
# `correspondence.tolerance` says how wide a necessity-sufficiency interval
# still counts as correspondence, measured in the units of X. With one
# condition that is one number. With several on different measurement scales
# -- budget in millions beside a seven-point scale, which is ordinary in this
# field -- one number cannot be right for both: "2" is a rounding error on one
# axis and a third of the range on the other.
#
# So a single non-zero value is accepted only when there is one condition, and
# refused otherwise with the names to use. Zero is exempt because zero is zero
# in every unit, and it is the default: the refusal is about an unstated scale,
# not about arithmetic.
#
# This is the same rule as the per-condition scope (D-19) and for the same
# reason: a quantity whose meaning depends on the axis may not be stated once
# for axes that differ.
.nsca_validate_correspondence <- function(value, x_names) {
  if (is.null(value)) {
    value <- 0
  }
  values <- suppressWarnings(as.numeric(value))
  if (length(values) == 0L || anyNA(values) || any(values < 0)) {
    stop(
      paste0(
        "'correspondence.tolerance' must be a non-negative number, or a ",
        "vector of such numbers named after the conditions."
      ),
      call. = FALSE
    )
  }

  labels <- names(value)
  if (!is.null(labels) && any(nzchar(labels))) {
    missing <- setdiff(x_names, labels)
    unknown <- setdiff(labels, x_names)
    if (length(missing) > 0L || length(unknown) > 0L) {
      stop(
        sprintf(
          paste0(
            "'correspondence.tolerance' must name every condition and no ",
            "others. Missing: %s. Unknown: %s."
          ),
          if (length(missing) == 0L) "none" else paste(missing, collapse = ", "),
          if (length(unknown) == 0L) "none" else paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(stats::setNames(unname(values[match(x_names, labels)]), x_names))
  }

  if (length(values) == 1L) {
    if (length(x_names) > 1L && values[[1L]] > 0) {
      stop(
        sprintf(
          paste0(
            "'correspondence.tolerance' is a width in the units of the ",
            "condition, and this analysis has %d conditions (%s), which need ",
            "not share a scale, so a single number may not suit all of them.\n",
            "\n",
            "One value per condition can be supplied by name:\n",
            "  correspondence.tolerance = c(%s)\n",
            "\n",
            "Alternatively, it can be left at 0 and the interval read directly: ",
            "zero is zero in every unit, which is why it is the default and is ",
            "exempt from this check."
          ),
          length(x_names), paste(x_names, collapse = ", "),
          paste(sprintf("%s = ", x_names), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(stats::setNames(rep(values[[1L]], length(x_names)), x_names))
  }

  if (length(values) != length(x_names)) {
    stop(
      sprintf(
        paste0(
          "'correspondence.tolerance' has %d values for %d conditions. ",
          "Expected one number, or one per condition named after the ",
          "conditions."
        ),
        length(values), length(x_names)
      ),
      call. = FALSE
    )
  }
  stats::setNames(values, x_names)
}

# The frontier techniques the engines understand. The counterpart of
# SCAtools' `.sca_known_ceilings`, duplicated for the reason every other
# validator here is: that one is private, and NSCA reaches NCA directly for
# the necessity side. A test pins both lists against the engine's own.
.nsca_known_ceilings <- c(
  "ce_fdh", "ce_vrs", "cr_fdh", "cr_vrs", "cols", "qr", "c_lp", "ols"
)

.nsca_validate_alpha <- function(value, label) {
  v <- suppressWarnings(as.numeric(value))
  if (length(v) != 1L || !is.finite(v) || v <= 0 || v >= 1) {
    stop(sprintf("'%s' must be a single number strictly between 0 and 1.",
                 label), call. = FALSE)
  }
  v
}

.nsca_validate_test_rep <- function(value, label) {
  v <- suppressWarnings(as.numeric(value))
  if (length(v) != 1L || !is.finite(v) || v < 0 || v != round(v)) {
    stop(sprintf("'%s' must be a single non-negative whole number.", label),
         call. = FALSE)
  }
  as.integer(v)
}

# Column selection that names the argument. `frame[x]` raises R's own
# "undefined columns selected", which mentions neither x nor the name that
# was not found.
.nsca_select_columns <- function(frame, selector, label) {
  if (is.numeric(selector)) {
    bad <- selector[selector < 1 | selector > ncol(frame)]
    if (length(bad) > 0L) {
      stop(sprintf("'%s' selects column %s, but the data has %d columns.",
                   label, paste(bad, collapse = ", "), ncol(frame)),
           call. = FALSE)
    }
    return(colnames(frame)[selector])
  }
  selector <- as.character(selector)
  missing <- setdiff(selector, colnames(frame))
  if (length(missing) > 0L) {
    stop(
      sprintf("'%s' names %s, which is not a column of the data. Available: %s.",
              label, paste(sprintf("'%s'", missing), collapse = ", "),
              paste(colnames(frame), collapse = ", ")),
      call. = FALSE
    )
  }
  duplicated_names <- unique(selector[duplicated(selector)])
  if (length(duplicated_names) > 0L) {
    stop(
      sprintf("'%s' names %s more than once.", label,
              paste(sprintf("'%s'", duplicated_names), collapse = ", ")),
      call. = FALSE
    )
  }
  selector
}

.nsca_validate_frame <- function(data) {
  if (is.null(data) || (!is.data.frame(data) && !is.matrix(data))) {
    stop("'data' must be a data frame or a matrix.", call. = FALSE)
  }
  frame <- as.data.frame(data)
  if (nrow(frame) == 0L || ncol(frame) == 0L) {
    stop("'data' has no rows or no columns.", call. = FALSE)
  }
  frame
}
