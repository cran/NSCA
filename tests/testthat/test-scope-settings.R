# The scope is the denominator both components share, so it decides whether
# d_nec and d_suf are on one scale at all. Everything that reports it, warns
# about it or refuses it is checked here, together with the Settings block.

sc_data <- function(n = 60, seed = 4) {
  set.seed(seed)
  x <- sort(runif(n))
  data.frame(X = x, Y = pmin(pmax(x + rnorm(n, 0, 0.08), 0), 1))
}

sc_two <- function(n = 60, seed = 4) {
  dat <- sc_data(n, seed)
  data.frame(A = dat$X, B = dat$X * 2, Y = dat$Y)
}

sc_lines <- function(fit) {
  strsplit(
    paste(utils::capture.output(suppressWarnings(print(fit))), collapse = "\n"),
    "\n", fixed = TRUE
  )[[1L]]
}

sc_line <- function(lines, label) {
  hit <- grep(paste0("^  ", label), lines, value = TRUE)
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}


# ---- nsca_scope() --------------------------------------------------------

test_that("nsca_scope() reports both scopes, one row per condition", {
  fit <- nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                       scope = c(0, 2, 0, 1))
  s <- nsca_scope(fit)

  expect_equal(nrow(s), 2L)
  expect_setequal(s$condition, c("A", "B"))
  expect_identical(unique(s$source), "declared")
  expect_false(any(s$adjusted))
  expect_true(all(is.finite(s$scope_area)))
})

test_that("nsca_scope() and sca_scope() agree on their column layout", {
  # The two packages are read side by side; a reader comparing a scope here
  # with one there should not have to translate.
  dat <- sc_data()
  nsca_cols <- names(nsca_scope(
    nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
  ))
  sca_cols <- names(SCAtools::sca_scope(
    SCAtools::sca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
  ))
  expect_identical(nsca_cols, sca_cols)
})

test_that("the joint table carries the shared denominator", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh",
                       scope = c(0, 1, 0, 1))
  tb <- nsca_table(fit)

  expect_true(all(c("scope_area", "scope_adjusted") %in% names(tb)))
  expect_equal(tb$scope_area[[1L]], 1)
  expect_false(tb$scope_adjusted[[1L]])
})

test_that("the component tables keep the same layout as each other", {
  # scope_area belongs to the joint row, not to the two component layers.
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  tables <- nsca_results(fit)
  expect_identical(names(tables$necessity), names(tables$sufficiency))
})


# ---- the engine's flag, and the widening -------------------------------

test_that("the engine's scope flag does not survive an NSCA call", {
  # NSCA reaches the engine twice: SCAtools clears the flag its own call
  # leaves, and the direct NCA call creates it again.
  if (exists("scope_warnings_nca", envir = globalenv(), inherits = FALSE)) {
    rm("scope_warnings_nca", envir = globalenv())
  }
  suppressWarnings(
    nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh",
                  scope = c(0.2, 0.8, 0.2, 0.8))
  )
  expect_false(exists("scope_warnings_nca", envir = globalenv(),
                      inherits = FALSE))
})

test_that("a widened scope is recorded and printed", {
  fit <- suppressWarnings(
    nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh",
                  scope = c(0.2, 0.8, 0.2, 0.8))
  )
  expect_true(nsca_scope(fit)$adjusted)
  expect_match(paste(sc_lines(fit), collapse = "\n"), "widened by the engines")
})


# ---- one outcome column, one theoretical range (D-19) -------------------

test_that("a per-condition scope with two outcome ranges is refused", {
  expect_error(
    nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  scope = list(c(0, 1, 0, 1), c(0, 2, 0, 10))),
    "more than one theoretical range"
  )
})

test_that("the refusal explains what would break and what would not", {
  err <- tryCatch(
    nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  scope = c(0, 1, 0, 1, 0, 2, 0, 10)),
    error = function(e) conditionMessage(e)
  )
  # The identity is the NSCA-specific consequence.
  expect_match(err, "admissible_region_share", fixed = TRUE)
  # And the two components are not the problem, which a reader will assume.
  expect_match(err, "two corners of one scatter plot", fixed = TRUE)
})

test_that("differing X ranges with one outcome range still work", {
  fit <- nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                       scope = list(c(0, 1, 0, 1), c(0, 2, 0, 1)))
  expect_s3_class(fit, "nsca_result")
  expect_equal(nrow(nsca_table(fit)), 2L)
})


# ---- the Settings block --------------------------------------------------

test_that("print() ends with a Settings block naming every choice", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  lines <- sc_lines(fit)

  expect_true(any(grepl("^Settings$", lines)))
  for (label in c("frontiers", "test.rep", "test.seed", "shared.test.rep",
                  "shared.seed", "alpha", "relevance", "scope", "na.action",
                  "thresholds", "convention", "geometry")) {
    expect_false(is.na(sc_line(lines, label)), info = label)
  }
})

test_that("the geometry tolerances are named, because they decide support", {
  # geometry_acceptable gates joint_support, and a verdict that rests on a
  # number the reader has never seen is not a reportable verdict.
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  geometry <- sc_line(sc_lines(fit), "geometry")
  expect_match(geometry, "max_overlap")
  expect_match(geometry, "max_reconstruction")
})

test_that("an argument the call named is not marked as a default", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 20, test.seed = 1)
  lines <- sc_lines(fit)
  expect_false(grepl("(default)", sc_line(lines, "test.rep"), fixed = TRUE))
  expect_true(grepl("(default)", sc_line(lines, "na.action"), fixed = TRUE))
})

test_that("nsca() records its own call so its defaults read as defaults", {
  fit <- suppressWarnings(nsca(sc_data(), "X", "Y", test.rep = 0))
  expect_identical(as.character(attr(fit, "nsca")$call[[1L]]), "nsca")
  expect_true(
    grepl("(default)", sc_line(sc_lines(fit), "frontiers"), fixed = TRUE)
  )
})

test_that("nsca_output() ends with the same block", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  printed <- paste(
    utils::capture.output(suppressWarnings(
      nsca_output(fit, plots = FALSE)
    )),
    collapse = "\n"
  )
  expect_match(printed, "Settings")
  expect_match(printed, "max_reconstruction")
})


# ---- correspondence.tolerance (F8) --------------------------------------

test_that("at the default tolerance the verdict columns are absent", {
  # At zero they would carry no judgment, and a binary verdict beside a
  # continuous interval is the column that gets read. The default output
  # reports the interval.
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  th <- nsca_thresholds(fit)

  expect_false("corresponds" %in% names(th))
  expect_false("correspondence_tolerance" %in% names(th))
  # The interval itself is still there, which is the point.
  expect_true("necessity_sufficiency_interval_actual" %in% names(th))
  expect_false(any(th$threshold_status == "near correspondence"))
})

test_that("declaring a tolerance brings the verdict columns", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  th <- nsca_thresholds(fit, correspondence.tolerance = 0.02)
  expect_true(all(c("corresponds", "correspondence_tolerance") %in% names(th)))
  expect_true(all(th$correspondence_tolerance >= 0.02))
})

test_that("a small crossing under a line technique reads as correspondence", {
  # An exact envelope cannot cross itself. A fitted line can, and the narrower
  # the band the more easily two separately fitted lines cross inside it -- so
  # a small crossing under cr_fdh is the signature of a band thinner than the
  # fitting error, which is what an almost-perfect relation looks like.
  # Treating it as a fault would discard exactly the strongest case.
  set.seed(11)
  n <- 80
  x <- sort(runif(n))
  # A very thin band: near-deterministic, so the two frontiers nearly coincide.
  dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(n, 0, 0.01), 0), 1))
  fit <- nsca_analysis(dat, "X", "Y", ceilings = "cr_fdh",
                       scope = c(0, 1, 0, 1))
  th <- nsca_thresholds(fit, ceiling = "cr_fdh",
                        correspondence.tolerance = 0.05)
  crossed <- !is.na(th$threshold_gap_actual) & th$threshold_gap_actual < 0 &
    abs(th$threshold_gap_actual) <= 0.05
  if (any(crossed)) {
    expect_true(all(th$corresponds[crossed]))
    expect_true(all(th$threshold_status[crossed] %in%
                      c("exact correspondence", "near correspondence")))
  }
  # And the rule holds by construction whatever the data did.
  usable <- is.finite(th$threshold_gap_actual)
  expect_equal(
    th$corresponds[usable],
    ifelse(abs(th$threshold_gap_actual[usable]) <= 0.05, TRUE,
           ifelse(th$threshold_gap_actual[usable] > 0, FALSE, NA))
  )
})

test_that("a crossing too large to be fitting error has no answer", {
  # Beyond the declared tolerance in the negative direction the floor really
  # is above the ceiling over a stretch, which no relation produces and no
  # tolerance excuses. `width_actual` is pmax(gap, 0) and would read that as a
  # perfect interval, so the column is built from the signed gap instead.
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  th <- nsca_thresholds(fit, correspondence.tolerance = 1e-9)
  overlapping <- th$threshold_status == "frontiers overlap"
  if (any(overlapping)) {
    expect_true(all(is.na(th$corresponds[overlapping])))
  }
  # An exact envelope should essentially never cross: CE-FDH's two frontiers
  # are a running minimum and a running maximum over the same points.
  expect_true(sum(overlapping, na.rm = TRUE) <= 1L)
})

test_that("a declared tolerance opens the near-correspondence band", {
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh")
  wide <- nsca_thresholds(fit, correspondence.tolerance = 1)
  expect_true("corresponds" %in% names(wide))
  expect_true(all(wide$correspondence_tolerance >= 1))
  # Everything estimable and compatible is within a tolerance of one whole
  # unit on a condition that spans [0, 1].
  usable <- !is.na(wide$corresponds)
  expect_true(all(wide$corresponds[usable]))
})

test_that("the tolerance can be overridden without refitting", {
  # Reporting the verdict at two or three tolerances is how a reader sees how
  # sensitive it is, so the override must not need a new fit.
  fit <- nsca_analysis(sc_data(), "X", "Y", ceilings = "ce_fdh",
                       correspondence.tolerance = 0)
  strict <- nsca_thresholds(fit)
  loose <- nsca_thresholds(fit, correspondence.tolerance = 1)
  expect_false("correspondence_tolerance" %in% names(strict))
  expect_true(all(loose$correspondence_tolerance >= 1))
})

test_that("one tolerance for several conditions on different scales is refused", {
  # A width in the units of X cannot be stated once for axes that need not
  # share a scale. Same rule as the per-condition scope, same reason.
  expect_error(
    nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  correspondence.tolerance = 0.05),
    "units of the condition"
  )
})

test_that("zero is exempt, because zero is zero in every unit", {
  fit <- nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                       correspondence.tolerance = 0)
  expect_s3_class(fit, "nsca_result")
})

test_that("a named tolerance, one per condition, is accepted", {
  fit <- nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                       correspondence.tolerance = c(A = 0.01, B = 0.02))
  stored <- attr(fit, "nsca")$correspondence.tolerance
  expect_equal(unname(stored[["A"]]), 0.01)
  expect_equal(unname(stored[["B"]]), 0.02)
})

test_that("a tolerance that names the wrong conditions is refused", {
  expect_error(
    nsca_analysis(sc_two(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  correspondence.tolerance = c(A = 0.01, C = 0.02)),
    "Unknown"
  )
})


# ---- a flat condition, the mirror of the flat outcome refusal -----------

test_that("an analysis in which no condition varies is refused", {
  flat <- data.frame(X = rep(0.5, 30), Y = runif(30))
  expect_error(nsca_analysis(flat, "X", "Y", ceilings = "ce_fdh"),
               "No condition varies")
})

test_that("the two degeneracies are refused the same way", {
  # A constant outcome was refused; a constant condition was not, and returned
  # NaN. Same geometry on the other axis, so it gets the same treatment.
  flat_x <- data.frame(X = rep(0.5, 30), Y = runif(30))
  flat_y <- data.frame(X = runif(30), Y = rep(0.5, 30))
  expect_error(nsca_analysis(flat_x, "X", "Y", ceilings = "ce_fdh"))
  expect_error(nsca_analysis(flat_y, "X", "Y", ceilings = "ce_fdh"))
})
