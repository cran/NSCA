# nsca_extract(): the single-value accessor.
#
# Two things matter here. It must never invent a number, and it must be honest
# about which engine a value came from. The nec:/suf: prefixes exist so that a
# reader can tell a joint quantity from a component one at the call site, and
# these tests hold that distinction in place.

extract_data <- function(n = 45, seed = 31) {
  set.seed(seed)
  x <- sort(runif(n))
  data.frame(
    X1 = x,
    X2 = runif(n),
    Y = pmin(pmax(x + rnorm(n, 0, 0.10), 0), 1)
  )
}

test_that("columns of the joint table are reachable by name", {
  fit <- nsca_analysis(extract_data(), "X1", "Y", ceilings = "ce_fdh")
  row <- nsca_table(fit)
  for (param in c("d_nec", "d_suf", "weakest_effect", "balanced_joint_effect", "joint_empty_space_coverage",
                  "degenerate", "admissible_region_share",
                  "overlap_share", "weaker_component", "significance_status",
                  "equivalence_class", "necessity_corner")) {
    expect_equal(nsca_extract(fit, param = param), row[[param]][[1L]],
                 info = param)
  }
})

test_that("the prefixes reach the engines and disagree with each other", {
  dat <- extract_data()
  fit <- nsca_analysis(dat, "X1", "Y", direction = "HH", ceilings = "ce_fdh")

  expect_equal(
    nsca_extract(fit, param = "nec:Effect size"),
    as.numeric(NCA::nca_extract(
      NCA::nca_analysis(dat, "X1", "Y", corner = 1, ceilings = "ce_fdh"),
      "X1", "ce_fdh", "Effect size"
    ))
  )
  expect_equal(nsca_extract(fit, param = "nec:Effect size"),
               nsca_extract(fit, param = "d_nec"))
  expect_equal(nsca_extract(fit, param = "suf:Effect size"),
               nsca_extract(fit, param = "d_suf"))
  # The two sides are different analyses; a prefix that silently read the same
  # object would make them agree.
  expect_false(isTRUE(all.equal(
    nsca_extract(fit, param = "nec:Ceiling zone"),
    nsca_extract(fit, param = "suf:Ceiling zone")
  )))
})

test_that("defaults pick the first condition and the first frontier", {
  fit <- nsca_analysis(extract_data(), c("X1", "X2"), "Y",
                       ceilings = c("ce_fdh", "cr_fdh"))
  expect_equal(nsca_extract(fit), nsca_extract(fit, x = "X1", ceiling = "ce_fdh"))
  expect_identical(nsca_extract(fit, param = "condition"), "X1")
  expect_identical(nsca_extract(fit, param = "frontier"),
                   "ce_fdh")
  # And the non-defaults really do select something else.
  expect_identical(nsca_extract(fit, x = "X2", param = "condition"), "X2")
  expect_identical(
    nsca_extract(fit, ceiling = "cr_fdh", param = "frontier"),
    "cr_fdh"
  )
})

test_that("the retired name for the technique column is refused by name", {
  # `ceiling` named this column until 0.5.1. The retirement policy is the
  # one every other renamed column follows: extract() refuses it and says
  # what replaced it, rather than resolving it quietly, because a caller who
  # gets an answer never learns the name moved. The duplicate column is
  # still there under legacy = TRUE for code that cannot be changed yet.
  fit <- nsca_analysis(extract_data(), c("X1", "X2"), "Y",
                       ceilings = c("ce_fdh", "cr_fdh"))
  expect_error(nsca_extract(fit, param = "ceiling"), "frontier")
  legacy <- nsca_table(fit, legacy = TRUE)
  expect_true("ceiling" %in% names(legacy))
  expect_identical(legacy$ceiling, legacy$frontier)

  # `frontier` is the current name in every table, not a retired one.
  expect_false("frontier" %in% names(nsca_legacy_names()))
})

test_that("an unknown parameter is refused rather than returned as NA", {
  fit <- nsca_analysis(extract_data(), "X1", "Y", ceilings = "ce_fdh")
  expect_error(nsca_extract(fit, param = "d_total"), "Unknown parameter")
  expect_error(nsca_extract(fit, param = "d_total"), "nec:")
  expect_error(nsca_extract(fit, x = "X9"), "Unknown condition")
  # Since 0.5.0 the technique is checked before the row is built, so the
  # message names it and lists what was estimated instead of reporting the
  # absent row. Same wording as nsca_test_plot().
  expect_error(nsca_extract(fit, ceiling = "cr_fdh"),
               "was not estimated")
})

test_that("an engine parameter the engine does not know returns NA, not junk", {
  # The prefixed path goes straight to the engine, which has its own parameter
  # vocabulary. A name it does not recognize must come back as NA so that a
  # typo cannot be mistaken for a measurement.
  fit <- nsca_analysis(extract_data(), "X1", "Y", ceilings = "ce_fdh")
  expect_true(is.na(nsca_extract(fit, param = "nec:No Such Parameter")))
  expect_true(is.na(nsca_extract(fit, param = "suf:No Such Parameter")))
})

test_that("extracted values match the table for every condition and frontier", {
  fit <- nsca_analysis(extract_data(), c("X1", "X2"), "Y",
                       ceilings = c("ce_fdh", "cr_fdh"))
  table <- nsca_table(fit)
  for (i in seq_len(nrow(table))) {
    condition <- table$condition[[i]]
    ceiling <- table$frontier[[i]]
    expect_equal(
      nsca_extract(fit, x = condition, ceiling = ceiling, param = "weakest_effect"),
      table$weakest_effect[[i]],
      info = paste(condition, ceiling)
    )
  }
})

test_that("a frontier that was not estimated is refused on every route", {
  # The nec: and suf: routes read the engines directly and used to return NA
  # for a technique the model never fitted, which reads as a missing
  # measurement rather than a wrong request.
  fit <- nsca_analysis(extract_data(), "X1", "Y", ceilings = "ce_fdh")
  expect_error(nsca_extract(fit, ceiling = "cr_fdh", param = "nec:Effect size"),
               "was not estimated")
  expect_error(nsca_extract(fit, ceiling = "cr_fdh", param = "suf:Effect size"),
               "was not estimated")
  expect_error(nsca_extract(fit, ceiling = c("ce_fdh", "cr_fdh")),
               "was not estimated")
  expect_true(is.finite(nsca_extract(fit, param = "nec:Effect size")))
})
