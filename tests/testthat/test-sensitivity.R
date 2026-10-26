# The scope is the denominator of both components and of every joint index
# built from them, so a series over boxes is the one diagnostic that says how
# much of a joint result was chosen rather than estimated.
#
# The direction is the thing worth pinning, because it runs opposite to the
# first guess and an earlier draft of the documentation had it backwards.

sens_data <- function(n = 80, seed = 11) {
  set.seed(seed)
  a <- sort(stats::runif(n) * 100)
  data.frame(
    A = a,
    B = stats::runif(n) * 100,
    Y = pmin(pmax(15 + 0.7 * a + stats::rnorm(n, 0, 7), 0), 100)
  )
}

sens_fit <- function(..., data = sens_data()) {
  suppressWarnings(suppressMessages(
    nsca_analysis(data, "A", "Y", ceilings = "ce_fdh", test.rep = 0, ...)
  ))
}


test_that("the series reports both components and both joint indices", {
  out <- nsca_scope_sensitivity(sens_fit())
  expect_true(all(c("d_nec", "d_suf", "weakest_effect",
                    "balanced_joint_effect", "admissible_region_share",
                    "geometry_acceptable") %in% names(out)))
  expect_identical(nrow(out), 4L)
  expect_identical(out$scope[[1L]], "observed range")
})

test_that("widening the box raises both components, not lowers them", {
  # The frontiers stay where the data is, so every square unit added beyond the
  # cloud falls on the empty side of one frontier or the other: the empty areas
  # grow faster than the box. This is the assertion an earlier draft of the
  # roxygen got backwards, which is why it is pinned rather than described.
  out <- nsca_scope_sensitivity(sens_fit())
  expect_true(all(diff(out$d_nec) > 0))
  expect_true(all(diff(out$d_suf) > 0))
  expect_true(all(diff(out$scope_area) > 0))
})

test_that("the admissible region falls toward zero as the box grows", {
  out <- nsca_scope_sensitivity(sens_fit())
  expect_true(all(diff(out$admissible_region_share) < 0))
  expect_true(out$admissible_region_share[[nrow(out)]] <
                out$admissible_region_share[[1L]] / 2)
})

test_that("the geometry flag turns FALSE on reconstruction, not on overlap", {
  # The mechanism, pinned because an earlier draft of the roxygen named the
  # wrong one. The two frontiers do not come close to crossing on a widened
  # box -- overlap_share *falls*. What gives way is reconstruction_error:
  # outside the observed range of X there are no data to fix a step frontier,
  # so the engines' extension and this package's reconstruction diverge and
  # the area identity stops closing.
  out <- nsca_scope_sensitivity(sens_fit())

  expect_true(out$geometry_acceptable[[1L]])
  expect_true(any(!out$geometry_acceptable))

  # The overlap never grows; the reconstruction error does, monotonically.
  expect_true(all(diff(out$overlap_share) <= 0))
  expect_true(all(diff(out$reconstruction_error) > 0))

  broken <- out[!out$geometry_acceptable, ]
  expect_true(all(broken$reconstruction_error > 0.02))
  expect_true(all(broken$overlap_share < 0.01))
})

test_that("a flagged row is unverified rather than known-wrong", {
  # The distinction the documentation makes, and the reason both diagnostic
  # columns are in the table: the areas on a flagged row may be correct, but
  # the package can no longer confirm that they describe the frontiers the
  # engines fitted. The components can still sit below 1 while the flag is
  # FALSE, which is exactly the case a "geometry failed" summary would hide.
  out <- nsca_scope_sensitivity(sens_fit())
  first_broken <- out[!out$geometry_acceptable, ][1L, ]
  expect_lt(first_broken$d_nec + first_broken$d_suf, 1)
})

test_that("the data does not change between rows, only the box", {
  # Every difference down a column is the scope. If the number of observations
  # moved, the series would be comparing two things at once.
  out <- nsca_scope_sensitivity(sens_fit())
  expect_length(unique(out$condition), 1L)
  expect_length(unique(out$frontier), 1L)
})

test_that("the box the model was fitted on is added as its own row", {
  # Without it the series answers a question the caller did not ask.
  fit <- sens_fit(scope = c(0, 100, 0, 100))
  out <- nsca_scope_sensitivity(fit, widen = c(0, 0.5))
  expect_true("as fitted (declared)" %in% out$scope)
  declared <- out[out$scope == "as fitted (declared)", ]
  expect_equal(c(declared$xmin, declared$xmax, declared$ymin, declared$ymax),
               c(0, 100, 0, 100))
})

test_that("an observed-scope model gets no declared row", {
  out <- nsca_scope_sensitivity(sens_fit(), widen = c(0, 0.5))
  expect_false("as fitted (declared)" %in% out$scope)
  expect_identical(nrow(out), 2L)
})

test_that("explicit boxes replace the widening series", {
  out <- nsca_scope_sensitivity(
    sens_fit(),
    scopes = list(tight = c(0, 100, 0, 100), wide = c(-50, 150, -50, 150))
  )
  expect_identical(out$scope, c("tight", "wide"))
  expect_true(out$d_nec[[2L]] > out$d_nec[[1L]])
})

test_that("the selection and the arguments are validated by name", {
  fit <- suppressWarnings(suppressMessages(
    nsca_analysis(sens_data(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  test.rep = 0)
  ))
  expect_error(nsca_scope_sensitivity(fit, x = c("A", "B")), "exactly one")
  expect_error(nsca_scope_sensitivity(fit, widen = -1), "non-negative")
  expect_error(nsca_scope_sensitivity(fit, scopes = list(c(1, 2))), "four")
  expect_error(nsca_scope_sensitivity(fit, ceiling = "ce_vrs"),
               "did not estimate")
})

test_that("it defaults to the first condition rather than refusing", {
  fit <- suppressWarnings(suppressMessages(
    nsca_analysis(sens_data(), c("A", "B"), "Y", ceilings = "ce_fdh",
                  test.rep = 0)
  ))
  expect_identical(unique(nsca_scope_sensitivity(fit)$condition), "A")
})


test_that("nsca_normalize() is the SCAtools function", {
  data <- sens_data()[1:5, ]
  expect_identical(
    nsca_normalize(data, scale = c(0, 1)),
    SCAtools::sca_normalize(data, scale = c(0, 1))
  )
})
