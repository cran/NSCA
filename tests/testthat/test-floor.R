# A permutation test with too few replications, or too few observations, can
# never return a p-value at or below the threshold. Its silence then looks
# exactly like a null result, so the package says so when the test is set up.

floor_data <- function(n = 30) {
  set.seed(3)
  x <- stats::runif(n)
  data.frame(X = x, Y = pmin(pmax(x + stats::rnorm(n, 0, 0.1), 0), 1))
}

warnings_from <- function(expr) {
  seen <- character(0)
  withCallingHandlers(
    suppressMessages(utils::capture.output(expr)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  seen
}

test_that("the smallest p-value follows each sequence's own rule", {
  # The engines report max((count + 1) / (R + 1), 1 / R); the shared sequence
  # (count + 1) / (R + 1). Both cap R at n!.
  expect_equal(NSCA:::.nsca_smallest_p(19), 1 / 19)
  expect_equal(NSCA:::.nsca_smallest_p(19, shared = TRUE), 1 / 20)
  expect_equal(NSCA:::.nsca_smallest_p(1000, 4), 1 / 24)
  expect_equal(NSCA:::.nsca_smallest_p(1000, 4, shared = TRUE), 1 / 25)
})

test_that("nsca_analysis() warns once when no component can be significant", {
  seen <- warnings_from(
    nsca_analysis(floor_data(), "X", "Y", ceilings = "ce_fdh",
                  test.rep = 19, test.seed = 1)
  )
  # SCAtools' own warning for the sufficiency side is replaced, not repeated.
  expect_length(grep("too few", seen), 1L)
  expect_match(seen[grep("too few", seen)], "neither component")
  expect_match(seen[grep("too few", seen)], "test.rep >= 20")

  expect_length(
    grep("too few", warnings_from(
      nsca_analysis(floor_data(), "X", "Y", ceilings = "ce_fdh",
                    test.rep = 20, test.seed = 1)
    )),
    0L
  )
})

test_that("the shared sequence is checked against its own floor", {
  # (count + 1) / (R + 1) reaches .05 at R = 19, one fewer than the engines.
  seen <- warnings_from(
    nsca_analysis(floor_data(), "X", "Y", ceilings = "ce_fdh",
                  shared.test.rep = 18, shared.seed = 1)
  )
  expect_match(seen[grep("too few", seen)], "shared.test.rep >= 19")
  expect_length(
    grep("too few", warnings_from(
      nsca_analysis(floor_data(), "X", "Y", ceilings = "ce_fdh",
                    shared.test.rep = 19, shared.seed = 1)
    )),
    0L
  )
})

test_that("nsca_power() warns when test.rep forces an estimated power of 0", {
  expect_warning(
    nsca_power(n = 30, rep = 1, test.rep = 19, seed = 1),
    "number of permutations rather"
  )
  expect_warning(
    nsca_power(n = 30, rep = 1, test.rep = 19, seed = 1, shared = TRUE),
    NA
  )
  expect_error(nsca_power(n = 30, rep = 1, test.rep = 0), "at least 1")
  expect_error(nsca_power(n = 30, rep = 1, test.rep = 20, shared = NA),
               "TRUE or FALSE")
})

test_that("adjustment says the critical values are single-test lines", {
  set.seed(4)
  n <- 60
  data <- data.frame(A = stats::runif(n), B = stats::runif(n))
  data$Y <- pmin(pmax(data$A + stats::rnorm(n, 0, 0.1), 0), 1)
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    data, c("A", "B"), "Y", ceilings = "ce_fdh", test.rep = 40, test.seed = 1
  )))
  expect_message(nsca_table(fit, p.adjust = "holm"), "single-test reference")
  adjusted <- suppressMessages(nsca_table(fit, p.adjust = "holm"))
  plain <- nsca_table(fit)
  expect_identical(adjusted$d_nec_critical, plain$d_nec_critical)
  expect_identical(adjusted$d_suf_critical, plain$d_suf_critical)
})

test_that("the adjustment floor follows the shared sequence when there is one", {
  # 2 / (38 + 1) is above .05 and 2 / (39 + 1) is not; the engines' rule,
  # 2 / R, would need 40.
  expect_warning(
    NSCA:::.nsca_warn_adjust_floor(2, 38, 0.05, "bonferroni", shared = TRUE),
    "shared.test.rep >= 39"
  )
  expect_silent(
    NSCA:::.nsca_warn_adjust_floor(2, 39, 0.05, "bonferroni", shared = TRUE)
  )
  expect_warning(
    NSCA:::.nsca_warn_adjust_floor(2, 39, 0.05, "bonferroni"),
    "test.rep >= 40"
  )
})
