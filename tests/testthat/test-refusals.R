# The NSCA half of the refusals added in 0.5.0. Same six inputs, and each one
# is worse here than in SCAtools because two engines act on it rather than one.

ref_data <- function(n = 40, seed = 3) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  data.frame(X = x, Y = pmin(pmax(x + stats::rnorm(n, 0, 0.08), 0), 1))
}


test_that("an unknown frontier technique is refused before either engine runs", {
  # It reached both engines, each fell back to its own defaults, and the
  # metadata still recorded the name that was asked for -- so the two
  # components were fitted with a technique the Settings block did not name.
  expect_error(
    nsca_analysis(ref_data(), "X", "Y", ceilings = "banana", test.rep = 0),
    "Unknown frontier technique"
  )
  expect_error(
    nsca_analysis(ref_data(), "X", "Y", ceilings = c("ce_fdh", "nope"),
                  test.rep = 0),
    "'nope'"
  )
})

test_that("the pinned technique list is the engine's own", {
  engine <- get("ceilings", envir = asNamespace("NCA"))
  ours <- get(".nsca_known_ceilings", envir = asNamespace("NSCA"))
  expect_setequal(ours, engine)
})

test_that("the two packages pin the same list", {
  # They must, or a technique accepted by one and refused by the other would
  # make the component tables unreachable from the joint one.
  expect_setequal(
    get(".nsca_known_ceilings", envir = asNamespace("NSCA")),
    get(".sca_known_ceilings", envir = asNamespace("SCAtools"))
  )
})

test_that("every valid frontier technique is still accepted", {
  for (technique in c("ce_fdh", "ce_vrs", "cr_fdh", "cr_vrs")) {
    expect_silent(suppressWarnings(suppressMessages(
      nsca_analysis(ref_data(), "X", "Y", ceilings = technique, test.rep = 0)
    )))
  }
})


test_that("both replication counts are validated", {
  # A string that parses to a number is coerced rather than refused.
    for (bad in list(-5, 1.5, NA, "not a number")) {
    expect_error(
      nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = bad),
      "test.rep"
    )
    expect_error(
      nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 0,
                    shared.test.rep = bad),
      "shared.test.rep"
    )
  }
})

test_that("a significance threshold outside (0, 1) is refused", {
  # Here it decides three columns at once: necessity_significant,
  # sufficiency_significant and the joint verdict built from them.
  for (bad in list(2, 0, 1, -0.5, NA)) {
    expect_error(
      nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 0,
                    test.p_threshold = bad),
      "test.p_threshold"
    )
  }
})


test_that("a reversed or zero-width declared scope is refused", {
  # The scope is the denominator of both components here, so a silently
  # repaired box moves d_nec, d_suf and every joint index at once.
  expect_error(
    nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 0,
                  scope = c(1, 0, 0, 1)),
    "minimum"
  )
  expect_error(
    nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 0,
                  scope = c(0, 1, 0.5, 0.5)),
    "no width"
  )
  expect_silent(suppressWarnings(suppressMessages(
    nsca_analysis(ref_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 0,
                  scope = c(0, NA, 0, 1))
  )))
})


test_that("a variable may not be both a condition and the outcome", {
  # Worse than in SCAtools: both frontiers land on the same diagonal, so the
  # two components and the admissible region between them are all fixed by the
  # geometry rather than by the data.
  expect_error(
    nsca_analysis(ref_data(), "X", "X", test.rep = 0),
    "both a condition and the outcome"
  )
})


test_that("a failed column selection names the argument and the name", {
  expect_error(nsca_analysis(ref_data(), "X", "Z", test.rep = 0), "'y' names")
  expect_error(nsca_analysis(ref_data(), "Z", "Y", test.rep = 0), "'x' names")
  expect_error(nsca_analysis(ref_data(), "Z", "Y", test.rep = 0), "Available")
  expect_error(nsca_analysis(ref_data(), c("X", "X"), "Y", test.rep = 0),
               "more than once")
  expect_error(nsca_analysis("hello", "X", "Y", test.rep = 0),
               "data frame or a matrix")
  expect_error(nsca_analysis(ref_data()[0, ], "X", "Y", test.rep = 0),
               "no rows")
})


test_that("the two packages refuse the same things the same way", {
  # One vocabulary for refusals as well as for columns: a caller who learns
  # SCAtools' messages should not have to learn a second set.
  data <- ref_data()
  cases <- list(
    list(args = list(data, "X", "Z"),                       match = "'y' names"),
    list(args = list(data, "Z", "Y"),                       match = "'x' names"),
    list(args = list(data, c("X", "X"), "Y"),               match = "more than once"),
    list(args = list(data, "X", "X"),                       match = "both a condition"),
    list(args = list("hello", "X", "Y"),                    match = "data frame or a matrix")
  )
  for (case in cases) {
    expect_error(do.call(SCAtools::sca_analysis, c(case$args, list(test.rep = 0))),
                 case$match, info = case$match)
    expect_error(do.call(nsca_analysis, c(case$args, list(test.rep = 0))),
                 case$match, info = case$match)
  }
})
