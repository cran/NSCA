# Bands for the components, and no bands for the joint indices.
#
# The asymmetry is the point. Dul's bands are somebody else's published
# convention and citing them is ordinary; a band for `balanced_joint_effect`
# would be invented by the same package that defines the index, which is a
# claim with nothing behind it.

mag_data <- function(n = 70, seed = 21) {
  set.seed(seed)
  t <- sort(stats::runif(n) * 100)
  data.frame(
    Training = t,
    Luck = stats::runif(n) * 100,
    Performance = pmin(pmax(20 + 0.6 * t + stats::rnorm(n, 0, 9), 0), 100)
  )
}

mag_fit <- function(...) {
  suppressWarnings(suppressMessages(
    nsca_analysis(mag_data(), c("Training", "Luck"), "Performance",
                  ceilings = "ce_fdh", test.rep = 100, test.seed = 1, ...)
  ))
}


test_that("nsca_magnitude() is the SCAtools function", {
  expect_identical(nsca_magnitude(), SCAtools::sca_magnitude())
  expect_identical(nsca_magnitude(c(0.05, 0.2, 0.4, 0.7)),
                   SCAtools::sca_magnitude(c(0.05, 0.2, 0.4, 0.7)))
})

test_that("the components get a band, beside their floor", {
  out <- nsca_compare(mag_fit())
  expect_true("magnitude" %in% names(out))
  expect_identical(out$magnitude, nsca_magnitude(out$effect_size))
  expect_true(all(c("effect_size_null", "effect_size_excess") %in% names(out)))
})

test_that("a band and a floor together say more than either alone", {
  skip_on_cran()
  # The case this exists for. A noise condition can be labeled `medium` on
  # the conventional scale while more than half of the number is what random
  # pairing produces at that sample size. Both numbers, one row.
  out <- nsca_compare(mag_fit())
  luck <- out[out$condition == "Luck" & out$component == "necessity", ]
  # Not the specific label -- that depends on the draw. What matters is
  # that a value carrying a conventional band can have most of itself
  # accounted for by the floor, and that both numbers are in the row.
  expect_false(is.na(luck$magnitude))
  expect_true(luck$effect_size_null > luck$effect_size / 2)
  expect_identical(luck$magnitude, nsca_magnitude(luck$effect_size))
})

test_that("print() points at the bands for the components", {
  skip_on_cran()
  text <- paste(utils::capture.output(suppressWarnings(print(mag_fit()))),
                collapse = "\n")
  expect_match(text, "nsca_magnitude()", fixed = TRUE)
  expect_match(text, "Dul (2016)", fixed = TRUE)
  expect_match(text, "d_nec_null", fixed = TRUE)
})

test_that("print() still refuses a band for the joint indices", {
  skip_on_cran()
  # The line the package will not cross, and it has to stay crossable only by
  # somebody publishing a convention, not by this package inventing one.
  text <- paste(utils::capture.output(suppressWarnings(print(mag_fit()))),
                collapse = "\n")
  expect_match(text, "joint indices have no such convention", fixed = TRUE)
  expect_match(text, "calibrate by", fixed = TRUE)
  expect_match(text, "simulation before naming a value", fixed = TRUE)
})

test_that("summary() says the same about the joint indices", {
  skip_on_cran()
  text <- paste(utils::capture.output(suppressWarnings(
    print(summary(mag_fit()))
  )), collapse = "\n")
  expect_match(text, "established magnitude benchmarks", fixed = TRUE)
})

test_that("the joint table columns are not labeled by the band function", {
  skip_on_cran()
  # nsca_magnitude() will happily label any number handed to it -- it is a
  # lookup, not a guard. What must not happen is the package doing it for the
  # joint indices on the user's behalf.
  out <- suppressWarnings(nsca_table(mag_fit()))
  expect_false(any(grepl("magnitude", names(out))))
})
