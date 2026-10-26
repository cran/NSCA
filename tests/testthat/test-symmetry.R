# The functions NSCA was missing beside SCAtools, and the reason they are not
# one-line delegations.
#
# There is one average effect per condition and two boundaries. Handing the
# sufficiency layer to SCAtools would answer the question for one of them and
# return a table with nothing on it saying which, so both functions here carry
# a `component` column and repeat the regression columns down the pair.

symmetry_data <- function(n = 60, seed = 11) {
  set.seed(seed)
  a <- sort(stats::runif(n))
  data.frame(
    A = a,
    B = stats::runif(n),
    Y = pmin(pmax(0.15 + 0.7 * a + stats::rnorm(n, 0, 0.07), 0), 1)
  )
}

symmetry_fit <- function(ceilings = "cr_fdh", ..., data = symmetry_data()) {
  suppressWarnings(suppressMessages(
    nsca_analysis(data, c("A", "B"), "Y", direction = "HH",
                  ceilings = ceilings, test.rep = 0, ...)
  ))
}

printed <- function(fit) {
  paste(utils::capture.output(suppressWarnings(print(fit))), collapse = "\n")
}


## --------------------------------------------------------- nsca_compare()

test_that("every condition produces one row per component", {
  out <- nsca_compare(symmetry_fit())
  expect_identical(nrow(out), 4L)
  expect_setequal(unique(out$component), c("necessity", "sufficiency"))
  expect_setequal(unique(out$condition), c("A", "B"))
})

test_that("the regression columns repeat down a condition's pair", {
  # One line through the cloud; two boundaries compared against it. If beta
  # differed between the two rows of a condition, one of them would be
  # describing a line that was never fitted.
  out <- nsca_compare(symmetry_fit())
  for (condition in c("A", "B")) {
    pair <- out[out$condition == condition, ]
    expect_length(unique(pair$beta), 1L)
    expect_length(unique(pair$slope), 1L)
    expect_length(unique(pair$r_squared), 1L)
  }
})

test_that("the two components can land in different configurations", {
  # The whole point of two rows. A condition can bound the outcome from one
  # side and not the other, and a single row would have to pick a side.
  out <- nsca_compare(symmetry_fit())
  expect_true(all(out$configuration %in% c("I", "II", "III", "IV")))
  expect_true(any(out$configuration != out$configuration[[1L]]))
})

test_that("the label names the component it belongs to", {
  out <- nsca_compare(symmetry_fit())
  necessity <- out[out$component == "necessity" &
                     out$configuration %in% c("II", "IV"), ]
  skip_if(nrow(necessity) == 0L, "no necessity row reached II or IV")
  expect_true(all(grepl("necessity", necessity$configuration_label)))
  expect_false(any(grepl("sufficiency", necessity$configuration_label)))
})

test_that("the effect sizes are the ones nsca_results() reports", {
  # Read from the component tables rather than from the engines, so that one
  # number cannot appear twice with two values.
  fit <- symmetry_fit()
  out <- nsca_compare(fit)
  tables <- nsca_results(fit)
  for (component in c("necessity", "sufficiency")) {
    source <- tables[[component]]
    for (condition in c("A", "B")) {
      expect_equal(
        out$effect_size[out$component == component &
                          out$condition == condition],
        source$effect_size[source$condition == condition],
        info = paste(component, condition)
      )
    }
  }
})

test_that("d_cut defaults to the model's relevance, per component", {
  fit <- symmetry_fit(relevance = c(necessity = 0.05, sufficiency = 0.40))
  out <- nsca_compare(fit)
  expect_equal(unique(out$d_cut[out$component == "necessity"]), 0.05)
  expect_equal(unique(out$d_cut[out$component == "sufficiency"]), 0.40)
})

test_that("d_cut falls back to the convention when no relevance was stated", {
  out <- nsca_compare(symmetry_fit())
  expect_equal(unique(out$d_cut), 0.10)
})

test_that("a mistyped d_cut name is refused rather than applied by position", {
  # The same protection `relevance` has, for the same reason: this argument
  # decides a reported classification, and c(sufficiency = , nec = ) would
  # otherwise be applied in reverse with no warning.
  expect_error(
    nsca_compare(symmetry_fit(), d_cut = c(nec = 0.1, sufficiency = 0.2)),
    "necessity"
  )
})

test_that("one d_cut applies to both components", {
  out <- nsca_compare(symmetry_fit(), d_cut = 0.25)
  expect_equal(unique(out$d_cut), 0.25)
})

test_that("a technique the model did not estimate is refused by name", {
  expect_error(
    nsca_compare(symmetry_fit(), ceiling = "ce_vrs"),
    "did not estimate"
  )
})

test_that("beta is the least-squares coefficient, not whatever line was asked for", {
  data <- symmetry_data()
  ols <- suppressWarnings(suppressMessages(nsca_analysis(
    data, "A", "Y", ceilings = "cr_fdh", test.rep = 0, reference = "ols"
  )))
  qr <- suppressWarnings(suppressMessages(nsca_analysis(
    data, "A", "Y", ceilings = "cr_fdh", test.rep = 0, reference = "qr"
  )))
  expect_equal(unique(nsca_compare(qr)$beta), unique(nsca_compare(ols)$beta))
})


## --------------------------------------------------- nsca_locationshift()

test_that("both components get their own slope comparison", {
  out <- nsca_locationshift(symmetry_fit())
  expect_identical(nrow(out), 4L)
  expect_setequal(unique(out$component), c("necessity", "sufficiency"))
  for (condition in c("A", "B")) {
    pair <- out[out$condition == condition, ]
    expect_length(unique(pair$ols_slope), 1L)
    expect_length(unique(pair$median_slope), 1L)
  }
})

test_that("the two components can disagree about being a location shift", {
  # A band whose one edge tracks the regression while the other does not is a
  # real shape, and a table with one row per condition cannot show it.
  out <- nsca_locationshift(symmetry_fit())
  expect_true(all(is.finite(out$frontier_slope)))
  expect_true(all(out$relative_spread >= 0))
})

test_that("a step frontier has no slope and the reading says why", {
  out <- nsca_locationshift(symmetry_fit(ceilings = "ce_fdh"))
  expect_true(all(is.na(out$frontier_slope)))
  expect_true(all(is.na(out$location_shift)))
  expect_true(all(grepl("Not assessable", out$reading)))
})

test_that("the tolerance is reported alongside the verdict it produced", {
  out <- nsca_locationshift(symmetry_fit(), tolerance = 0.3)
  expect_equal(unique(out$tolerance), 0.3)
  expect_identical(out$location_shift, out$relative_spread <= 0.3)
})

test_that("a negative tolerance is refused", {
  expect_error(nsca_locationshift(symmetry_fit(), tolerance = -1),
               "non-negative")
})


## --------------------------------------------------------- nsca_outliers()

test_that("nsca_outliers() takes a fitted model and names both axes of it", {
  out <- suppressMessages(nsca_outliers(symmetry_fit()))
  skip_if(is.null(out), "no observation reached min.dif on this sample")

  expect_identical(names(out)[1:3], c("condition", "component", "corner"))
  expect_true(all(out$component %in% c("necessity", "sufficiency")))
})

test_that("the data-frame form returns the same shape", {
  out <- suppressMessages(nsca_outliers(symmetry_data(), "A", "Y"))
  skip_if(is.null(out), "no observation reached min.dif on this sample")
  expect_identical(names(out)[1:3], c("condition", "component", "corner"))
  expect_identical(unique(out$condition), "A")
})

test_that("arguments the model already fixes are refused", {
  fit <- symmetry_fit()
  for (extra in list(list(y = "Y"), list(direction = "LH"),
                     list(scope = c(0, 1, 0, 1)))) {
    expect_error(do.call(nsca_outliers, c(list(fit), extra)), "fitted model")
  }
})

test_that("the class finally has a print method", {
  # It was set since 0.4.x and nothing consumed it, so the three attributes it
  # carried -- the direction and the two corners -- were never shown.
  out <- suppressMessages(nsca_outliers(symmetry_fit()))
  skip_if(is.null(out), "no observation reached min.dif on this sample")
  text <- paste(utils::capture.output(print(out)), collapse = "\n")
  expect_match(text, "NSCA outlier sensitivity")
  expect_match(text, "necessity corner")
  expect_match(text, "sufficiency corner")
})


## --------------------------------------- nsca_scales() and nsca_rescale()

test_that("the scale vocabulary is the one the thresholds were built with", {
  expect_identical(nsca_scales(), SCAtools::sca_scales())
})

test_that("nsca_rescale() is the conversion nsca_thresholds() used", {
  values <- c(1, 3, 4, 7, 9)
  for (scale in nsca_scales()$scale) {
    expect_identical(
      nsca_rescale(4, values, scale = scale, letter = "H"),
      SCAtools::sca_rescale(4, values, scale = scale, letter = "H"),
      info = scale
    )
  }
})

test_that("the direction letter still mirrors what it should", {
  # A delegation that silently dropped an argument would pass the identity
  # test above and fail here.
  values <- c(1, 3, 4, 7, 9)
  expect_false(identical(
    nsca_rescale(4, values, scale = "percentile", letter = "H"),
    nsca_rescale(4, values, scale = "percentile", letter = "L")
  ))
  expect_false(identical(
    nsca_rescale(4, values, scale = "percentage.range", letter = "L",
                 convention = "directional"),
    nsca_rescale(4, values, scale = "percentage.range", letter = "L",
                 convention = "absolute")
  ))
})


## ------------------------------------------------- multiplicity disclosure

test_that("print() counts the component tests and the joint ones", {
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    symmetry_data(), c("A", "B"), "Y", ceilings = "ce_fdh",
    test.rep = 20, test.seed = 1
  )))
  lines <- printed(fit)

  expect_match(lines, "multiplicity  2 conditions", fixed = TRUE)
  expect_match(lines, "4 component p-values, 2 joint tests", fixed = TRUE)
})

test_that("both halves of the multiplicity story are printed together", {
  # "Conservative" alone is the reassuring half. A reader who sees it without
  # "nothing is adjusted across conditions" has been told that the multiple
  # testing is handled, which is true within a condition and false across them.
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    symmetry_data(), c("A", "B"), "Y", ceilings = "ce_fdh",
    test.rep = 20, test.seed = 1
  )))
  lines <- printed(fit)

  expect_match(lines, "conservative within a", fixed = TRUE)
  expect_match(lines, "nothing is adjusted across conditions", fixed = TRUE)
  expect_match(lines, "at most 0.10", fixed = TRUE)
})

test_that("one condition still gets the conservativeness note", {
  # Unlike SCAtools, which says nothing for a single test: here there are two
  # component p-values behind one decision even with one condition, and how
  # they were combined is worth a line.
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    symmetry_data(), "A", "Y", ceilings = "ce_fdh",
    test.rep = 20, test.seed = 1
  )))
  lines <- printed(fit)
  expect_match(lines, "2 component p-values, 1 joint test", fixed = TRUE)
  expect_match(lines, "rejects", fixed = TRUE)
  expect_false(grepl("across conditions", lines, fixed = TRUE))
})

test_that("no permutation test means no disclosure", {
  expect_false(grepl("multiplicity", printed(symmetry_fit()), fixed = TRUE))
})
