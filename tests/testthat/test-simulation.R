# nsca_random() and nsca_power() (F5).

test_that("nsca_random() leaves both corners empty", {
  set.seed(1)
  dat <- nsca_random(n = 200, d_nec = 0.35, d_suf = 0.35)
  fit <- nsca_analysis(dat, "X", "Y", direction = "HH", ceilings = "ce_fdh",
                       scope = c(0, 1, 0, 1))
  row <- nsca_table(fit)

  # A band leaves an empty region on both sides, so both components have
  # something to estimate. NCA::nca_random() fills one corner only.
  expect_true(row$d_nec[[1L]] > 0)
  expect_true(row$d_suf[[1L]] > 0)
})

test_that("larger population effect sizes give larger components", {
  set.seed(2)
  narrow <- nsca_random(n = 300, d_nec = 0.45, d_suf = 0.45)
  wide <- nsca_random(n = 300, d_nec = 0.28, d_suf = 0.28)

  d <- function(dat) {
    row <- nsca_table(nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
                                    scope = c(0, 1, 0, 1)))
    row$d_nec[[1L]] + row$d_suf[[1L]]
  }
  expect_true(d(narrow) > d(wide))
})

test_that("the two components can differ", {
  # Power planning for an asymmetric claim needs an asymmetric band.
  set.seed(3)
  dat <- nsca_random(n = 400, d_nec = 0.21, d_suf = 0.48)
  row <- nsca_table(nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
                                  scope = c(0, 1, 0, 1)))
  expect_false(isTRUE(all.equal(row$d_nec[[1L]], row$d_suf[[1L]],
                                tolerance = 0.05)))
})

test_that("the band leaves exactly the effect sizes it was given", {
  # Closed form where the line runs from the left edge to the top edge,
  # root-finding elsewhere; both must give back the requested areas.
  cases <- list(c(0.10, 0.10, 1, 1), c(0.32, 0.08, 1, 1), c(0.10, 0, 2, 1),
                c(0.30, 0.20, 0.2, 0.2), c(0.05, 0.40, 0, 0.5))
  for (case in cases) {
    lines <- NSCA:::.nsca_simulation_band(case[[1]], case[[2]], case[[3]], case[[4]])
    truth <- NSCA:::.nsca_band_truth(lines)
    expect_equal(c(truth$d_nec, truth$d_suf), case[1:2], tolerance = 1e-9,
                 info = paste(case, collapse = ", "))
  }
  # Slope 1, d = 0.10: the necessity line is y = (1 - sqrt(0.2)) + x.
  expect_equal(NSCA:::.nsca_simulation_band(0.10, 0.10, 1, 1)$upper,
               c(1 - sqrt(0.2), 1))
})

test_that("effect sizes that leave no band are refused", {
  expect_error(nsca_random(n = 50, d_nec = 0.6, d_suf = 0.6), "leave no band")
})

test_that("every direction is the HH pattern reflected", {
  set.seed(4)
  for (direction in c("HH", "LH", "HL", "LL")) {
    dat <- nsca_random(n = 150, direction = direction, d_nec = 0.38,
                       d_suf = 0.38)
    fit <- nsca_analysis(dat, "X", "Y", direction = direction,
                         ceilings = "ce_fdh", scope = c(0, 1, 0, 1))
    row <- nsca_table(fit)
    expect_true(row$d_nec[[1L]] > 0, info = direction)
    expect_true(row$d_suf[[1L]] > 0, info = direction)
  }
})

test_that("an invalid band is refused before any fitting happens", {
  expect_error(nsca_random(n = 50, d_nec = 1), "[0, 1)", fixed = TRUE)
  expect_error(nsca_random(n = 50, d_suf = -0.1), "[0, 1)", fixed = TRUE)
  expect_error(nsca_random(n = 50, slope_nec = -1), "at least 0")
  expect_error(nsca_random(n = 2), "at least 3")
})


test_that("nsca_power() reports three rates", {
  # Deliberately tiny: this test is about the shape of the answer, not about
  # the value. Every replication fits two engines and runs two permutation
  # tests, so a realistic design belongs in a script, not a test suite.
  out <- suppressWarnings(
    nsca_power(n = 60, d_nec = 0.38, d_suf = 0.38, rep = 3, test.rep = 20, seed = 1)
  )

  expect_equal(nrow(out), 1L)
  for (column in c("power_necessity", "power_sufficiency", "power_joint",
                   "power_if_independent", "rejection_phi", "replications",
                   "shortfall")) {
    expect_true(column %in% names(out), info = column)
  }
  # The joint test requires both components, so it can never beat either.
  expect_true(out$power_joint <= out$power_necessity)
  expect_true(out$power_joint <= out$power_sufficiency)
  expect_equal(out$power_if_independent,
               out$power_necessity * out$power_sufficiency)
})

test_that("the rates are over the replications that produced a verdict", {
  # Dividing by `rep` would count a failed replication as a non-rejection,
  # which biases every rate down and hides the failures at the same time.
  out <- suppressWarnings(
    nsca_power(n = 60, d_nec = 0.38, d_suf = 0.38, rep = 3, test.rep = 20, seed = 1)
  )
  expect_equal(out$replications + out$shortfall, 3L)
})

test_that("nsca_power() restores the caller's random stream", {
  set.seed(77)
  before <- runif(1)
  set.seed(77)
  invisible(suppressWarnings(
    nsca_power(n = 60, d_nec = 0.38, d_suf = 0.38, rep = 2, test.rep = 20, seed = 5)
  ))
  expect_equal(before, runif(1))
})

test_that("a malformed band fails before the first fit", {
  expect_error(
    nsca_power(n = 60, d_nec = 0.6, d_suf = 0.6, rep = 1, test.rep = 10),
    "leave no band"
  )
})


# ---- what the band declares, against what a sample realizes (1e) --------

test_that("nsca_power() reports the declared effect sizes and the realized ones", {
  out <- suppressWarnings(
    nsca_power(n = 60, d_nec = 0.32, d_suf = 0.32, rep = 3, test.rep = 20, seed = 1)
  )
  for (column in c("d_nec_declared", "d_suf_declared",
                   "d_nec_realized", "d_suf_realized")) {
    expect_true(column %in% names(out), info = column)
  }
  expect_true(out$d_nec_declared > 0 && out$d_nec_declared < 1)
  expect_true(out$d_suf_declared > 0 && out$d_suf_declared < 1)
})

test_that("the realized effect size overshoots the declared one", {
  skip_on_cran()
  # Measured, not reasoned: the estimate is biased *upward*. The FDH ceiling
  # is the running maximum of the sample, which sits below the population
  # boundary, so the empty region above it is larger than the true one --
  # and the same argument runs the other way for the floor. Both empty areas
  # are therefore overstated, and the overstatement shrinks as n grows.
  out <- suppressWarnings(
    nsca_power(n = c(40, 200), d_nec = 0.32, d_suf = 0.32, rep = 4, test.rep = 20, seed = 2)
  )
  expect_true(all(out$d_nec_realized >= out$d_nec_declared))
  expect_true(all(out$d_suf_realized >= out$d_suf_declared))
  # And it converges: the larger sample is closer to the truth.
  gap <- out$d_nec_realized - out$d_nec_declared
  expect_true(gap[[2L]] < gap[[1L]])
})

test_that("the declared effect sizes do not depend on the direction", {
  # The other three directions are the same band reflected, and a reflection
  # moves an area without changing it.
  hh <- suppressWarnings(
    nsca_power(n = 50, d_nec = 0.36, d_suf = 0.36, rep = 1, test.rep = 10, seed = 3)
  )
  ll <- suppressWarnings(
    nsca_power(n = 50, direction = "LL", d_nec = 0.36, d_suf = 0.36, rep = 1,
               test.rep = 10, seed = 3)
  )
  expect_equal(hh$d_nec_declared, ll$d_nec_declared)
  expect_equal(hh$d_suf_declared, ll$d_suf_declared)
})


# ---- dependence, measured rather than assumed (2e) ---------------------

test_that("rejection_phi replaces the difference from the product", {
  skip_on_cran()
  out <- suppressWarnings(
    nsca_power(n = 60, d_nec = 0.38, d_suf = 0.38, rep = 4, test.rep = 20, seed = 4)
  )
  expect_true("rejection_phi" %in% names(out))
  expect_false("dependence_gap" %in% names(out))
  # A correlation, so it is bounded, or NA when one test never varied.
  expect_true(is.na(out$rejection_phi) ||
                (out$rejection_phi >= -1 && out$rejection_phi <= 1))
})

test_that("phi is NA when a component rejected always or never", {
  expect_true(is.na(NSCA:::.nsca_rejection_phi(c(TRUE, TRUE), c(TRUE, FALSE))))
  expect_true(is.na(NSCA:::.nsca_rejection_phi(logical(0), logical(0))))
})

test_that("phi is 1 when the two components reject on exactly the same draws", {
  nec <- c(TRUE, TRUE, FALSE, FALSE)
  expect_equal(NSCA:::.nsca_rejection_phi(nec, nec), 1)
  expect_equal(NSCA:::.nsca_rejection_phi(nec, !nec), -1)
})


# ---- failed replications are redrawn, capped, reported (9d) ------------

test_that("the shortfall is reported rather than the failures being silent", {
  out <- suppressWarnings(
    nsca_power(n = 60, d_nec = 0.38, d_suf = 0.38, rep = 3, test.rep = 20, seed = 5)
  )
  expect_true(all(c("replications", "attempts", "shortfall") %in% names(out)))
  expect_equal(out$shortfall, 3L - out$replications)
  # Redrawing means attempts may exceed the replications obtained, never the
  # other way round.
  expect_true(out$attempts >= out$replications)
})


# ---- nsca_powerplot() ---------------------------------------------------

# Computed once: every test below reads the same table, and each run costs
# twelve fits with their permutation tests.
pp_cache <- new.env()
pp_table <- function() {
  if (is.null(pp_cache$table)) {
    pp_cache$table <- suppressWarnings(
      nsca_power(n = c(40, 80), d_nec = 0.245, d_suf = 0.245, rep = 3, test.rep = 20, seed = 8)
    )
  }
  pp_cache$table
}

test_that("nsca_powerplot() draws one curve per component", {
  skip_on_cran()
  p <- nsca_powerplot(pp_table())
  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  expect_equal(length(levels(p$data$component)), 3L)
})

test_that("the joint curve is drawn first, so the others show through", {
  # Reversed in 0.5.0 after looking at the figure. Joint used to be the last
  # level and was drawn on top of the two components it is built from -- and
  # at high power all three coincide, so a plot with three legend entries
  # showed one curve. Joint is now the first level: drawn underneath, with
  # the two component curves dashed over it, and it leads the legend because
  # it is the curve a study is powered on.
  p <- nsca_powerplot(pp_table())
  expect_true("joint (both reject)" %in% levels(p$data$component))
  expect_identical(levels(p$data$component)[[1L]], "joint (both reject)")
})

test_that("coincident power curves stay readable", {
  # The reason for the order above: distinct linetypes and shapes, so three
  # curves at 1.00 are still three curves.
  p <- nsca_powerplot(pp_table())
  expect_true(any(vapply(p$scales$scales,
                         function(s) "linetype" %in% s$aesthetics,
                         logical(1L))))
  expect_true(any(vapply(p$scales$scales,
                         function(s) "shape" %in% s$aesthetics,
                         logical(1L))))
})

test_that("the independence benchmark is a separate unkeyed layer", {
  # Its distance from the joint curve is the dependence made visible, and it
  # is not a fourth test, so it stays out of the legend.
  with_ref <- nsca_powerplot(pp_table(), reference = TRUE)
  without <- nsca_powerplot(pp_table(), reference = FALSE)
  expect_true(length(with_ref$layers) > length(without$layers))
})

test_that("the target line can be turned off", {
  tbl <- pp_table()
  expect_true(
    length(nsca_powerplot(tbl, target = 0.8)$layers) >
      length(nsca_powerplot(tbl, target = NULL)$layers)
  )
})

test_that("a subset of the rates may be drawn", {
  p <- nsca_powerplot(pp_table(), show = "joint", reference = FALSE)
  expect_equal(length(levels(p$data$component)), 1L)
})

test_that("a table that is not a power table is refused by name", {
  tbl <- pp_table()
  expect_error(nsca_powerplot(tbl[, c("n", "power_joint")]),
               "Missing power-result column")
})

test_that("the caption carries the overshoot caveat", {
  # The curves inherit it from nsca_power(): the band is filled uniformly and
  # the realized effect sizes overshoot at small n.
  expect_match(nsca_powerplot(pp_table())$labels$caption, "overshoot")
})


# ---- reproducibility and the supplementary tests (0.5.2) ----------------

test_that("nsca_power() reproduces with a seed", {
  # Until 0.5.2 the engines' tests reseeded the session from the clock in the
  # middle of the simulation, so every replication after the first drew its
  # data from an unseeded stream.
  first <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.32, d_suf = 0.32, rep = 2, test.rep = 20, seed = 9)
  )
  second <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.32, d_suf = 0.32, rep = 2, test.rep = 20, seed = 9)
  )
  expect_identical(first, second)
})

test_that("shared = TRUE adds the power of the two supplementary tests", {
  plain <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.32, d_suf = 0.32, rep = 2, test.rep = 20, seed = 9)
  )
  expect_true(all(is.na(plain$power_weakest_perm)))
  expect_true(all(is.na(plain$power_sum_perm)))

  shared <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.32, d_suf = 0.32, rep = 2, test.rep = 20, seed = 9,
               shared = TRUE)
  )
  for (column in c("power_weakest_perm", "power_sum_perm")) {
    expect_true(is.finite(shared[[column]]), info = column)
    expect_true(shared[[column]] >= 0 && shared[[column]] <= 1, info = column)
  }
  # The weakest effect rejects only when the smaller area is extreme, which
  # the joint test also needs of both; neither bound is guaranteed, so only
  # the shape is checked here.
  expect_identical(names(shared), names(plain))
})

test_that("nsca_powerplot() draws the supplementary rates only when asked", {
  skip_on_cran()
  shared <- suppressWarnings(
    nsca_power(n = c(30, 50), d_nec = 0.32, d_suf = 0.32, rep = 2, test.rep = 20, seed = 9,
               shared = TRUE)
  )
  p <- nsca_powerplot(shared, show = c("joint", "weakest", "sum"))
  expect_equal(length(levels(p$data$component)), 3L)
  expect_match(p$labels$caption, "supplementary")
  expect_equal(length(levels(nsca_powerplot(shared)$data$component)), 3L)
  expect_error(nsca_powerplot(pp_table(), show = "sum"), "shared = TRUE")
})

test_that("the two sides can be planned with different effect sizes", {
  out <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.32, d_suf = 0.08, rep = 1, test.rep = 20,
               seed = 1)
  )
  expect_equal(nrow(out), 1L)
  expect_equal(c(out$d_nec_declared, out$d_suf_declared), c(0.32, 0.08))
})

test_that("a slope changes the boundary, not the declared effect size", {
  steep <- suppressWarnings(
    nsca_power(n = 30, d_nec = 0.10, slope_nec = 2, rep = 1, test.rep = 20,
               seed = 1)
  )
  expect_equal(steep$d_nec_declared, 0.10)
  expect_equal(NSCA:::.nsca_simulation_band(0.10, 0.10, 2, 1)$upper[[1L]],
               1 - sqrt(2 * 2 * 0.10))
})
