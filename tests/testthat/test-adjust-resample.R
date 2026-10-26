# The NSCA halves of the three items pulled back into 0.5.0. Each differs from
# its SCAtools counterpart in one way, and that way is what is tested here:
# two component families instead of one, two floors instead of one, and four
# resampled quantities instead of one because the joint indices inherit both
# overstatements.

ar_data <- function(n = 60, seed = 11) {
  set.seed(seed)
  a <- sort(stats::runif(n))
  data.frame(
    A = a,
    B = stats::runif(n),
    C = stats::runif(n),
    Y = pmin(pmax(0.15 + 0.7 * a + stats::rnorm(n, 0, 0.07), 0), 1)
  )
}

# The fit is deterministic (test.seed = 1), so a call with the same arguments
# returns the model already built instead of rerunning 100 permutations for
# each of three conditions. Most tests below ask for the same two fits.
ar_fit <- local({
  cache <- list()
  function(x = c("A", "B", "C"), ..., data = ar_data()) {
    key <- list(x = x, dots = list(...), data = data)
    for (entry in cache) {
      if (identical(entry$key, key)) {
        return(entry$fit)
      }
    }
    fit <- suppressWarnings(suppressMessages(
      nsca_analysis(data, x, "Y", ceilings = "ce_fdh",
                    test.rep = 100, test.seed = 1, ...)
    ))
    cache[[length(cache) + 1L]] <<- list(key = key, fit = fit)
    fit
  }
})


## ------------------------------------------------- multiple-testing p.adjust

test_that("the default changes nothing", {
  fit <- ar_fit()
  expect_identical(names(nsca_table(fit)),
                   names(nsca_table(fit, p.adjust = "none")))
  expect_false("p_nec_adjusted" %in% names(nsca_table(fit)))
})

test_that("each component family is adjusted separately", {
  fit <- ar_fit()
  raw <- nsca_table(fit)
  out <- nsca_table(fit, p.adjust = "BH")

  expect_equal(out$p_nec_adjusted, stats::p.adjust(raw$p_nec, method = "BH"))
  expect_equal(out$p_suf_adjusted, stats::p.adjust(raw$p_suf, method = "BH"))
})

test_that("the joint p is rebuilt from the pair, not adjusted as a maximum", {
  # p.adjust(p_nsca_iut) is a different and wrong number. The IUT rejects only
  # when both components reject, so the coherent adjusted decision is
  # max(p_nec_adjusted, p_suf_adjusted) -- and only that agrees with the
  # component columns printed beside it.
  #
  # Pinned by the definition rather than by a difference from p.adjust(p_iut).
  # The two constructions coincide on many samples -- they differ only when the
  # ranking of a component across conditions differs from the ranking of the
  # maximum -- so asserting that they disagree would pass or fail with the
  # data. What must always hold is the identity below.
  fit <- ar_fit()
  out <- nsca_table(fit, p.adjust = "BH")

  expect_equal(out$p_nsca_iut_adjusted,
               pmax(out$p_nec_adjusted, out$p_suf_adjusted))
  # And therefore never below either component, which is what makes the joint
  # decision the more conservative of the two.
  expect_true(all(out$p_nsca_iut_adjusted >= out$p_nec_adjusted - 1e-12))
  expect_true(all(out$p_nsca_iut_adjusted >= out$p_suf_adjusted - 1e-12))
})

test_that("an adjustment can only make a decision more conservative", {
  fit <- ar_fit()
  raw <- nsca_table(fit)
  out <- nsca_table(fit, p.adjust = "bonferroni")

  expect_true(all(out$p_nec_adjusted >= raw$p_nec - 1e-12))
  expect_true(all(out$p_suf_adjusted >= raw$p_suf - 1e-12))
  # Nothing that failed before can pass now.
  gained <- which(!isTRUE(raw$joint_support) & isTRUE(out$joint_support))
  expect_length(gained, 0L)
})

test_that("every significance and support column is re-derived", {
  fit <- ar_fit()
  out <- nsca_table(fit, p.adjust = "bonferroni")

  expect_identical(out$necessity_significant, out$p_nec_adjusted <= 0.05)
  expect_identical(out$sufficiency_significant, out$p_suf_adjusted <= 0.05)
  expect_identical(out$nsca_significant,
                   out$necessity_significant & out$sufficiency_significant)
  expect_false(any(out$significance_status == "supported" &
                     !out$nsca_significant))
})

test_that("the relevance half of support survives an adjustment", {
  skip_on_cran()
  fit <- ar_fit(relevance = c(necessity = 0.2, sufficiency = 0.2))
  raw <- nsca_table(fit)
  out <- nsca_table(fit, p.adjust = "bonferroni")
  expect_identical(raw$d_nec, out$d_nec)
  expect_identical(raw$d_suf, out$d_suf)
})

test_that("an unknown method is refused", {
  expect_error(nsca_table(ar_fit(), p.adjust = "banana"), "Unknown p.adjust")
})


## ----------------------------------------------------- the two independence floors

test_that("the joint table reports a floor for each component", {
  out <- nsca_table(ar_fit("A"))
  expect_true(all(c("d_nec_null", "d_suf_null",
                    "d_nec_excess", "d_suf_excess") %in% names(out)))
  expect_equal(out$d_nec_excess, out$d_nec - out$d_nec_null)
  expect_equal(out$d_suf_excess, out$d_suf - out$d_suf_null)
})

test_that("each floor is the mean of that component's own null", {
  fit <- ar_fit("A")
  out <- nsca_table(fit)
  expect_equal(out$d_nec_null,
               mean(fit$necessity$tests[["A"]][["ce_fdh"]]$data))
  expect_equal(out$d_suf_null,
               mean(fit$sufficiency$tests[["A"]][["ce_fdh"]]$data))
})

test_that("the component tables carry the same two columns", {
  # They come from two different places -- the necessity rows are built here,
  # the sufficiency rows come from SCAtools -- so a column present on one side
  # and absent on the other is the failure mode worth pinning.
  tables <- nsca_results(ar_fit("A"))
  for (side in c("necessity", "sufficiency")) {
    expect_true(
      all(c("effect_size_null", "effect_size_excess") %in% names(tables[[side]])),
      info = side
    )
  }
})

test_that("with no test there is no floor", {
  out <- nsca_table(suppressWarnings(suppressMessages(
    nsca_analysis(ar_data(), "A", "Y", ceilings = "ce_fdh", test.rep = 0)
  )))
  expect_true(is.na(out$d_nec_null))
  expect_true(is.na(out$d_suf_excess))
})


## -------------------------------------------------------- nsca_resample()

test_that("both components and both joint indices are resampled", {
  skip_on_cran()
  fit <- ar_fit("A", data = ar_data(n = 100))
  out <- nsca_resample(fit, sizes = c(0.4, 0.8), reps = 15, seed = 1)

  expect_setequal(unique(out$quantity),
                  c("d_nec", "d_suf", "weakest_effect",
                    "balanced_joint_effect"))
  expect_identical(nrow(out), 8L)
  expect_setequal(unique(out$m), c(40L, 80L))
})

test_that("the joint index inherits both overstatements", {
  skip_on_cran()
  # The reason all four are reported. A small sample inflates each component,
  # and balanced_joint_effect is a geometric mean of two inflated numbers, so
  # it does not cancel -- it carries both.
  fit <- ar_fit("A", data = ar_data(n = 120))
  out <- nsca_resample(fit, sizes = c(0.3, 0.7), reps = 20, seed = 1)
  joint <- out[out$quantity == "balanced_joint_effect", ]
  joint <- joint[order(joint$m), ]

  expect_true(all(joint$subsample_median > joint$observed))
  expect_true(diff(joint$subsample_median) < 0)
})

test_that("the observed column is the full-sample value of that quantity", {
  fit <- ar_fit("A")
  full <- nsca_table(fit)
  out <- nsca_resample(fit, sizes = 0.5, reps = 10, seed = 1)
  for (q in c("d_nec", "d_suf", "weakest_effect", "balanced_joint_effect")) {
    expect_equal(out$observed[out$quantity == q], full[[q]], info = q)
  }
})

test_that("no confidence interval is produced here either", {
  out <- nsca_resample(ar_fit("A"), sizes = 0.5, reps = 8, seed = 1)
  expect_false(any(grepl("conf|ci_|lower|upper", names(out))))
})

test_that("a full-size draw is taken once", {
  out <- nsca_resample(ar_fit("A"), sizes = 1, reps = 30, seed = 1)
  expect_true(all(out$reps == 1L))
  expect_true(all(is.na(out$subsample_sd)))
  expect_equal(out$subsample_median, out$observed)
})

test_that("the caller's random stream is given back", {
  fit <- ar_fit("A")
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(nsca_resample(fit, sizes = 0.5, reps = 4, seed = 3))
  expect_identical(before, stats::runif(1))
})

test_that("more than one condition is refused", {
  expect_error(nsca_resample(ar_fit(), x = c("A", "B")), "exactly one")
})


## ------------------------------------------- the replication floor warning

test_that("the adjustment floor names the count that is needed", {
  # The arithmetic of the test below without its eight fitted conditions, so
  # that it also runs where the slower test is skipped. Under BH a condition
  # that is the only strong one among 8 is multiplied by 8: 8 / 50 cannot
  # reach .05, and 8 / 160 can.
  expect_warning(NSCA:::.nsca_warn_adjust_floor(8L, 50, 0.05, "BH"),
                 "resolution of the test")
  expect_warning(NSCA:::.nsca_warn_adjust_floor(8L, 50, 0.05, "BH"),
                 "test.rep >= 160")
  expect_warning(NSCA:::.nsca_warn_adjust_floor(8L, 159, 0.05, "BH"),
                 "test.rep >= 160")
  expect_silent(NSCA:::.nsca_warn_adjust_floor(8L, 160, 0.05, "BH"))
})

test_that("an adjustment that cannot reject says so", {
  skip_on_cran()
  # Found by the adversarial pass. The engines' p cannot fall below
  # 1/test.rep, and a condition that is the only strong one among k is
  # multiplied by at least k, so unless k/test.rep <= alpha it cannot reach
  # alpha however strong the evidence. Eight conditions on a low count looks
  # exactly like a null result and is not one.
  set.seed(4)
  n <- 80
  data <- data.frame(Y = stats::runif(n))
  for (i in 1:8) data[[paste0("C", i)]] <- stats::runif(n)
  data$C1 <- pmin(pmax(data$Y + stats::rnorm(n, 0, 0.1), 0), 1)

  thin <- suppressWarnings(suppressMessages(nsca_analysis(
    data, paste0("C", 1:8), "Y", ceilings = "ce_fdh",
    test.rep = 50, test.seed = 1
  )))
  expect_warning(suppressMessages(nsca_table(thin, p.adjust = "BH")),
                 "resolution of the test")
  # 8 / 159 is just above .05 and 8 / 160 is not: until 0.5.2 the warning
  # used 1/(test.rep + 1) and asked for one replication too few.
  expect_warning(suppressMessages(nsca_table(thin, p.adjust = "BH")),
                 "test.rep >= 160")

  # Enough replications, no warning -- and the one real condition survives.
  thick <- suppressWarnings(suppressMessages(nsca_analysis(
    data, paste0("C", 1:8), "Y", ceilings = "ce_fdh",
    test.rep = 400, test.seed = 1
  )))
  out <- suppressMessages(nsca_table(thick, p.adjust = "BH"))
  expect_identical(sum(out$joint_support, na.rm = TRUE), 1L)
})

test_that("the floor warning is silent when nothing is adjusted", {
  set.seed(4)
  n <- 60
  data <- data.frame(A = stats::runif(n), B = stats::runif(n))
  data$Y <- pmin(pmax(data$A + stats::rnorm(n, 0, 0.1), 0), 1)
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    data, c("A", "B"), "Y", ceilings = "ce_fdh", test.rep = 20, test.seed = 1
  )))
  expect_silent(nsca_table(fit))
})

test_that("the floor warning is silent when no test was run", {
  set.seed(4)
  n <- 60
  data <- data.frame(A = stats::runif(n), B = stats::runif(n))
  data$Y <- pmin(pmax(data$A + stats::rnorm(n, 0, 0.1), 0), 1)
  fit <- suppressWarnings(suppressMessages(nsca_analysis(
    data, c("A", "B"), "Y", ceilings = "ce_fdh", test.rep = 0
  )))
  expect_silent(nsca_table(fit, p.adjust = "BH"))
})
