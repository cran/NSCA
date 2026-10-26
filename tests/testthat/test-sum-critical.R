# 0.5.2: p_sum_perm, the critical values beside every null, the Monte Carlo
# margin of every p-value, and NCA's sampling rule in the shared sequence.

sc_data <- function(n = 40, noise = 0.1, seed = 8) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  data.frame(X = x, Y = pmin(pmax(x + stats::rnorm(n, 0, noise), 0), 1))
}

sc_shared <- function(data = sc_data(), reps = 49, seed = 1, ...) {
  suppressWarnings(nsca_analysis(data, "X", "Y", ceilings = "ce_fdh",
                                 shared.test.rep = reps, shared.seed = seed,
                                 ...))
}

test_that("p_sum_perm is the permutation p-value of component_sum", {
  fit <- sc_shared()
  row <- nsca_table(fit)
  draws <- attr(fit, "nsca")$shared$sum[, 1L, 1L]
  expect_equal(row$p_sum_perm,
               (1 + sum(draws >= row$component_sum)) / (length(draws) + 1))
  expect_gte(row$p_sum_perm, 1 / (length(draws) + 1))
})

test_that("without a shared sequence the two supplementary tests are NA", {
  row <- nsca_table(suppressWarnings(nsca_analysis(
    sc_data(), "X", "Y", ceilings = "ce_fdh", test.rep = 49, test.seed = 1
  )))
  for (column in c("p_weakest_perm", "p_sum_perm", "p_weakest_perm_accuracy",
                   "p_sum_perm_accuracy", "component_sum_critical",
                   "weakest_effect_critical")) {
    expect_true(is.na(row[[column]]), info = column)
  }
  # The component tests are the engines' own, and so are their criticals.
  expect_true(is.finite(row$d_nec_critical))
  expect_true(is.finite(row$d_suf_critical))
})

test_that("each critical value is the 1 - alpha quantile of its own null", {
  fit <- sc_shared()
  row <- nsca_table(fit)
  shared <- attr(fit, "nsca")$shared
  q <- function(v) unname(stats::quantile(v, 0.95))
  expect_equal(row$d_nec_critical, q(shared$necessity[, 1L, 1L]))
  expect_equal(row$d_suf_critical, q(shared$sufficiency[, 1L, 1L]))
  expect_equal(row$weakest_effect_critical, q(shared$weakest[, 1L, 1L]))
  expect_equal(row$component_sum_critical, q(shared$sum[, 1L, 1L]))
  # The sum has its own critical value, not the sum of the two.
  expect_false(isTRUE(all.equal(row$component_sum_critical,
                                row$d_nec_critical + row$d_suf_critical)))
})

test_that("floors and critical values come from the same draws as p", {
  # Shared sequence only: the engines ran no test, and the floor must still
  # be there, read from the shuffles the p-values came from.
  fit <- sc_shared()
  row <- nsca_table(fit)
  shared <- attr(fit, "nsca")$shared
  expect_equal(row$d_nec_null, mean(shared$necessity[, 1L, 1L]))
  expect_equal(row$d_suf_null, mean(shared$sufficiency[, 1L, 1L]))
})

test_that("the engine route reads criticals and accuracy from the engines", {
  fit <- suppressWarnings(nsca_analysis(sc_data(), "X", "Y",
                                        ceilings = "ce_fdh", test.rep = 49,
                                        test.seed = 1))
  row <- nsca_table(fit)
  tables <- nsca_results(fit)
  expect_equal(row$p_nec_accuracy, tables$necessity$p_accuracy)
  expect_equal(row$p_suf_accuracy, tables$sufficiency$p_accuracy)
  expect_equal(row$d_nec_critical, tables$necessity$effect_size_critical)
  expect_equal(row$d_suf_critical, tables$sufficiency$effect_size_critical)
  expect_true(is.finite(tables$necessity$effect_size_critical))
})

test_that("shared accuracies follow NCA's p-accuracy formula", {
  fit <- sc_shared(test.p_confidence = 0.95)
  row <- nsca_table(fit)
  reps <- dim(attr(fit, "nsca")$shared$necessity)[[1L]]
  margin <- function(p) stats::qnorm(0.975) * sqrt(p * (1 - p) / reps)
  expect_equal(row$p_nec_accuracy, margin(row$p_nec))
  expect_equal(row$p_sum_perm_accuracy, margin(row$p_sum_perm))
})

test_that("a count above n! is lowered to n! and every ordering is used", {
  skip_on_cran()
  tiny <- data.frame(X = c(0.1, 0.3, 0.5, 0.7, 0.9),
                     Y = c(0.2, 0.3, 0.6, 0.8, 0.9))
  expect_message(
    fit <- suppressWarnings(nsca_analysis(tiny, "X", "Y", ceilings = "ce_fdh",
                                          shared.test.rep = 200,
                                          shared.seed = 1)),
    "Lowered shared.test.rep to 120"
  )
  shared <- attr(fit, "nsca")$shared
  expect_equal(dim(shared$necessity)[[1L]], 120L)
  expect_true(attr(shared, "exhaustive"))
  expect_identical(attr(fit, "nsca")$shared.test.rep, 120L)
  row <- nsca_table(fit)
  expect_identical(row$p_nec_accuracy, 0)
  expect_identical(row$p_sum_perm_accuracy, 0)

  # Every ordering, so the null no longer depends on the seed.
  other <- suppressMessages(suppressWarnings(nsca_analysis(
    tiny, "X", "Y", ceilings = "ce_fdh", shared.test.rep = 200,
    shared.seed = 99
  )))
  expect_equal(sort(shared$sum[, 1L, 1L]),
               sort(attr(other, "nsca")$shared$sum[, 1L, 1L]))
})

test_that("up to 720 permutations the orderings are distinct", {
  skip_on_cran()
  # Five observations have 120 orderings; asking for exactly that many must
  # use each once, which free sampling almost never does (it would repeat
  # about 44 of them).
  five <- data.frame(X = c(0.1, 0.3, 0.5, 0.7, 0.9),
                     Y = c(0.2, 0.3, 0.6, 0.8, 0.9))
  expect_no_message(
    fit <- suppressWarnings(nsca_analysis(five, "X", "Y", ceilings = "ce_fdh",
                                          shared.test.rep = 120,
                                          shared.seed = 3))
  )
  expect_true(attr(attr(fit, "nsca")$shared, "exhaustive"))
})

test_that("p.adjust adjusts p_sum_perm with the other families", {
  skip_on_cran()
  data <- sc_data()
  data$X2 <- stats::runif(nrow(data))
  fit <- suppressWarnings(nsca_analysis(data, c("X", "X2"), "Y",
                                        ceilings = "ce_fdh",
                                        shared.test.rep = 49,
                                        shared.seed = 1))
  out <- nsca_table(fit, p.adjust = "holm")
  expect_equal(out$p_sum_perm_adjusted,
               stats::p.adjust(out$p_sum_perm, method = "holm"))
  expect_identical(match("p_sum_perm_adjusted", names(out)),
                   match("p_sum_perm_accuracy", names(out)) + 1L)
})

test_that("the plots draw a dashed critical line on every panel", {
  skip_on_cran()
  fit <- sc_shared()
  p <- nsca_shared_plot(fit)
  dashed <- Filter(function(l) inherits(l$geom, "GeomVline") &&
                     identical(l$aes_params$linetype, "dashed"), p$layers)
  expect_length(dashed, 1L)
  expect_equal(nrow(dashed[[1L]]$data), 4L)
  expect_equal(nlevels(droplevels(ggplot2::ggplot_build(p)$plot$data$component)),
               4L)

  engine <- suppressWarnings(nsca_analysis(sc_data(), "X", "Y",
                                           ceilings = "ce_fdh", test.rep = 49,
                                           test.seed = 1))
  q <- nsca_test_plot(engine)
  dashed <- Filter(function(l) inherits(l$geom, "GeomVline") &&
                     identical(l$aes_params$linetype, "dashed"), q$layers)
  expect_length(dashed, 1L)
  expect_equal(sort(dashed[[1L]]$data$critical),
               sort(c(nsca_table(engine)$d_nec_critical,
                      nsca_table(engine)$d_suf_critical)))
})

test_that("print says which p-value decides and lists the supplementary tests", {
  out <- paste(utils::capture.output(print(sc_shared())), collapse = "\n")
  expect_match(out, "joint_support uses p_nsca_iut")
  expect_match(out, "p_sum_perm      the two corners together")
  expect_match(out, "supplementary tests on the same shuffles")
})
