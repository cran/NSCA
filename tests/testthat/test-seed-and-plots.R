# test.seed, and the plotting that had no counterpart in this package until
# 0.5.0.

sp_data <- function(n = 60, seed = 6) {
  set.seed(seed)
  x <- sort(runif(n))
  data.frame(X = x, Y = pmin(pmax(x + rnorm(n, 0, 0.08), 0), 1))
}


# ---- test.seed -----------------------------------------------------------

test_that("test.seed reproduces both component p-values", {
  # p_nec and p_suf decide p_nsca_iut and therefore joint_support. Without a
  # seed the engines draw from the clock and the verdict moves between runs.
  fit <- function() {
    nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                  test.rep = 40, test.seed = 7)
  }
  a <- nsca_table(fit())
  b <- nsca_table(fit())

  expect_equal(a$p_nec, b$p_nec)
  expect_equal(a$p_suf, b$p_suf)
  expect_equal(a$p_nsca_iut, b$p_nsca_iut)
})

test_that("one seed gives the two components different sequences", {
  # Seeding both engines with the same number would hand them the same
  # shuffle: each derives its whole sequence from a single draw off the
  # caller's stream. That is the shared-sequence design, reached by accident
  # and with p_source still reporting two independent runs.
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 60, test.seed = 7)
  nec <- as.numeric(fit$necessity$tests[["X"]][["ce_fdh"]]$data)
  suf <- as.numeric(fit$sufficiency$tests[["X"]][["ce_fdh"]]$data)

  expect_equal(length(nec), length(suf))
  expect_false(isTRUE(all.equal(nec, suf)))
})

test_that("test.seed does not reseed the caller's stream", {
  # Built first: sp_data() calls set.seed() itself, so generating it inside
  # the block would reset the stream and test the helper rather than the fit.
  dat <- sp_data()
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
                          test.rep = 20, test.seed = 3))
  expect_equal(before, runif(1))
})

test_that("a non-integer test.seed is refused", {
  expect_error(
    nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                  test.rep = 10, test.seed = "one"),
    "test.seed"
  )
})

test_that("the seed is recorded so print() can report it", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 10, test.seed = 42)
  expect_equal(attr(fit, "nsca")$test.seed, 42)
})


# ---- nsca_test_plot() ----------------------------------------------------

test_that("nsca_test_plot() draws both component nulls", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 40, test.seed = 1)
  p <- nsca_test_plot(fit)

  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  # One panel per component.
  expect_equal(length(levels(droplevels(p$data$component))), 2L)
})

test_that("the subtitle carries the intersection-union p-value", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 40, test.seed = 1)
  p <- nsca_test_plot(fit)
  expect_match(p$labels$subtitle, "p_nsca_iut", fixed = TRUE)
})

test_that("the caption says these two panels are not a joint null", {
  # Two independent shuffles cannot give a null for a statistic of both
  # components; that is what nsca_shared_plot() is for.
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 40, test.seed = 1)
  expect_match(nsca_test_plot(fit)$labels$caption, "nsca_shared_plot",
               fixed = TRUE)
})

test_that("an untested model says test.rep, not 'no such frontier'", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh")
  expect_error(nsca_test_plot(fit), "test.rep", fixed = TRUE)
})

test_that("an unestimated frontier is named as such", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 20, test.seed = 1)
  expect_error(nsca_test_plot(fit, ceiling = "cr_fdh"), "not estimated")
})


# ---- one technique draws no linetype key --------------------------------

test_that("a single technique maps no linetype, two still do", {
  # With one technique the key distinguishes nothing and ggplot draws its
  # sample line in black, which reads as a third component beside a color
  # legend where every color is one of the two.
  line_aes <- function(p) {
    layers <- Filter(function(l) inherits(l$geom, "GeomLine"), p$layers)
    if (length(layers) == 0L) character(0) else names(layers[[1L]]$mapping)
  }

  one <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh")
  two <- nsca_analysis(sp_data(), "X", "Y", ceilings = c("ce_fdh", "cr_fdh"))

  expect_false("linetype" %in% line_aes(nsca_plot(one)))
  expect_true("linetype" %in% line_aes(nsca_plot(two)))
})


# ---- plot text is wrapped and the caption is left-aligned ---------------

test_that("captions are left-aligned in both plot families", {
  # ggplot2 right-aligns a caption, and a wrapped right-aligned block has a
  # ragged left edge that no reader tracks back to.
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh")
  expect_equal(nsca_plot(fit)$theme$plot.caption$hjust, 0)

  sca_fit <- SCAtools::sca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh")
  expect_equal(SCAtools::sca_plot(sca_fit)$theme$plot.caption$hjust, 0)
})

test_that("annotate = FALSE removes the corner labels", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh")
  with_labels <- length(nsca_plot(fit)$layers)
  without <- length(nsca_plot(fit, annotate = FALSE)$layers)
  # Two corners, so two labels.
  expect_equal(with_labels - without, 2L)
})


# ---- nsca_output() -------------------------------------------------------

test_that("nsca_output() returns the pieces it drew", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       test.rep = 20, test.seed = 1)
  out <- suppressWarnings(
    utils::capture.output(
      result <- nsca_output(fit, thresholds = TRUE, plots = FALSE,
                            tests = FALSE)
    )
  )
  expect_true(is.list(result))
  expect_true("thresholds" %in% names(result))
})


# ---- listwise deletion reaches the engines ------------------------------

test_that("listwise fits both components on the same rows", {
  # Until 0.5.0 na.action changed this package's bookkeeping and nothing the
  # engines reported, because each engine drops incomplete pairs per condition
  # inside its own loop.
  dat <- sp_data()
  dat$Z <- dat$X
  dat$X[c(3, 9)] <- NA
  dat$Z[c(4, 11)] <- NA

  fit <- nsca_analysis(dat, c("X", "Z"), "Y", ceilings = "ce_fdh",
                       na.action = "listwise")
  counts <- nsca_table(fit)$observations
  expect_equal(length(unique(counts)), 1L)
})


# ---- qr as a reference line (F3) ----------------------------------------

test_that("nsca_reference() delegates and reports both lines", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       reference = c("ols", "qr"))
  ref <- nsca_reference(fit)

  expect_equal(nrow(ref), 2L)
  expect_setequal(ref$line, c("ols", "qr"))
  expect_equal(ref$tau[ref$line == "qr"], 0.5)
})

test_that("it is the same line SCAtools would report", {
  # Delegated rather than reimplemented: the sufficiency layer is an
  # sca_result fitted on exactly these observations.
  dat <- sp_data()
  fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh", reference = "ols")
  expect_equal(
    nsca_reference(fit)$slope,
    SCAtools::sca_reference(fit$sufficiency)$slope
  )
})

test_that("the plot names both unkeyed lines", {
  fit <- nsca_analysis(sp_data(), "X", "Y", ceilings = "ce_fdh",
                       reference = c("ols", "qr"))
  caption <- nsca_plot(fit)$labels$caption
  expect_match(caption, "OLS central tendency")
  expect_match(caption, "quantile regression")
})
