# Why the test of an absent corner may reject when the opposite corner is
# empty (NSCA 0.5.2). For each design, 200 samples (2000 for the unrelated
# reference, whose rejection rate is compared with 0.05), with 99 permutations
# each. Columns:
#   d_suf     observed sufficiency effect size (observed scope)
#   null      mean of its random-pairing distribution
#   critical  its 1 - alpha quantile
#   reject    share of samples with p_suf <= 0.05 (and p_nec for necessity)
#   corner    cases with X > 0.7 and Y < 0.3, observed and under random
#             pairing (mean of 50 shuffles)
# In the first design, the sufficiency corner is absent in the population.
# About 30 minutes on one core.
suppressPackageStartupMessages(library(NSCA))
# The engines reseed R's random number stream from the clock when a
# permutation test finishes; drawing each test seed from the stream keeps the
# runs reproducible.
one <- function(dat) {
  seed <- sample.int(.Machine$integer.max, 1L)
  invisible(utils::capture.output(fit <- suppressWarnings(suppressMessages(
    nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh", test.rep = 99,
                  test.seed = seed)))))
  tab <- suppressWarnings(nsca_table(fit))
  corner <- function(x, y) sum(x > 0.7 & y < 0.3)
  c(d_suf = tab$d_suf, null = tab$d_suf_null, critical = tab$d_suf_critical,
    reject_suf = tab$p_suf <= 0.05, reject_nec = tab$p_nec <= 0.05,
    corner_observed = corner(dat$X, dat$Y),
    corner_shuffled = mean(replicate(50, corner(dat$X, sample(dat$Y)))))
}
report <- function(label, generate, samples = 200) {
  set.seed(7)
  r <- t(replicate(samples, one(generate())))
  cat(sprintf(paste0("%-34s d_suf %.3f  null %.3f  critical %.3f  ",
                     "reject suf %.3f  reject nec %.3f  ",
                     "corner cases %.2f observed vs %.2f shuffled\n"),
              label, mean(r[, "d_suf"]), mean(r[, "null"]),
              mean(r[, "critical"]), mean(r[, "reject_suf"]),
              mean(r[, "reject_nec"]), mean(r[, "corner_observed"]),
              mean(r[, "corner_shuffled"])))
}
report("only necessity, n = 30",
       function() nsca_random(30, d_nec = 0.245, d_suf = 0))
report("unrelated uniform, n = 30",
       function() data.frame(X = stats::runif(30), Y = stats::runif(30)),
       samples = 2000)
report("correlated normal r = .5, n = 60", function() {
  x <- stats::rnorm(60)
  data.frame(X = x, Y = 0.5 * x + sqrt(0.75) * stats::rnorm(60))
})
