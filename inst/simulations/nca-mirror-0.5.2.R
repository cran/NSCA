# The same tendency in the test of the NCA package (NSCA 0.5.2). Only
# sufficiency holds in the data (the upper-left corner is not empty in the
# population), and NCA::nca_analysis() tests the upper-left corner. 200
# samples per row (2000 for the unrelated reference), with 99 permutations
# each. About 20 minutes on one core. Columns:
#   d         observed NCA effect size (observed scope)
#   null      mean of its random-pairing distribution
#   critical  NCA's threshold value at alpha = 0.05
#   reject    share of samples with p <= 0.05
#   corner    cases with X < 0.3 and Y > 0.7, observed and under random
#             pairing (mean of 50 shuffles); high_y: share of Y above 0.7
suppressPackageStartupMessages({
  library(NCA)
  library(NSCA)
})
# NCA reseeds R's random number stream from the clock when its permutation
# test finishes. Running it under a seed drawn from the stream, and restoring
# the stream afterward, keeps the runs reproducible.
with_seed <- function(seed, expr) {
  previous <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", previous, envir = globalenv()))
  set.seed(seed)
  expr
}
one <- function(dat) {
  seed <- sample.int(.Machine$integer.max, 1L)
  invisible(utils::capture.output(m <- suppressWarnings(with_seed(seed,
    nca_analysis(dat, "X", "Y", ceilings = "ce_fdh", test.rep = 99)))))
  test <- m$tests$X$ce_fdh
  corner <- function(x, y) sum(x < 0.3 & y > 0.7)
  c(d = test$observed, null = mean(test$data), critical = test$threshold.value,
    reject = test$p_value <= 0.05,
    corner_observed = corner(dat$X, dat$Y),
    corner_shuffled = mean(replicate(50, corner(dat$X, sample(dat$Y)))),
    high_y = mean(dat$Y > 0.7))
}
report <- function(label, generate, samples = 200) {
  set.seed(7)
  r <- t(replicate(samples, one(generate())))
  cat(sprintf(paste0("%-28s d %.3f  null %.3f  critical %.3f  reject %.3f  ",
                     "corner cases %.2f observed vs %.2f shuffled  high_y %.2f\n"),
              label, mean(r[, "d"]), mean(r[, "null"]), mean(r[, "critical"]),
              mean(r[, "reject"]), mean(r[, "corner_observed"]),
              mean(r[, "corner_shuffled"]), mean(r[, "high_y"])))
}
report("only sufficiency, n = 30",
       function() nsca_random(30, d_nec = 0, d_suf = 0.245))
report("only sufficiency, n = 60",
       function() nsca_random(60, d_nec = 0, d_suf = 0.245))
report("unrelated uniform, n = 30",
       function() data.frame(X = stats::runif(30), Y = stats::runif(30)),
       samples = 2000)
