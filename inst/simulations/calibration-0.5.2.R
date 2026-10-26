suppressPackageStartupMessages({
  library(NSCA)
  library(parallel)
})
# Type I error rate of each p-value in the joint table under complete random
# pairing: independent X and Y, n = 30, 99 shared permutations, nominal
# alpha 0.05, 2500 simulations. About 20 minutes on four cores.
# The simulations run in parallel where the platform allows it. The results do
# not depend on this, because each design, or each simulation, has its own
# seed.
cores <- if (.Platform$OS.type == "windows") 1L else 4L
sims <- 2500
n <- 30
reps <- 99
one <- function(s) {
  set.seed(50000 + s)
  dat <- data.frame(X = stats::runif(n), Y = stats::runif(n))
  fit <- suppressWarnings(nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
                                        shared.test.rep = reps,
                                        shared.seed = s))
  row <- nsca_table(fit)
  c(nec = row$p_nec, suf = row$p_suf, iut = row$p_nsca_iut,
    weakest = row$p_weakest_perm, sum = row$p_sum_perm)
}
t0 <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(sims), one, mc.cores = cores))
elapsed <- difftime(Sys.time(), t0, units = "mins")
se <- sqrt(0.05 * 0.95 / sims)
cat(sprintf("simulations %d, n %d, %d shared permutations, alpha 0.05 (SE %.5f), %.1f min\n",
            sims, n, reps, se, as.numeric(elapsed)))
for (k in colnames(res)) {
  r <- mean(res[, k] <= 0.05)
  ci <- r + c(-1, 1) * 1.96 * sqrt(r * (1 - r) / sims)
  ks <- if (k %in% c("nec", "suf", "weakest", "sum")) {
    suppressWarnings(stats::ks.test(res[, k], "punif")$p.value)
  } else NA
  cat(sprintf("  %-8s reject %.4f   95%% CI [%.4f, %.4f]%s\n", k, r, ci[1], ci[2],
              if (is.na(ks)) "" else sprintf("   KS p = %.2f", ks)))
}
cat(sprintf("  any of iut, weakest, sum: %.4f\n",
            mean(apply(res[, c("iut", "weakest", "sum")] <= 0.05, 1, any))))
