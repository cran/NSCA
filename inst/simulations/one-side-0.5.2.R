# Designs in which one side holds and the other is absent (NSCA 0.5.2).
# Because the claim that both sides hold is false in these designs,
# power_joint, the rejection rate of the absent side's test,
# power_weakest_perm and power_sum_perm are false-positive rates; the
# rejection rate of the side that holds is its power. 200 replications per
# row, 99 shared permutations. Monte Carlo SE of a rate near 0.05: about 0.015.
# About 25 minutes on two cores.
suppressPackageStartupMessages({
  library(NSCA)
  library(parallel)
})
# The designs run in parallel where the platform allows it. The results do not
# depend on this, because each design, or each simulation, has its own seed.
cores <- if (.Platform$OS.type == "windows") 1L else 2L
cols <- c("n", "d_nec_declared", "d_suf_declared", "power_necessity",
          "power_sufficiency", "power_joint", "power_weakest_perm",
          "power_sum_perm")
designs <- list(
  "only necessity holds (d_nec = 0.245, d_suf = 0)" = c(d_nec = 0.245, d_suf = 0),
  "only sufficiency holds (d_nec = 0, d_suf = 0.245)" = c(d_nec = 0, d_suf = 0.245)
)
out <- mclapply(designs, function(d) {
  nsca_power(n = c(30, 60), d_nec = d[["d_nec"]], d_suf = d[["d_suf"]],
             rep = 200, test.rep = 99, seed = 1, shared = TRUE)
}, mc.cores = cores)
for (label in names(out)) {
  cat(label, "\n")
  print(out[[label]][, cols], digits = 3, row.names = FALSE)
  cat("\n")
}
