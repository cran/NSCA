# Power of the joint test and of the two supplementary tests (NSCA 0.5.2).
# Three population designs, two sample sizes, 100 replications each, and 99
# shared permutations per replication. About 15 minutes on three cores.
suppressPackageStartupMessages({
  library(NSCA)
  library(parallel)
})
# The designs run in parallel where the platform allows it. The results do not
# depend on this, because each design, or each simulation, has its own seed.
cores <- if (.Platform$OS.type == "windows") 1L else 3L
cols <- c("n", "d_nec_declared", "d_suf_declared", "power_necessity",
          "power_sufficiency", "power_joint", "power_weakest_perm",
          "power_sum_perm")
designs <- list(
  "A both corners equally empty, d = 0.165" = c(d_nec = 0.165, d_suf = 0.165),
  "B one corner clearly empty, the other barely: d = 0.245 and 0.011" =
    c(d_nec = 0.245, d_suf = 0.011),
  "C both corners moderately empty, d = 0.10" = c(d_nec = 0.10, d_suf = 0.10)
)
out <- mclapply(designs, function(d) {
  nsca_power(n = c(15, 30), d_nec = d[["d_nec"]], d_suf = d[["d_suf"]],
             rep = 100, test.rep = 99, seed = 1, shared = TRUE)
}, mc.cores = cores)
for (label in names(out)) {
  cat(label, "\n")
  print(out[[label]][, cols], digits = 3, row.names = FALSE)
  cat("\n")
}
cat("Monte Carlo SE of a rate near 0.5 with 100 replications: about 0.05\n")
