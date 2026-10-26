library(NSCA)

set.seed(42)

# A relationship that is close to deterministic, in which both components are
# expected to hold.
n <- 120
x <- sort(runif(n))
tight <- data.frame(X = x, Y = pmin(pmax(x + rnorm(n, 0, 0.08), 0), 1))

fit <- nsca_analysis(
  tight, x = "X", y = "Y",
  direction = "HH",
  ceilings = c("ce_fdh", "cr_fdh"),
  test.rep = 1000
)

fit                      # NCA, SCA, and joint NSCA results
summary(fit)
nsca_plot(fit)           # both frontiers, all techniques, band shaded

# Dual thresholds: the X ranges in which an outcome level is unattainable,
# possible, or guaranteed.
nsca_thresholds(fit, ceiling = "ce_fdh")
nsca_thresholds(fit, ceiling = "ce_fdh", scale = "percentile")

# As the relationship is loosened, the three summaries change in the same
# direction but not by the same amount. All three are retained because they
# address different questions.
for (w in c(0.05, 0.20, 0.40)) {
  noisy <- data.frame(X = x, Y = pmin(pmax(x + runif(n, -w, w), 0), 1))
  row <- nsca_table(nsca_analysis(noisy, "X", "Y", ceilings = "ce_fdh"))
  cat(sprintf(
    "noise %.2f:  weakest_effect = %.3f  balanced_joint_effect = %.3f  joint_empty_space_coverage = %.3f\n",
    w, row$weakest_effect, row$balanced_joint_effect, row$joint_empty_space_coverage
  ))
}

# The identity that all area quantities satisfy, and the residual used to
# check it.
row <- nsca_table(fit)[1L, ]
cat(sprintf(
  "\nidentity: 1 - d_nec - d_suf + O = %.6f   admissible_region_share = %.6f\n",
  1 - row$d_nec - row$d_suf + row$overlap_share, row$admissible_region_share
))
cat(sprintf("reconstruction_error = %.2e   geometry_acceptable = %s\n",
            row$reconstruction_error, row$geometry_acceptable))

# For a perfect diagonal, the bounds hold with equality rather than only as
# upper limits.
for (nn in c(20L, 40L, 200L)) {
  d <- nsca_table(nsca_analysis(
    data.frame(X = seq(0, 1, length.out = nn), Y = seq(0, 1, length.out = nn)),
    "X", "Y", ceilings = "ce_fdh"
  ))
  cat(sprintf(
    "n = %3d:  O = %.5f (1/(n-1) = %.5f)  weakest = %.5f (bound %.5f)  balanced = %.5f (bound %.5f)\n",
    nn, d$overlap_share, 1 / (nn - 1),
    d$weakest_effect, 0.5 + d$overlap_share / 2,
    d$balanced_joint_effect, 1 + d$overlap_share
  ))
}

# The shared permutation sequence and the support decision. Both engines use
# the same shuffle in each replication, so the null distribution behind
# p_weakest_perm consists of values that a single dataset could produce.
shared <- nsca_analysis(
  tight, "X", "Y", ceilings = "ce_fdh",
  shared.test.rep = 199, shared.seed = 1,
  relevance = c(necessity = 0.10, sufficiency = 0.10)
)
print(nsca_table(shared)[c(
  "p_source", "p_nec", "p_suf", "p_nsca_iut", "p_weakest_perm",
  "necessity_supported", "sufficiency_supported", "geometry_acceptable",
  "joint_support"
)])

# p_weakest_perm serves as a sensitivity analysis rather than as a replacement
# for the intersection-union rule: its null hypothesis is random pairing alone,
# which is a proper subset of the union null "not necessary or not sufficient".

# Retired column names remain available.
print(nsca_legacy_names())
legacy <- nsca_table(fit, legacy = TRUE)
stopifnot(identical(legacy$d_nsca, legacy$weakest_effect))

# The three summaries differ under curvature. Each curve below is an exact
# bijection, so each is perfectly necessary and sufficient.
# joint_empty_space_coverage reflects this; the other two summaries indicate how
# unequal the two empty spaces are.
grid <- seq(0.001, 1, length.out = 200)
for (power in c(1, 2, 3, 4)) {
  row <- nsca_table(nsca_analysis(
    data.frame(X = grid, Y = grid^power), "X", "Y", ceilings = "ce_fdh"
  ))
  cat(sprintf(
    "Y = X^%d:  d_nec = %.3f  d_suf = %.3f  weakest_effect = %.3f  balanced_joint_effect = %.3f  joint_empty_space_coverage = %.3f\n",
    power, row$d_nec, row$d_suf, row$weakest_effect, row$balanced_joint_effect, row$joint_empty_space_coverage
  ))
}

# The same two components under different degrees of compensation. p = -Inf
# gives the minimum, p = 0 the geometric mean and p = 1 half the sum. The sum
# does not distinguish a balanced pair from a one-sided one; this is the
# objection to it, and it concerns degree rather than kind.
for (p in c(-Inf, -1, 0, 1)) {
  cat(sprintf(
    "p = %5s:  (0.40, 0.40) -> %.3f    (0.70, 0.10) -> %.3f\n",
    format(p), nsca_joint(0.40, 0.40, p = p), nsca_joint(0.70, 0.10, p = p)
  ))
}

# A constant outcome is rejected with an error rather than scored, because it
# has no scope of its own.
flat <- data.frame(X = seq(0, 1, length.out = 60), Y = rep(0.5, 60))
try(nsca_analysis(flat, "X", "Y", ceilings = "ce_fdh"))

# With an imposed scope, the analysis runs, but it issues a warning and is
# flagged as degenerate, because an uninformative line at the center of the
# scope maximizes every index simultaneously. When the same line is moved off
# center, the result changes, although the data are equally uninformative.
nsca_table(suppressWarnings(nsca_analysis(
  flat, "X", "Y", ceilings = "ce_fdh", scope = c(0, 1, 0, 1)
)))[c("d_nec", "d_suf", "weakest_effect", "balanced_joint_effect", "joint_empty_space_coverage", "degenerate")]
low <- data.frame(X = flat$X, Y = rep(0.02, 60))
nsca_table(suppressWarnings(nsca_analysis(
  low, "X", "Y", ceilings = "ce_fdh", scope = c(0, 1, 0, 1)
)))[c("d_nec", "d_suf", "weakest_effect", "balanced_joint_effect", "joint_empty_space_coverage", "degenerate")]

# The other horizontal case: a flat frontier over a filled scope. Neither side
# has an empty space, so each index is at its floor and nothing is flagged.
set.seed(11)
filled <- data.frame(X = runif(400), Y = runif(400))
nsca_table(nsca_analysis(filled, "X", "Y", ceilings = "ce_fdh"))[
  c("d_nec", "d_suf", "weakest_effect", "balanced_joint_effect", "joint_empty_space_coverage", "degenerate")
]

# The corner mapping, for reference.
nsca_direction_map()
