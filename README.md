# NSCA

Necessary and Sufficient Condition Analysis (NSCA) in R.

## What the method is

A condition is **necessary** for an outcome when the outcome does not occur
without it, and **sufficient** when the presence of the condition guarantees
the outcome. A condition that is both necessary and sufficient defines a
dividing line at each level of the outcome: below the line, that level cannot
be reached; above it, that level is guaranteed.

Each half of this claim implies a specific prediction about the data. When
the condition X is plotted against the outcome Y, necessity implies that one
corner of the scatter plot is empty, and sufficiency implies that the
diagonally opposite corner is empty. For the claim "high X is necessary and
sufficient for high Y", there should be no cases with low X and high Y
(upper-left) and none with high X and low Y (lower-right).

NSCA estimates both empty corners from the same data, with the same scope and
the same frontier technique. It measures each corner, tests each against
random pairing of X and Y, and assesses the two jointly. The band between the
two frontiers is the **admissible region**, the set of combinations that
neither claim rules out. The narrower this region, the closer the relationship
is to a strict necessary-and-sufficient relationship.

| Direction | Claim | Necessity corner | Sufficiency corner |
|:---:|---|---|---|
| `HH` | high X is necessary and sufficient for high Y | upper-left | lower-right |
| `LH` | low X is necessary and sufficient for high Y | upper-right | lower-left |
| `HL` | high X is necessary and sufficient for low Y | lower-left | upper-right |
| `LL` | low X is necessary and sufficient for low Y | lower-right | upper-left |

## What it is for

- **Testing claims that a condition is both required and sufficient.** Many
  theories describe a condition as a prerequisite that, once met, also brings
  about the outcome. NSCA tests both components on the same data and reports
  joint support only when each component is supported.
- **Locating the dividing line.** For each outcome level, the dual-threshold
  table reports the X level that the outcome requires, the X level that
  guarantees it, and the interval between them, in which the outcome is
  possible but not guaranteed.
- **Assessing how closely the relationship approximates a
  necessary-and-sufficient one.** The joint table reports the effect size of
  each component and three joint summaries of the two.
- **Complementing regression analysis.** `nsca_compare()` reports the
  regression slope beside both effect sizes, so that a study can indicate
  whether the data display an average effect, a condition effect, or both.
- **Planning a study.** `nsca_power()` estimates the number of observations
  needed to detect both components.

## Installation

```r
install.packages("NSCA")
```

This also installs NCA, which estimates the necessity side, and SCAtools,
which estimates the sufficiency side.

## How to use it

```r
library(NSCA)

# Simulated data lying in a band, so both corners are empty
set.seed(1)
dat <- nsca_random(n = 100, d_nec = 0.3, d_suf = 0.3)

# 1. Fit the claim "high X is necessary and sufficient for high Y"
fit <- nsca_analysis(dat, x = "X", y = "Y", direction = "HH",
                     test.rep = 1000)
fit                     # all three layers, with every setting the result rests on

# 2. Read the results
nsca_results(fit)       # the necessity, sufficiency and joint tables
nsca_table(fit)         # the joint table and the joint-support decision
nsca_thresholds(fit)    # required level, guaranteeing level and the interval
nsca_plot(fit)          # both frontiers and the admissible region

# 3. Set it beside regression
nsca_compare(fit)
```

`nsca()` runs the same analysis, including the permutation tests, and draws
the plot in one step. `nsca_output()` prints the tables and draws the plots
typically needed for a report:

```r
nsca(dat, x = "X", y = "Y", direction = "HH")
nsca_output(fit, thresholds = TRUE, tests = TRUE)
```

### Checking the result

| Function | What it shows |
|---|---|
| `nsca_test_plot()` | both permutation distributions beside the observed effect sizes |
| `nsca_shared_plot()` | the shared permutation distributions of both sides, their minimum and their sum |
| `nsca_scope()` | the scope by which both effect sizes were divided |
| `nsca_scope_sensitivity()` | how both components and the joint summaries change with the scope |
| `nsca_resample()` | how much they vary across subsamples |
| `nsca_outliers()` | which observations determine each frontier |
| `nsca_locationshift()` | whether each frontier is a level shift or a threshold |
| `nsca_magnitude()` | conventional magnitude bands for a component effect size |

## Reading the results

**The two components.** `d_nec` and `d_suf` are the empty areas as
proportions of the scope, ranging from 0 to 1. `d_nec_null` and `d_suf_null`,
reported beside them, are the values that random pairing of the same data
produces. `d_nec_critical` and `d_suf_critical` are the values that random
pairing exceeds only 5% of the time, that is, the critical values of the two
components in a single test. When several conditions are analyzed,
`nsca_table(fit, p.adjust = "holm")` adjusts the p-values for multiple
testing; the critical values remain single-test reference lines.

**Three joint summaries.**

| Column | What it asks |
|---|---|
| `weakest_effect` | how strong the weaker of the two claims is |
| `balanced_joint_effect` | how strong the two are together, with higher values when they are balanced |
| `joint_empty_space_coverage` | how much of the scope the two claims together rule out |

`nsca_joint()` computes the full family of summaries between these, from no
compensation between the components to full compensation.

**The p-values.** Each p-value addresses a different question, and
`joint_support` uses the first:

| Column | What it asks |
|---|---|
| `p_nsca_iut` | is each corner emptier than under random pairing? `max(p_nec, p_suf)` |
| `p_weakest_perm` | is the smaller empty area larger than random pairing produces? |
| `p_sum_perm` | are the two corners together emptier than under random pairing? |

`p_nsca_iut` is significant only when both component tests reject, so
combining them requires no correction. The other two require
`shared.test.rep`, which applies the same shuffles to both sides. A single
strongly empty corner may by itself make `p_sum_perm` small. We suggest
choosing the p-value that addresses the research question before reading the
results.

**Joint support.** `joint_support` is `TRUE` when both components are
significant, both reach the relevance threshold if one was set with
`relevance =`, and the two frontiers have an acceptable geometry. The
individual flags beside it indicate which requirements each row meets.

**Dual thresholds.** For a high-X direction, each outcome level divides the X
axis into three parts:

| X | The outcome level is |
|---|---|
| below `necessity_threshold` | not attainable |
| between the two thresholds | possible but not guaranteed |
| above `sufficiency_threshold` | guaranteed |

For a low-X direction, the inequalities are reversed.
`necessity_sufficiency_interval` is the width of the middle part;
`threshold_status` reports `exact correspondence` where this width is zero.

## Good practice

- **Choice of direction.** We suggest choosing the direction on theoretical
  grounds and then inspecting a plot of the data. Each direction predicts one
  empty corner on each side, which assumes that the relationship is monotone.
- **Reporting the scope.** Every effect size is an area divided by the scope,
  so it is advisable to state whether the scope is the observed range or a
  theoretical one, and on what grounds. `nsca_scope()` and the printed
  `Settings` block record it.
- **Reporting the frontier technique and the number of permutations.** Both
  sides use the same technique; when several techniques are specified, each is
  reported in its own row. At least 1000 permutations are commonly used.
- **Specifying a relevance threshold in advance.** If a minimum effect size
  matters for the conclusion, it can be fixed before the analysis with
  `relevance =`.
- **Calibrating the joint summaries to the design.** Conventional bands are
  available for each component effect size (`nsca_magnitude()`); for the joint
  summaries, `nsca_power()` and `nsca_resample()` indicate which values the
  design can produce.

An ordinary least-squares or median regression line can be added with
`reference = "ols"` or `"qr"`. It is drawn for comparison and reported by
`nsca_reference()`; it is not part of either component.

## Learn more

- `?NSCA` describes the method and lists its references.
- `citation("NSCA")` gives the references to cite.
- `nsca_terms()` lists the vocabulary used in the tables, and
  `nsca_legacy_names()` maps names used by earlier versions to the current
  ones.
- `system.file("examples", "necessary-and-sufficient.R", package = "NSCA")`
  provides a worked example.
- SCAtools, which NSCA uses for the sufficiency side, can also be used on its
  own for Sufficient Condition Analysis.

Source and issue tracker: <https://github.com/youngchanresearcher/NSCA>

The Traditional Chinese version of this README is installed with the package:

```r
file.show(system.file("docs", "README_zh-TW.md", package = "NSCA"))
```
