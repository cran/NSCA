# NSCA 0.5.2

Users upgrading from 0.4.5 should also read *Upgrading from 0.4.5* under 0.5.1
below. This release does not change any value that 0.5.1 reported from
`nsca_analysis()`, except the random-pairing floors noted under *Fixed* and,
in small samples, the shared-sequence p-values noted under the sampling rule.
`nsca_random()` and `nsca_power()` take new arguments and a new default
design, and the tables of `nsca_power()` are now reproducible.

## New: `p_sum_perm`, a test of the two corners combined

`nsca_table()` gains `p_sum_perm`, the permutation p-value of `component_sum =
d_nec + d_suf`. It tests whether the two corners together are emptier than
under random pairing. It is computed from the shared permutation sequence, in
which both components come from the same shuffle, and is `NA` without
`shared.test.rep`.

It addresses a different question from `p_nsca_iut`. A single strongly empty
corner can make the sum large, so data in which only necessity holds may yield
a small `p_sum_perm`. `p_nsca_iut` asks the same question of each corner
separately; it remains the p-value that `joint_support` uses, unchanged.
`joint_empty_space_coverage` differs from `component_sum` only by
`overlap_share`, which is zero unless the two frontiers cross, so the test
applies to it as well.

`nsca_table(p.adjust =)` adjusts it in the same way as `p_weakest_perm`, into
`p_sum_perm_adjusted`. `nsca_shared_plot()` gains a fourth panel for the sum,
and `nsca_output()` prints `p_nsca_iut` beside the joint table, together with
`p_weakest_perm` and `p_sum_perm` when a shared sequence was run.

## New: a critical value for each null distribution

The joint table gains `d_nec_critical`, `d_suf_critical`,
`component_sum_critical` and `weakest_effect_critical`: the `1 - alpha` quantile
of each statistic's null distribution, that is, the value that random pairing
exceeds in a proportion `alpha` of permutations. Each statistic has its own
critical value. The critical value of the sum is not the sum of the two
component critical values, because the two components rarely reach their upper
tails in the same shuffle. The last two columns require a shared sequence.

The component tables of `nsca_results()` and `nsca_compare()` gain
`effect_size_critical`, the same quantity for one component, as reported by
SCAtools 0.5.2. `nsca_test_plot()` and `nsca_shared_plot()` draw each critical
value as a dashed line beside the observed value.

## New: the Monte Carlo margin of each p-value

The joint table gains `p_nec_accuracy`, `p_suf_accuracy`,
`p_weakest_perm_accuracy` and `p_sum_perm_accuracy`: the half-width of an
interval around each p-value at `test.p_confidence`, computed as NCA computes
its `p-accuracy`. With the engines' tests they are the engines' own values;
with a shared sequence they are computed from it, and they are zero when every
ordering of the outcome was evaluated, since the p-value is then exact.

## Changed: the shared sequence draws orderings as NCA does

The shared sequence now follows the rule of NCA's own permutation test:

- with 720 permutations or fewer, every ordering of the outcome is distinct;
- a count larger than the `n!` orderings of a sample of size `n` is lowered to
  `n!`, with a message, so that every ordering is used once and the p-values
  are exact.

The draws come from the same random stream as before, and a repeated ordering
is redrawn rather than evaluated twice. With the same `shared.seed`, results
are therefore identical to those of 0.4.5 and 0.5.1 whenever no ordering would
have repeated; in practice, this applies to samples of about ten observations
or more. In our checks on a sample of 40, at 199 and at 999 permutations,
0.4.5, 0.5.1 and 0.5.2 gave the same p-values to the last digit. Only small
samples, or a count above `n!`, can give different p-values, and these are the
cases the rule is intended to improve. `shared.test.rep` in the model's
settings records the requested count.

## Fixed: the random-pairing floors come from the same draws as the p-values

`d_nec_null` and `d_suf_null` in the joint table were previously read from the
engines' own tests even when the p-values beside them came from the shared
sequence, and were `NA` when only the shared sequence was run. Each null
summary in a row (floor, excess and critical value) is now read from the
distribution from which its p-value was obtained.

## Changed: `print()` states the question each p-value addresses

The decision block labels each p-value with the question it addresses and
states that `joint_support` uses `p_nsca_iut`:

```
p_nsca_iut      each corner is emptier than random pairing: max(p_nec, p_suf)
p_weakest_perm  the smaller empty area is larger than random pairing produces
p_sum_perm      the two corners together are emptier than random pairing
```

When a shared sequence was run, the multiplicity line in `Settings` lists the
supplementary tests it adds beyond the joint decisions, and recommends that
the p-value addressing the research question be chosen before the results
are read.

## Checked: the false-positive rate of each p-value

Simulation setting: independent X and Y, n = 30, 99 shared permutations,
nominal alpha 0.05, 2500 simulations (SE 0.0044):

```
p_nec            reject .0492   95% CI [.0407, .0577]
p_suf            reject .0564   95% CI [.0474, .0654]
p_nsca_iut       reject .0068   95% CI [.0036, .0100]
p_weakest_perm   reject .0580   95% CI [.0488, .0672]
p_sum_perm       reject .0520   95% CI [.0433, .0607]

at least one of p_nsca_iut, p_weakest_perm, p_sum_perm: .0772
```

In these simulations, each single test rejected at a rate whose interval
contains .05, and the intersection-union test was conservative, as expected
from its construction. The last line illustrates why we suggest choosing the
p-value before reading the results: taking whichever of the three is
significant raised the false-positive rate from at most .05 to about .077 in
this setting.

At n = 30, a step-frontier effect size takes few distinct values, so the
component p-values are discrete, and a test of their uniformity reflects the
ties more than the calibration; the rejection rate at alpha is therefore the
more relevant property. The sum, which takes more distinct values, also passed
a uniformity test (KS p = .23).

The script and its result are installed with the package, in
`system.file("simulations", package = "NSCA")`.

## New: the power of the supplementary tests

`nsca_power(shared = TRUE)` draws one shared permutation sequence of
`test.rep` shuffles per replication, takes all five p-values from it, and adds
`power_weakest_perm` and `power_sum_perm`: how often `p_weakest_perm` and
`p_sum_perm` reject. Both columns are `NA` with `shared = FALSE`, the default,
which keeps the engines' separate tests. `nsca_powerplot(show =)` accepts
`"weakest"` and `"sum"`; neither is drawn unless requested.

In our simulations, when both corners were empty, the supplementary tests
rejected more often than the joint test. When only one corner was strongly
empty, `power_sum_perm` remained high and `power_weakest_perm` exceeded the
rate at which the weaker corner's own test rejected, because both are tested
against random pairing of the two sides at once rather than against one
component holding without the other. This is one reason why `joint_support`
uses `p_nsca_iut`; the help page of `nsca_power()` describes this behavior.
The simulation scripts and their results are installed in the same folder as
the calibration.

## Documentation: what the permutation tests compare

`?NSCA` gains a section on what the permutation tests compare. Every p-value
in the package compares the data with random pairing of the same values, so a
small p-value indicates that a corner is emptier than random pairing would
leave it, rather than that it is empty in the population. When one corner is
strongly empty and the other is not, random pairing can fill the second corner
more than the data do, and its test may reject more often than the nominal
level; a plain correlation can have the same effect. The section offers some
ways of taking this into account, presented as considerations rather than
rules. In addition:

- the help page of `nsca_table()` states that `p_nsca_iut` is level-alpha
  against random pairing and only as reliable as the two component tests it
  combines; it previously gave the level without that qualification;
- it also states that, on one shared sequence, `p_weakest_perm` cannot exceed
  the p-value of the component whose observed area is smaller;
- the help page of `nsca_power()` describes how to check a design by
  simulating it with one side absent, how precise the simulated rates are,
  and the simulated design itself: the condition is uniform and, at each value
  of the condition, the outcome is uniform between the two boundaries. It
  previously stated that the band was filled uniformly.

## Fixed: `nsca_power(seed =)` results are now reproducible

The engines' permutation tests reseed the session from the clock when they
finish. As a result, every replication after the first previously drew its
data from an unseeded stream, and two runs with the same `seed` gave different
tables. Each replication now draws a data seed and a test seed from the seeded
stream before anything else runs. The engines' progress lines are no longer
printed during a simulation, and `test.rep = 0`, under which no replication
can reject, is now refused.

## Changed: `nsca_random()` and `nsca_power()` describe the band by its effect sizes

Both functions now take `d_nec` and `d_suf`, the population effect sizes of the
two corners, and `slope_nec` and `slope_suf`, the slopes of the two boundary
lines on the unit square; the intercepts follow. They replace `lower`,
`upper`, `centre` and `width`, which described the same lines by their
intercepts and required the effect sizes to be computed by hand. Neither
function has been on CRAN, so no released code depends on the old arguments.

- The default is 0.10 on each side with slope 1, matching the default effect
  size of `SCAtools::sca_power()`. The previous default band left about 0.38
  on each side.
- `d_suf = 0` or `d_nec = 0` gives a pattern with that corner absent, which is
  how `nsca_power()` can check a design for false positives (see its help
  page).
- Under the old arguments the default `width` was read as a band even when
  `lower` and `upper` were given, so the two-line form was refused unless
  `width = NULL` was passed as well. This no longer occurs.

## New: a warning when a test cannot reach its threshold

A permutation p-value cannot fall below `1 / test.rep` from the engines, or
`1 / (shared.test.rep + 1)` from the shared sequence, nor below the bound set
by the `n!` orderings of a small sample. `nsca_analysis()` now warns when the
floor of the p-values its tables will hold is above `test.p_threshold`, and
`nsca_power()` warns when it is above `alpha`: in that case neither component
could be significant, and every simulated power would be 0 regardless of the
data. The warning names the cause and, where more permutations would help, the
count needed. SCAtools' own warning for the sufficiency side is not repeated.

The examples of `nsca_power()` and `nsca_powerplot()` used `test.rep = 19`,
which cannot reject at 0.05 and therefore showed a power of 0. They now use
`test.rep = 39`.

## Changed: adjustment leaves the critical values unchanged, with a message

`nsca_table(p.adjust =)` re-derives every decision from the adjusted p-values
but leaves the `*_critical` columns unchanged: they are the critical values of
single tests at `alpha`, serving as reference lines rather than as the
threshold that applies to an adjusted test. An effect size that exceeds its
critical value can therefore fail to reach significance after adjustment. The
table now prints a message to this effect whenever it adjusts, and its help
page explains why there is no adjusted counterpart for the step-wise methods.
The help page also no longer advises passing `p_nsca_iut` to
`stats::p.adjust()`; the documentation of the `p.adjust` argument explains why
that is not the appropriate quantity.

## Fixed: the warning for an adjustment that cannot reject

The warning that `nsca_table(p.adjust =)` gives when no condition can reach
alpha after adjustment has been corrected in three respects:

- it took the smallest p-value as `1 / (test.rep + 1)`, whereas the engines'
  smallest p-value is `1 / test.rep`, and therefore recommended one
  replication fewer than needed (159 instead of 160 for eight conditions);
- it used `test.rep` even when the p-values came from the shared sequence; it
  now uses that sequence's own count and floor;
- it stated that no condition could reach alpha under every method. This
  holds when only one condition is strong; under BH, Hochberg and Hommel,
  several equally strong conditions can share a smaller adjusted p-value.

Both cases are now computed with `stats::p.adjust()` itself, and the warning
states which one applies.

## Changed: American English throughout

Help pages, messages, printed output and plot captions now use American
spelling. Two columns of `nsca_power()` are renamed accordingly, to
`d_nec_realized` and `d_suf_realized`; `nsca_power()` has not been on CRAN,
so no released code reads the old names.

## Documentation: wording of help pages and messages

Help pages, warnings, messages and printed notes have been revised for a more
measured tone. Several warning texts change as a result, so code that matches
the text of a warning may need updating. No computed value changes.

## Other changes

- The scripts and results of the simulations described in this release are
  installed in `system.file("simulations", package = "NSCA")`, with a README
  giving the run time of each.
- `Imports: SCAtools (>= 0.5.2)`, for `effect_size_critical`.
- `nsca_terms()` lists "component sum" and "critical value".
- The installed example script `necessary-and-sufficient.R` referred to
  `joint_empty_zone_coverage`, a column renamed in 0.4.0, and stopped with
  "undefined columns selected". It now uses `joint_empty_space_coverage`.
- The error for a `scope` that gives the outcome more than one range no longer
  also raises a warning about an unused `sprintf()` argument.
- The area identity in the help page of `nsca_table()` uses the current name,
  `admissible_region_share`, in place of `data_zone_share`.
- The slowest integration tests call `testthat::skip_on_cran()`, so that the
  check time on CRAN stays within its limit. They run when the environment
  variable `NOT_CRAN` is `"true"`, as under `devtools::test()`. A unit test of
  the adjustment floor, which one of them covered, runs in every setting.

# NSCA 0.5.1

NSCA 0.5.0 was not released on CRAN, and was developed alongside the CRAN
releases 0.4.4 and 0.4.5 rather than after them. This release brings the 0.5.0
changes below to CRAN users of 0.4.5 and includes all fixes made in 0.4.4 and
0.4.5. We suggest reading *Upgrading from 0.4.5* first; it collects in one
place every change that can break a script written for that release.

## Changed: the technique column is `frontier` in every table

The column naming the frontier technique was `frontier` in the component tables
of `nsca_results()` and `ceiling` in the joint tables. 0.5.0 renamed both to
`frontier_technique`. It is now `frontier` everywhere: in `nsca_table()`,
`nsca_thresholds()`, both component tables, `nsca_compare()`,
`nsca_locationshift()`, `nsca_power()`, `nsca_resample()` and
`nsca_scope_sensitivity()`. The component tables therefore keep the name they
had in 0.4.5, and only the joint tables change.

`ceiling` was removed from the joint tables because the same row used "ceiling"
for the boundary type in `necessity_boundary`; `frontier` identifies which
fitted frontier the row reports, with `frontier_accuracy` beside it.
`nsca_table(legacy = TRUE)` and `nsca_thresholds(legacy = TRUE)` still append
`ceiling` as an exact copy. The plot legend reads "Frontier". The arguments
`ceilings` and `ceiling` are unchanged.

SCAtools 0.5.1 makes the same choice, so the sufficiency table it returns and
the necessity table built here share the column.

## Carried over from 0.4.4 and 0.4.5

- `URL` and `BugReports` point to <https://github.com/youngchanresearcher/NSCA>,
  and both READMEs install from CRAN with `install.packages("NSCA")`, as in
  0.4.4 and 0.4.5. 0.5.0 still carried the old repository name and installed a
  0.4.3 source tarball.
- `nsca_extract()` refuses a frontier that was not estimated on every route.
  0.5.0 checked this only on the joint route, so `nec:` and `suf:` returned
  `NA`, which could be read as a missing measurement rather than an invalid
  request.
- `nsca_extract()` has a runnable example instead of one in `\dontrun{}` that
  used a `fit` object it did not create.
- 0.5.0 had independently arrived at the other 0.4.5 fixes: outcome bounds
  taken per condition when a condition has missing values, a named `relevance`
  reordered as numbers regardless of its type, and a seed that leaves no random
  stream behind in a fresh session. The 0.4.5 tests that cover them are now
  part of the suite, to guard against regressions.

## Terminology: Sufficient Condition Analysis

The sufficiency side is named *Sufficient Condition Analysis*. `print()`
already used this name in 0.5.0; the README, its Traditional Chinese
translation and `citation("NSCA")` now do as well.

The READMEs also used two joint-table names retired in 0.5.0,
`joint_empty_zone_coverage` and `joint_support_status`; they now use
`joint_empty_space_coverage` and `significance_status`, and their tables of
retired names list the 0.5.x renames.

## Documentation

- `Description` in DESCRIPTION, the text on the CRAN page, now states what the
  method is, what it is for and how it can be useful. The origins of the
  method, with its authors and references, are described on a new package help
  page, `?NSCA`, which also covers how to use the package, joint support and
  how to read the result.
- The README and its Traditional Chinese translation are rewritten around the
  same questions, with a worked example. The design notes they previously
  contained are now in the help pages, and the version history is kept only in
  this file.
- Every exported function now has a runnable example. Examples that run a
  permutation test use 99 permutations so that they finish within seconds;
  the example of `nsca()`, which is intended to show the default of 1000
  permutations, is wrapped in `\donttest{}`, as in 0.4.5. The example that
  0.4.5 gave `nsca_reference()` had been lost when the help pages were
  regenerated; it has been restored.
- `citation("NSCA")` and the help pages cite the two papers that set out the
  method: Chan (in press), *From Necessity to Sufficiency and
  Necessity-and-Sufficiency: Extending Condition Analysis in Degree*, and Chan
  (in press), *Mean Effects Are Not Sufficiency Effects: How Sufficient
  Condition Analysis Complements Regression Analysis*, both University of
  Leeds.
- `nsca_legacy_names()` stated that `nsca_extract()` accepts retired names with
  a warning and that retired arguments still work. Neither has applied since
  0.5.0. The page now describes the current behavior and lists the 0.5.x
  renames;
  `legacy` in `nsca_table()` and `nsca_thresholds()` is documented as
  appending the names of 0.3.0, 0.4.0 and 0.5.x.
- `nsca_extract()` no longer offers retired names in its documentation or in
  its error message for an unknown parameter.

## Upgrading from 0.4.5

This section lists every change that can make a script written for 0.4.5 fail
or give a different number. Each is explained in the corresponding section of
0.5.0 below.

- In the component tables of `nsca_results()`, `analysis` is now `component`
  and `empty_zone_area` is now `empty_space_area`. These tables have no legacy
  columns.
- In `nsca_table()` and `nsca_thresholds()`, `ceiling` is now `frontier`,
  `joint_empty_zone_coverage` is now `joint_empty_space_coverage`, and
  `joint_support_status` is now `significance_status`. `legacy = TRUE` appends
  the old names as exact copies.
- `nsca_extract()` refuses retired column names. 0.4.5 warned and returned the
  value under the new name.
- `boundary` is removed from `nsca_thresholds()`; `inequality` replaces it.
- Six inputs that were previously accepted are now refused: an unknown
  frontier technique, a `test.p_threshold` outside (0, 1), a negative
  `test.rep`, a reversed `scope`, a `scope` with no width, and the same
  variable as condition and outcome.
- One outcome can no longer be given two different theoretical ranges across
  conditions.
- `Imports: SCAtools (>= 0.5.1)`.

# NSCA 0.5.0

*Not released on CRAN; its changes reached CRAN in 0.5.1.*

Two of the changes below alter reported numbers, and one of them could reverse
a support decision. We suggest reading the *Fixed* sections before upgrading
an ongoing analysis.

## New: sca_magnitude() / nsca_magnitude(), conventional bands with citations

These functions return the bands of Dul (2016) for an NCA effect size, and the
label of the band into which a value falls:

```r
sca_magnitude()                       # the table
sca_magnitude(c(0.04, 0.18, 0.35))    # "small" "medium" "large"
```

`sca_compare()` and `nsca_compare()` gain a `magnitude` column.

**Relation to an earlier position.** The package was already applying Dul's
lower edge as the default `d_cut` and Cohen's as the default `beta_cut`, with
both conventions recorded only in a code comment that users would not see,
while the joint table printed "magnitude benchmarks are not established". In
our view, a convention is better either applied openly or not applied at all
than applied without disclosure.

The bands are a convention to be cited, in the same way that a reader of a
regression cites Cohen. Reporting "d = 0.35, large by the conventional bands
(Dul, 2016)" is an ordinary citation. The table records this in a `source`
column, not only in the help page, because a table may be copied into a
script whereas a help page is not.

### Cohen's bands, and the offset between the two conventions

`beta_cut = 0.30` follows Cohen (1988), and, like Dul's convention, it had been
recorded only in a code comment. Both are now exposed:

```r
sca_magnitude(scale = "beta")     # Cohen's bands for a standardized coefficient
```

`sca_compare()` and `nsca_compare()` gain `beta_magnitude` beside `magnitude`.

The two conventions share their edges, but their labels are offset by one
band. Both use .10, .30 and .50. Cohen's "medium" runs from .30 to .50; Dul's
"large" runs from .30 to .50.

```
value    as a coefficient (Cohen)   as an effect size (Dul)
0.05     below small                small
0.15     small                      medium
0.35     medium                     large
0.60     large                      very large
```

This explains why the two defaults differ: `beta_cut = 0.30` and
`d_cut = 0.10` are each the lower edge of "medium" in their own tradition, and
these are different numbers even though the edges are shared. A reader
comparing a beta and a d in the same row of `sca_compare()` may need this
information, which a comment in the source did not provide.

### The joint indices still have no bands, for a different reason

`weakest_effect`, `balanced_joint_effect` and `joint_empty_space_coverage` are
this package's own constructions. Citing a published convention differs from
proposing a new one for one's own index, and the latter would be a claim
without supporting evidence. `print()` now distinguishes the two cases instead
of stating that no benchmarks exist for either.

### Interpreting a band alongside the floor

A band indicates where a number sits on a conventional scale. It does not
indicate how much of the number the geometry alone would have produced;
`effect_size_null` addresses that. A `d_nec` of 0.10 labeled "medium" at a
floor of 0.055 is a different finding from the same 0.10 at a floor of 0.007,
and both numbers are now reported in the same row of `nsca_compare()`.

## Fixed: the reference line is now labeled

A gray line ran through the middle of the point cloud, and the figure
identified it only in one sentence of the caption. This was reported by a
user reading a plot.

Keeping the line out of the *Component* legend remains appropriate: coloring an
average effect as a component would place it in the same key as an
empty-space frontier, and the two are different kinds of object. However,
leaving it out of that key had resulted in it having no key at all.

It is now labeled on the line itself: `OLS`, or `QR tau=0.5` for a quantile
reference. ggplot2 allows one scale per aesthetic, and both relevant scales are
already in use, so a third legend would require repurposing color or linetype
in a way that reintroduces the conflation the original decision avoided. A
label placed on the line itself avoids both problems.

The labels are staggered along x, because an OLS line and a median line nearly
coincide on symmetric data, and two labels at the same position would overlap,
with one obscuring the other.

This change applies to both packages.

## Checked: calibration of the permutation tests

The Type I error rate was estimated by simulation. Simulation setting:
independent X and Y, n = 50, `test.rep = 199`, nominal alpha 0.05, 2500
simulations (SE 0.00436):

```
necessity      reject .0520   95% CI [.0433, .0607]   KS p = .27
sufficiency    reject .0528   95% CI [.0440, .0616]   KS p = .32
IUT (joint)    reject .0080   95% CI [.0045, .0115]
```

In these simulations, both component p-value distributions were consistent
with uniformity, not only correct at the cut. A first run of 500 simulations
put sufficiency at .0740 (z = +2.46), which suggested a possible
anti-conservative bias; at 2500 simulations it was .0528, consistent with
simulation noise.

The intersection-union test rejected well below alpha and did not pass a
uniformity test. This is expected from its design: it rejects only when both
components reject, so its level is at most alpha by construction.
`joint_support` is therefore conservative, and a null result from it should be
interpreted with this in mind.

The calibration of 0.5.2, installed in
`system.file("simulations", package = "NSCA")`, supersedes this one.

## Fixed: regression frontiers are drawn within the box they are measured in

`cr_fdh` and the other line techniques were drawn across the full x range
without clipping, so near a corner the line crossed the gray dashed scope
boundary and continued beyond it.

The reported areas were not affected: the band integral clips the frontier to
the scope, which is why the identity held and `fit` and `frontier_accuracy`
were 1. However, a figure in which the frontier leaves the box whose area is
reported may appear erroneous to a reader comparing the number with the
picture, and a reader cannot readily verify that the discrepancy is only
cosmetic.

The drawn line is now clipped in the same way as in the integral, so that it
shows the boundary of the region that is measured.

## Fixed: long plot titles are now wrapped within the figure

This was found by rendering and inspecting the figures; the test suite checks
only whether a plot is a `ggplot` and therefore could not detect it.

Title, subtitle and caption were all wrapped at 95 characters. However, a title
is bold and about one and a half times the size of a caption, so with a long
variable name the first line overflowed the panel and was clipped mid-word:

```
NSCA: High Hours of structured professional development training per year
is necessa|
```

In `nsca_shared_plot()`, the title was not wrapped at all and extended beyond
the figure. Titles now wrap at a width chosen for their own size, determined by
rendering rather than by calculation from point sizes.

## Fixed: overlapping power curves

`nsca_powerplot()` drew the joint curve last, on top of the two component
curves from which it is built. At high power, which is the usual case, all
three coincide, so the figure showed one line above a legend naming three.

The joint curve is now drawn first, underneath, with the two component curves
dashed and dotted over it. It also comes first in the legend, because it is
the curve on which a study's power is based. The x axis now breaks at the
sample sizes that were simulated; it previously broke at 50/75/100 for a study
run at n = 30/60/120.

## Changed: the plot legend is titled "Technique"

It was previously titled `Frontier`, which is now the name of the boundary
itself. The tables call the estimator `frontier_technique`, and the figures are
now consistent with them.

## New: the configuration label is reported with its floor

`sca_compare()` and `nsca_compare()` gain `effect_size_null` and
`effect_size_excess`.

In one example, a pure-noise condition at n = 70 was classified as
configuration II, "necessity-dominant", with `d_nec = 0.104` against the
conventional cut of 0.10, while its independence floor was 0.055; about half
of that effect size therefore reflected the envelope fitting sample extremes.
The label was defensible, but the reading could be misleading.

```
condition   component  effect_size  effect_size_null  effect_size_excess  label
Training    necessity      0.36825           0.07451             0.29374  mean and necessity
Luck        necessity      0.10418           0.05464             0.04955  necessity-dominant
Luck        sufficiency    0.03725           0.08698            -0.04974  neither
```

`Luck`'s excess is below the cut on one side and negative on the other. The
cut is unchanged, since applying it to the excess would amount to introducing
a new criterion, but the number that qualifies it is now reported in the same
row.

## Internal: SCAtools' NAMESPACE is again generated by roxygen

It had been maintained by hand, so each new `@export` also had to be added to
it manually, and `sca_resample()` was nearly released without being exported.
Before roxygen was given responsibility for the file, the tags and the file
were checked to be in agreement: 31 exports and 5 S3 methods, with nothing
dropped and nothing missing.

## Fixed: the last remaining error message from the engine is replaced

```
before:  The length of scope needs to be equal to the length of x
         or
after:   'scope' has 6 values. With 3 conditions it must have 4 --
         c(xmin, xmax, ymin, ymax), one box for all of them -- or 12,
         four per condition in the order the conditions were given.
```

A list of the wrong length, or with an entry that is not four numbers, is now
refused in the same way.

## New: a warning when an adjustment cannot reject

This was identified during adversarial testing. In an example with eight
conditions (one with a real effect and seven pure noise) and `test.rep = 50`,
the real condition is supported without adjustment. Under BH, no condition is
supported, and this is not because the evidence was weak.

The limit follows from arithmetic. A permutation p-value cannot fall below
`1 / (test.rep + 1)`, and the most significant of `k` conditions is multiplied
by `k` under both Bonferroni and Benjamini-Hochberg. Therefore no condition can
reach `alpha` unless

```
k / (test.rep + 1) <= alpha,   that is   test.rep >= k / alpha - 1
```

At `alpha = 0.05`, this allows at most 2 conditions with 50 replications, 5
with 100, 10 with 200 and 50 with 1000. When eight conditions are adjusted at
a low replication count and no condition is supported, the result reflects
`test.rep` rather than the data, and the output previously gave no indication
of this.

```
Warning: With 8 conditions and test.rep = 50, the smallest p-value the
permutation test can produce is 0.0196, and 'BH' multiplies the most
significant of 8 by 8. No condition can reach alpha = 0.05 however strong the
evidence, so an absence of support here is arithmetic rather than a finding.
Use test.rep >= 159.
```

This is a warning rather than an error, because the adjusted p-values are
still correct and a caller may want them; however, the absence of a rejection
should not then be read as evidence. No warning is given when nothing is
adjusted or when no test was run.

## Breaking: six previously accepted inputs are now refused

Adversarial testing of 0.5.0 supplied both packages with degenerate data,
extreme scales and arguments of the wrong shape. Of sixteen invalid inputs,
seven were accepted. The six that are now refused are described below; three
of the seven produced a fit whose own description of itself was inaccurate,
which a user could not detect from the output.

### An unknown frontier technique reached the engine

This was the most consequential case.

```r
sca_analysis(data, "X", "Y", ceilings = "banana")   # used to work
```

The name was not validated. `"banana"` was passed to the engine, which fell
back to its own defaults and fitted ce_fdh and cr_fdh, while the metadata still
recorded `"banana"`; `print()` therefore showed `frontiers  banana` above a
table of two techniques that had not been requested. A misspelled `"ce_fhd"`
had the same effect, without any message.

Such names are now refused, and the message lists the available names. A test
checks the list against the engine's own, so an NCA release that adds or
removes a technique causes a test failure rather than going unnoticed.

### A significance threshold outside (0, 1)

`test.p_threshold = 2` was accepted and stored. Every row was returned with
`significant = TRUE`, and the `Settings` block printed `alpha  2` without
comment.

### A negative replication count

`test.rep = -5` was accepted and stored, and printed in the `Settings` block
beside a `p_value` of `NA`.

### A scope reversed or with no width

`scope = c(1, 0, 1, 0)` was accepted. The engine adjusts such a scope without a
message, widening whatever it is given until the data fit, so the analysis was
correct while `sca_scope()` reported `declared_xmin 1, declared_xmax 0`, which
is not a scope that could be reported in a methods section.
`c(0.5, 0.5, 0, 1)` was likewise widened back to the observed range.

`scope = c(0, NA, 0, 1)` is not affected: `NA` means "no theoretical bound on
this side", which is the documented way to declare only part of a box.

### A variable used as both condition and outcome

`sca_analysis(data, "X", "X")` returned `d = 0.53`. A variable plotted against
itself lies on a perfect diagonal, whose empty corner is half the scope
regardless of the data; the result is a plausible-looking number that is an
artifact of the geometry. In NSCA the problem is compounded, because both
frontiers fall on the same diagonal.

### Error messages now name the argument

`frame[x]` raised R's own `undefined columns selected`, which mentions neither
the argument nor the name that was not found. A caller with eight conditions
had to locate the misspelled name by elimination.

```
before:  undefined columns selected
after:   'x' names 'Trainng', which is not a column of the data.
         Available: Training, Experience, Support, Performance.
```

Both packages now give the same message for the same error, so the messages
are consistent across the two packages.

## New: nsca_scope_sensitivity(), and nsca_normalize()

SCAtools had `sca_scope_sensitivity()`; NSCA had `nsca_scope()`, which reports
the box, but no function reporting the effect of changing it. The gap matters
more here: in a one-sided analysis the scope is the denominator of one number,
whereas here it is the denominator of both components, and therefore of every
joint index built from them.

```r
nsca_scope_sensitivity(fit)
```

The direction of the effect may be counterintuitive. Widening the box does not
dilute the empty areas. The frontiers remain where the data are, so every
square unit added beyond the cloud falls on the empty side of one frontier or
the other: the empty areas grow faster than the box, both components rise, and
`admissible_region_share` falls toward 0. In the limit, the data occupy a small
region of a large box, and almost all of that box is excluded by one claim or
the other.

```
widen  d_nec  d_suf  admissible  overlap_share  reconstruction_error  geometry
 0.00   .353   .471        .177         .00008                .00005  TRUE
 0.25   .462   .515        .079         .00003                .05555  FALSE
 0.50   .526   .555        .044         .00002                .12502  FALSE
 1.00   .595   .608        .020         .00001                .22221  FALSE
```

The geometry check also fails early, for a specific reason. This is not
because the frontiers cross: `overlap_share` decreases, to one part in a
hundred thousand. Rather, `reconstruction_error` exceeds its 0.02 tolerance:
outside the observed range of X there are no data to determine a step
frontier, so the engines' extension of it and this package's own
reconstruction diverge, and the area identity no longer holds.
`reconstruction_error` was designed to detect this situation.

A flagged row should therefore be regarded as unverified rather than as
incorrect. Both diagnostic columns are included in the table for this reason:
a report of "geometry failed" that does not state which criterion failed may
be interpreted as whichever failure the reader expected.

No column reports a change from a baseline, for the reason given in
[SCAtools::sca_scope_sensitivity()]: the choice of baseline row rests with the
caller, and the table does not make that choice on the caller's behalf.

`nsca_normalize()` delegates to `SCAtools::sca_normalize()`, so an analysis
does not need to attach the other package to rescale its inputs. Note that the
effect sizes are not invariant to a monotone rescaling of either variable.

## Breaking: the component tables use `component` instead of `analysis`

The same value previously appeared under two names. `nsca_results()$necessity`
reported `analysis = "necessity"`, while `nsca_compare()`,
`nsca_locationshift()` and `nsca_outliers()` used `component`, and the joint
table already carried `weaker_component` and `component_sum` for the same
concept. Three tables and the documentation used one term; two tables used the
other.

`nsca_results()$necessity$component` and `nsca_table()$weaker_component` now
use the same vocabulary, so the two can be joined without a translation step.

No legacy duplicate column is available: that mechanism applies to
`nsca_table()` and `nsca_thresholds()`, and this column belongs to neither.

## New: nsca_table(p.adjust =), adjusted per component family

```r
nsca_table(fit, p.adjust = "BH")
```

This differs from the SCAtools counterpart in one important respect: a
condition here is decided by two component tests combined, so there are two
families rather than one. Each component's p-values are adjusted over the
conditions, and the joint p-value is rebuilt from the adjusted pair:

```
p_nsca_iut_adjusted = max(p_nec_adjusted, p_suf_adjusted)
```

This is not the same as `p.adjust(p_nsca_iut)`. The two give different
numbers, and only the first agrees with the component columns printed beside
it. The intersection-union test rejects only when both components reject, so
the coherent adjusted decision is that both adjusted p-values are significant
at `alpha`, just as the unadjusted decision requires both raw p-values to be.

Controlling each component family at `alpha` keeps the joint decision at a
level of at most `alpha`, since the joint decision is the more conservative of
the two. The IUT was already conservative in this direction, so the adjustment
reduces power but does not increase false rejections.

The argument adds `p_nec_adjusted`, `p_suf_adjusted`, `p_nsca_iut_adjusted`,
`p_weakest_perm_adjusted` and `p_adjust_method`, and re-derives every
significance and support column. The relevance part of `necessity_supported`
and `sufficiency_supported` is unchanged. `"none"` is the default and leaves
the table unchanged.

## New: the effect size produced by random pairing alone, on each side

`nsca_table()` gains `d_nec_null`, `d_suf_null`, `d_nec_excess` and
`d_suf_excess`; the component tables gain `effect_size_null` and
`effect_size_excess`.

Both frontiers are fitted to the sample's own extremes, so both components
tend to be overstated in a small sample. The joint indices are built from the
two, so they inherit both overstatements rather than canceling them.
`weakest_effect` takes the smaller component and is inflated along with it;
`balanced_joint_effect` is a geometric mean of two inflated numbers. Reporting
the two floors beside the two components allows this to be checked.

The floor is the mean of that component's own permutation null distribution:
under independence the true empty area is zero, so the mean of the
distribution the test already computes corresponds to the overstatement at
this `n`, scope and technique. It is obtained from the permutation
distribution rather than from a model.

`*_excess` is a comparison, not a bias-corrected estimator: the bias under
independence is not, in general, the bias under the alternative. A condition
that is pure noise can show a positive `d_nec` and a negative `d_nec_excess`,
and that combination is itself informative. Both are `NA` when nothing was
tested.

## New: nsca_resample(), reporting four quantities without a confidence interval

`nsca_resample()` refits one condition on random subsamples of several sizes
and reports the distribution of both components and both joint indices at
each size.

```
              quantity   m  observed  subsample_median
                 d_nec  36    0.4664            0.5182
                 d_suf  36    0.3470            0.4064
 balanced_joint_effect  36    0.8046            0.9230
                 d_nec  72    0.4664            0.4846
                 d_suf  72    0.3470            0.3712
 balanced_joint_effect  72    0.8046            0.8549
```

The joint indices are included because they are affected as well: in this
example, `balanced_joint_effect` is 0.12 higher at `m = 36` than in the full
sample. Resampling only the components would leave the reader to infer how the
main quantity this package reports behaves.

No confidence interval is reported, for the reasons set out in
[SCAtools::sca_resample()]: the n-out-of-n bootstrap is inconsistent for
boundary estimators, subsampling restores consistency but is biased upward by
construction, and correcting for that requires a convergence rate that this
package does not assert.

The results may be compared with `d_nec_null` and `d_suf_null`, which measure
the same overstatement in a different way: by random pairing at the full sample
size rather than by fitting the observed pairs at a smaller one.

## Breaking: the technique column is `frontier_technique` everywhere

This change addresses two issues. The joint table called the column `ceiling`
and the component tables called it `frontier`, so the same value had two names
depending on the table. In addition, in the joint table the name conflicted
with values in the same row:

```
ceiling   necessity_boundary   sufficiency_boundary
ce_fdh    ceiling              floor
```

The column named `ceiling` held `ce_fdh`, while the column holding the value
`ceiling` had a different name. `ceiling` means "Y bounded from above", which
is what `necessity_boundary` reports and what the glossary defines, so the term
is now reserved for that meaning.

```
nsca_results()$necessity    frontier -> frontier_technique
nsca_results()$sufficiency  frontier -> frontier_technique
nsca_table()                ceiling  -> frontier_technique
nsca_thresholds()           ceiling  -> frontier_technique
nsca_compare()              frontier -> frontier_technique
nsca_locationshift()        frontier -> frontier_technique
nsca_power()                frontier -> frontier_technique
```

After this, `ceiling` means only "Y bounded from above", `frontier` means only
the estimated boundary itself, and `frontier_technique` means the estimator
that produced it. The name also matches the argument it comes from:
`ceilings = "ce_fdh"` produces `frontier_technique = "ce_fdh"`.

The arguments `ceilings` and `ceiling` are unchanged. They are NCA's names,
and this package passes them through.

`ceiling` is listed in [nsca_legacy_names()], so `nsca_table(legacy = TRUE)`
and `nsca_thresholds(legacy = TRUE)` append it as a duplicate column holding
the same values as the new column. `nsca_extract(param = "ceiling")` refuses
the retired name and states what replaced it, as it does for every other
retired column here: a caller who simply receives a value would not learn that
the name has changed. A test covers both behaviors.

`frontier` is not included in that map. It was the component tables' name, and
the component tables have no legacy mechanism; listing it would append a
column to the joint table under a name that table never used. For the same
reason, `empty_zone_area` has been absent from the map since 0.5.0.

`nsca_outliers()` keeps its `ceiling` column. That column is NCA's own output
and has a different meaning: `"X"` marks a removed observation that was a
peer, and the neighboring `scope` column marks one whose removal changed the
scope. For `cols`, `qr`, `c_lp` and `custom`, the engine sets the first
unconditionally, since those techniques have no peers to test against; the
column therefore carries no information there, and a more specific name would
imply more than the engine determines.

## New: nsca_compare() and nsca_locationshift()

SCAtools had both functions; NSCA had neither, so the only way to place the
average effect beside a condition claim was to extract `model$sufficiency` and
run the SCAtools function on it. That approach answers the question for only
one of the two components and returns a table that does not indicate which.

Both new functions return one row per condition and component. There is one
average effect and there are two boundaries, and the interpretation depends on
which boundary is being compared: a flat slope with an empty upper-left corner
and a flat slope with an empty lower-right corner are different findings. A
single row would have to select one without indicating which.

```r
fit <- nsca_analysis(data, x = c("A", "B"), y = "Y", direction = "HH",
                     ceilings = "cr_fdh")
nsca_compare(fit)        # beta, effect_size, configuration, per component
nsca_locationshift(fit)  # ols / median / frontier slopes, per component
```

`beta` and the slope columns repeat down each condition's pair, because they
describe the same line.

`nsca_compare()` classifies each pair into the four configurations of the
condition analysis in degree framework. The numbering follows that framework,
where it is stated for sufficiency; the necessity rows apply the same two
questions to the other corner. This reflects the argument on which NCA was
proposed (a condition can be necessary without raising the expected outcome,
which is configuration II read on the necessity side) rather than a new
typology, but the numerals are borrowed, and the documentation states this.
`configuration_label` names the component, so a `necessity` row reads
`"mean and necessity"`, not `"mean and sufficiency"`.

`d_cut` takes its default from the model's `relevance`, per component, and
accepts `c(necessity = , sufficiency = )`. A name that is not one of the two is
an error rather than a signal to fall back on position, for the same reason
that `relevance` gained this check in this release:
`c(sufficiency = 0.2, nec = 0.1)` would otherwise be applied in reverse without
a warning, and this argument determines a reported classification.

`nsca_locationshift()` reports whether each frontier is the central tendency
shifted as a whole to one edge of the data. The two components can give
different answers, and that asymmetry (a band one edge of which tracks the
regression while the other does not) cannot be seen in a table with one row per
condition. Its two central-tendency slopes are computed by
`SCAtools::sca_locationshift()` rather than here, for the same reason that
`nsca_reference()` delegates: two implementations of a median slope can
diverge over time, whereas a single implementation cannot.

## New: nsca_scales() and nsca_rescale()

`nsca_thresholds()` reports on a scale chosen by `threshold.x` and
`threshold.y`, and the vocabulary for that choice was available only in
SCAtools. A reader who wanted to know what `"percentile", "directional"` meant,
or to place a value of their own on the same scale before comparing it with a
reported threshold, had to be aware of the other package and attach it.

Both functions delegate to SCAtools. They make the conversions available from
NSCA and ensure that the scale of an NSCA threshold and the scale a caller
converts to are the same, because a single function is used rather than two
that currently agree.

## New: nsca_outliers() takes a fitted model

It previously took `data`, `x`, `y` and `direction`, which meant restating
information already contained in the model, with the possibility of restating
it differently. This matters more here than in SCAtools: a different
`direction` could place the two components on corners the fit did not use,
which would invalidate a two-sided table.

```r
nsca_outliers(fit)
```

Every condition is examined unless `x` selects a subset. `y`, `direction` and
`scope` supplied alongside a model are refused, with an error naming the
argument.

The strict column check described under *Changed: nsca_outliers() refuses a
mismatch instead of repairing it* still applies, to the same case: the two
sides of one condition, given identical arguments, are expected to return
identical columns; otherwise the engine's interface has changed. Across
conditions, the model form instead takes the union of the columns, because
there a difference is a property of the fit rather than of the engine: a
declared scope may be given per condition, so one condition can use a declared
box while another uses the observed range, and `NCA::nca_outliers()` drops its
`scope` column for the first and keeps it for the second.

### The table now names its condition, and the class has a print method

Both forms gain a leading `condition` column, so the model form's rows are
distinguishable and the two forms return the same structure.

`print.nsca_outliers()` has been added. The class had been set since 0.4.x
without a corresponding method, so the three attributes it carried (the
direction and the two corners) were not displayed.

## New: print() reports how many tests the table contains

Each condition involves two component tests, combined into one decision by the
intersection-union test, and the output did not previously report how many of
either there were. The `Settings` block now reports both:

```
  multiplicity  3 conditions: 6 component p-values, 3 joint tests
                the intersection-union test is conservative within a
                condition; nothing is adjusted across conditions
                were the 3 joint tests independent, P(at least one
                false positive) at alpha 0.05 would be at most 0.14
```

Both statements are reported together. Within a condition, multiplicity is
already accounted for: the IUT rejects only when both components reject, so
`p_nsca_iut` is the larger of the two and its level is at most `alpha`.
Across conditions, no adjustment is made. Reporting "conservative" without
"nothing is adjusted across conditions" would give only the reassuring part of
the statement, which holds for one condition but not for three.

The package reports the multiplicity rather than adjusting for it, for the
reason given in the SCAtools entry: the appropriate correction depends on the
design, not on the data. Applying `stats::p.adjust()` to `p_nsca_iut` takes
one line once the family is known.

Unlike in SCAtools, a single-condition model still receives a line, because
even then two component p-values underlie one decision, and how they were
combined is relevant. When no permutation test was run, the block is omitted.

## New: test.seed, for reproducible component p-values

SCAtools gained `test.seed` in this release. NSCA had no equivalent, so the
two numbers that determine `p_nsca_iut`, and therefore `joint_support`,
changed from run to run and could not be fixed; `shared.seed` fixes only the
shared sequence, which is off by default.

`nsca_analysis(test.seed = )` fixes both. It uses the seed once, to draw two
component seeds, rather than passing the same number to both engines. Each
engine derives its whole permutation sequence from a single draw from the
caller's stream, so identical starting states give identical shuffles. Seeding
both engines with one number would therefore compute the two component
p-values on a single sequence while `p_source` still reported two independent
engine runs; this would amount to the shared-sequence design without its
supporting machinery or labeling. Two drawn seeds keep the components
independent and make both reproducible from a single number.

`nsca()` does not take it, by design. It is the four-argument convenience
function, and its `Settings` block now reports `test.seed  not fixed (default)`
so that this is visible to the reader.

## Changed: nsca() takes test.seed, and reports when the p-values are not fixed

`nsca()` runs both component tests by default and previously had no way to fix
them, so `p_nec`, `p_suf` and therefore `p_nsca_iut` and `joint_support`
varied between runs. It now takes `test.seed`, which is used once to draw two
component seeds, and `print()` notes beneath the results when the p-values are
not fixed.

## Changed: the correspondence columns are opt-in

`corresponds` and `correspondence_tolerance` appear only when a non-zero
`correspondence.tolerance` was declared. At zero they added no information
(`corresponds` marked exactly the rows that `threshold_status` already labeled
`"exact correspondence"`), while still placing a binary verdict beside a
continuous interval, and readers tend to attend to the verdict. The default
output reports the interval; a verdict is added when a width has been chosen
and can be justified.

## Changed: nsca_outliers() refuses a mismatch instead of repairing it

The two sides are given the same `scope`, `k`, `min.dif` and `condensed`, so
both take the same branch of `NCA::nca_outliers()` and return the same columns.
An earlier draft kept the intersection when they differed; this guarded
against a situation that could not be shown to occur, and did so by dropping a
column without notice. It is now an error naming both column sets: if the two
ever disagree, the engine's interface has changed, and this is reported as an
error, as the engine contract check already does elsewhere.

## New: nsca_random() and nsca_power()

`NCA::nca_random()` covers one corner. A necessary-and-sufficient pattern
requires two, on diagonally opposite corners, so the cloud lies inside a band
rather than above or below one line. Neither engine provides this, so
`nsca_random()` generates it.

The band takes two forms, which are not interchangeable. Two lines, `lower`
and `upper`, each `c(intercept, slope)`, allow the two component effect sizes
to differ, as power planning for an asymmetric claim requires: a relation can
be strongly sufficient and weakly necessary. A `centre` line plus a `width` is
the shorthand for the symmetric case. All data are generated for `"HH"` and
then reflected, because every direction is `"HH"` with one or both axes
mirrored; one generator plus a reflection keeps the four directions exact
mirror images, rather than relying on four implementations that may diverge.

`nsca_power()` reports three rates: `power_necessity`, `power_sufficiency` and
`power_joint`. The joint rate is the one relevant to powering a
necessary-and-sufficient study, and it cannot exceed either component rate,
since the intersection-union test requires both components.

`power_if_independent` is the product of the two component rates, and
`rejection_phi` measures the departure from it: the correlation of the two
rejection indicators across replications, which is zero under independence,
positive when the components tend to reject together, negative when they trade
off, and `NA` when one of them rejected in every replication or in none. An
earlier draft reported `power_joint - power_nec * power_suf` under the name
`dependence_gap`. The difference between the joint rate and a benchmark that
is not expected to hold is not a dependence measure; the 2x2 table of
rejections has a correlation, which is such a measure.

Failed replications are redrawn up to ten times per sample size, and any
shortfall is reported. Redrawing alone could introduce bias: if the
replications that fail are those in which the band is tightest, replacing them
without disclosure would bias every rate. The cap and the reported shortfall
address this. This is the same three-part rule that the shared permutation
sequence already uses.

`d_nec_declared` and `d_suf_declared` are the effect sizes the band implies;
`d_nec_realised` and `d_suf_realised` are the means obtained. They tend to
exceed the declared values: the CE-FDH ceiling is the running maximum of the
sample and sits below the population boundary, so the empty region above it is
larger than the true one, and the same argument applies in reverse to the
floor. In our simulations, at width 0.4 the declared 0.32 was realized as about
0.43 at n = 40 and 0.37 at n = 200.

The implications for power are not asserted. Both the observed statistic and
the permutation null, which the same estimator computes on shuffled data, are
inflated; the extent to which the two cancel is an empirical question that
this function does not address.

## New: qr as a reference line, and reference.tau

This passes through the corresponding SCAtools change. `reference = "qr"`
draws a quantile regression at `reference.tau`, default `0.5`, beside the two
frontiers; `nsca_plot()` draws OLS dotted and `qr` long-dashed and names both
in the caption with the quantile.

`nsca_reference()` now delegates to `SCAtools::sca_reference()` instead of
reimplementing the fit. The sufficiency layer of an NSCA result is an
`sca_result` fitted on exactly these observations, so the line that package
computes is the line this package would compute; a single implementation gives
one answer and one place to maintain it. It gains `line` and `tau` columns and
returns one row per condition and line.

## New: nsca_outliers()

`nsca_outliers()` identifies which single observations support the necessity
frontier and which support the sufficiency frontier. These are usually
different cases (necessity rests on the top-left of the cloud, sufficiency on
the bottom-right), so a single list would obscure which claim is fragile. The
`component` column indicates which.

It handles two behaviors of the engine: the engine returns `NULL` and prints a
message, rather than raising an error, when it finds nothing; and it drops its
`scope` column when a theoretical scope was declared, so the two sides can
return different columns.

It is documented as a fragility check rather than a data-cleaning step. An
observation that moves a frontier behaves as the method expects, since
frontiers are determined by extreme cases; identifying one is not in itself a
reason to remove it.

## New: correspondence.tolerance, and the correspondence columns

`correspondence.tolerance` sets how wide a necessity-sufficiency interval may
be while still counting as threshold correspondence. In a perfectly
necessary-and-sufficient relation the two thresholds coincide at every outcome
level, so the interval is zero; the tolerance states how far from zero is
still close enough. `nsca_thresholds()` gains `corresponds` and
`correspondence_tolerance`, and `threshold_status` gains a
`"near correspondence"` reading.

The default is `0`, which is not a recommendation; the package does not make
one. The width at which a relation ceases to be necessary-and-sufficient is a
substantive claim about the subject matter, and a package-supplied default
would be the kind of unjustified tolerance the framework cautions against. At
`0` nothing changes: `corresponds` marks the same rows that `threshold_status`
already labeled `"exact correspondence"`.

A crossing is now interpreted by magnitude, not by sign, which corrects what
this package previously stated. `threshold_status` has always reported a
negative gap as `"frontiers overlap"`, and the documentation described it as
"mutually incompatible ... a modeling problem". That description holds for an
exact envelope (CE-FDH's two frontiers are a running minimum and a running
maximum over the same observations, so the floor cannot rise above the
ceiling) but not for a fitted line. CR-FDH, CR-VRS, `qr` and `cols` regress a
straight line through the peers, and the narrower the band, the more easily
two separately fitted lines cross inside it. A small crossing under a line
technique is therefore characteristic of a band thinner than the fitting
error, which is what a nearly perfect necessary-and-sufficient relation looks
like; it indicates very close correspondence rather than a misspecified model.

`corresponds` is accordingly `TRUE` whenever `abs(gap)` is within the declared
tolerance, in whichever direction the estimate fell; `FALSE` for a genuine
admissible interval; and `NA` only for a crossing too large to be attributed
to fitting error. `threshold_status` gains the matching band, so a small
crossing is now reported as `"near correspondence"` rather than as an overlap.

This tolerance is kept separate from the floating-point tolerance used to
decide whether two computed numbers are equal. The latter has a numerical
answer; the former does not.

`nsca_thresholds()` accepts it as an override, so the verdict can be reported
at two or three tolerances without refitting; this allows a reader to see how
much of the conclusion depends on the chosen value.

It is a width in the units of the condition, so a single non-zero value is
accepted only when there is one condition; otherwise it is refused, and the
message gives the names to use. Zero is exempt, because zero is zero in every
unit. This is the same rule as for the per-condition scope, for the same
reason: a quantity whose meaning depends on the axis cannot be stated once for
axes that differ.

## New: nsca_powerplot()

`nsca_power()` returns three rejection rates per sample size, which are
difficult to compare at a glance in a table. `nsca_powerplot()` draws them as
curves.

The joint curve is drawn heaviest, because it is the curve on which a study's
power is based: the intersection-union test requires both components, so a
design is adequately powered when the joint curve reaches the target, not when
either component does.

`power_if_independent` is drawn as a dashed reference. It is the joint power
the two tests would have if they rejected independently, which they cannot do,
so the vertical distance between it and the joint curve shows the dependence:
above it, the components tend to reject together and a study sized on the
product is larger than necessary; below it, they trade off. It is omitted from
the legend because it is not a fourth test.

A dotted line marks a conventional target, 0.80 by default, which can be
removed. The caption states the caveat that the curves inherit from
`nsca_power()`: the band is filled uniformly, and the realized effect sizes
tend to exceed the declared ones at small n, so the curves should be read
alongside `d_nec_declared` and `d_nec_realised` rather than alone.

It is the counterpart of `SCAtools::sca_powerplot()`, which draws one curve
because a sufficiency analysis has one test.

## Fixed: a constant condition is refused, as a constant outcome already was

A constant outcome has been refused since 0.4.x: it leaves the scope with no
height, so every empty area would be measured against zero. A constant
condition has the same effect on the other axis and was not checked. It
returned a fit in which `scope_area` was `0`, every effect size was `NaN` and
every p-value `NA`: the same degeneracy was refused on one axis and reported
as a result on the other.

It is now refused, and the message names the alternative: a theoretical
`scope` with distinct `xmin` and `xmax` gives the box a width so that the
areas are defined, and the message notes that the result then depends on where
the constant lies within that scope.

The refusal applies only when every condition is constant, in which case the
analysis cannot report anything. One constant condition among several leaves
the others estimable, and stopping would discard them; such rows keep their
coverage warning and their `NaN`, which is the appropriate result for a single
row but not for the whole call.

This was found during a calibration and stress check of properties that the
test suite does not cover (`stress-0.5.0.R`), which also indicated that the
permutation test is correctly sized.

## Fixed: na.action = "listwise" is now applied before the engine is called

`na.action` was applied in the function that computes the axis bounds, the
coverage and the data the plots draw, and the full data frame was then passed
to the engine. The engine drops incomplete pairs per condition inside its own
loop (`p_create_loop_data()`), so every effect size was fitted on that
condition's own rows regardless of the setting of `na.action`.

`"listwise"` therefore changed the denominators but not the numerators, which
is the inconsistency it is intended to prevent, and the documented behavior
(every condition fitted on the same rows) did not occur. Complete data were
not affected, which is why this went unnoticed.

The rule is now applied to the data frame before it is passed to the engine,
on the selected columns only: a missing value in a column that is not part of
the analysis does not concern the analysis, and dropping its row would reduce
the sample without notice. Listwise deletion that leaves no complete row now
produces an error naming `"pairwise"` as the alternative, rather than an engine
failure.

`"pairwise"`, the default, is unchanged.

## Breaking: one outcome column can no longer be given two theoretical ranges

SCAtools adds the same refusal in this release; the reason is more pressing
here.

Both components are areas over the same rectangle, which places `d_nec` and
`d_suf` on one scale and allows the joint indices to combine them, and the
admissible region between them is measured against that rectangle as well.
With a per-condition `scope` whose Y halves differ, the engines compute each
condition's components over that condition's own box (they read
`scope[[id.x]]`), while `.nsca_band_for()` clips and divides by the box this
package builds from the first segment. The identity

    1 - d_nec - d_suf + overlap_share = admissible_region_share

then no longer holds, `reconstruction_error` reports it, and
`geometry_acceptable` becomes FALSE. The check works as intended, but the
message a user sees is "the area identity does not close", which gives no
indication that the cause is a scope specified in two different ways.

This concerns conditions, not the two components. Necessity and sufficiency
are two corners of one scatter plot: `nsca_analysis()` passes both engines the
same `x`, the same `y` and the same `scope`, and only `corner` differs.
`corner` controls `flip.y`, which reverses the order of the outcome levels and
leaves the box unchanged, so the two components cannot differ in their scope.

## New: nsca_test_plot() draws both component nulls in one figure

SCAtools has had `sca_test_plot()` since 0.2.0. NSCA had no counterpart, so a
model fitted with `test.rep > 0` carried two permutation distributions but
offered no way to plot either: `nsca_output(tests = TRUE)` drew the
sufficiency null through SCAtools and printed the engine's text output for the
necessity null, so the two halves of a two-sided analysis were reported in two
different forms.

`nsca_test_plot()` draws both, necessity above sufficiency, each with the
observed value marked and `p_nsca_iut = max(p_nec, p_suf)` in the subtitle.
`nsca_output(tests = TRUE)` now uses it, and still prints the engine's text
afterward for the p-accuracy, for which the plot has no space.

It is not a substitute for `nsca_shared_plot()`, and the caption states this.
These are the engines' own tests, based on two independent shuffles. Each is a
valid null for its own component, and the intersection-union p-value built
from them is valid because it rejects only when both components reject. Two
independent shuffles cannot, however, provide a null for a statistic that uses
both components at once, such as `min(d_nec, d_suf)`: the pairs behind the two
histograms were not drawn from one joint distribution. That still requires
`shared.test.rep > 0`.

## New: print() separates the researcher's decisions from the package's

`nsca()` takes four arguments and reports two engines, three joint indices, an
admissible region, a permutation test and a geometry check, so almost every
number it prints depends on a choice made by the package. `print()` now ends
with a `Settings` block listing all of these choices, with `(default)` beside
each value the call did not supply, including the two geometry tolerances that
determine `geometry_acceptable` and therefore `joint_support`. A verdict that
depends on a value the reader has not seen is difficult to report.

`nsca()` records its own call for this purpose, since it forwards all of its
arguments to `nsca_analysis()` by name. `nsca_output()` ends with the same
block.

## New: nsca_scope() reports the scope both components were divided by

`d_nec` and `d_suf` are each an empty area over the same scope area, which is
why they are on one scale and why the joint indices can combine them.
`nsca_scope()` returns the declared bounds, the bounds the engines used, the
resulting area and whether the two differ, one row per condition. It is the
counterpart of `SCAtools::sca_scope()`.

`nsca_table()` gains `scope_area` and `scope_adjusted` beside `observations`.
Code that selects columns by position rather than by name will need updating.

## Fixed: the engines' scope widening is recorded, and NCA's global flag removed

The NCA engine widens a declared scope until it contains every observation, so
a declared bound that excludes a data point is replaced by the data point and
the denominator of both components changes while the declared scope does not.
`SCAtools::sca_analysis()` runs first in every NSCA fit and now warns about
this per condition and per axis; NSCA records the same in
`metadata$scope_adjustment`, in the `Settings` block and in `nsca_scope()`,
without repeating the warning.

The engine also writes a suppression flag named `scope_warnings_nca` into the
caller's global environment and does not remove it. SCAtools clears the one its
own engine call creates, but NSCA calls `NCA::nca_analysis()` directly
afterward, which creates it again; this flag is now also cleared, unless it
already existed before the fit.

## Changed: plots wrap their text and no longer obscure the empty corners

* Titles, subtitles and captions are wrapped, and captions are left-aligned.
  ggplot2 wraps none of them, and a right-aligned wrapped caption has a ragged
  left edge.
* The two corner labels are translucent with a thin border. Each lies inside a
  corner that the claim states is empty, which is where a counterexample would
  appear. `annotate = FALSE` removes both.
* With a single frontier technique, no linetype legend is drawn. With one
  technique, the key distinguished nothing, and ggplot drew its sample line in
  black; beside a color legend in which each color denotes one of the two
  components, a black key could be read as a nonexistent third component.

## Documentation

`qr.tau` documents the inversion (`tau` for an upper empty corner, `1 - tau`
for a lower one), that `frontier_accuracy` for `qr` is `qr.tau` by construction
rather than evidence of fit, and that an `NA` `fit` is ordinary for `qr`. The
README now states that the four directions assume a monotone relationship,
which is what makes each expected empty space a single corner; a curvilinear or
U-shaped relationship has no single expected corner on either side, so neither
component represents the intended claim, and the joint indices would combine
answers to questions that were not posed.

## Fixed: threshold-table alignment is checked, and three further cases handled

The necessity and sufficiency tables are given the same outcome levels, but the
engine reverses them for a lower corner and an NSCA model always has one corner
of each, so the two are paired by matching the actual level. That match was not
checked: a failure left NA in every threshold and turned every row into
"not estimable", which cannot be distinguished from a genuinely unattainable
target. It now stops with an error stating how many levels could not be
paired.

Three further cases that previously produced no message:

* `.nsca_actual_column()` has always reported whether it read the engine's
  exact values or fell back to its printed text, rounded to three decimals,
  but no code read that flag. SCAtools warns in the same situation, and the
  issue is more consequential here: "exact correspondence" is decided against
  a tolerance near 1e-8, so at three decimals a genuinely zero interval is
  reported as an admissible one.
* `geometry.max_overlap` and `geometry.max_reconstruction` had no validator,
  unlike `relevance` and `shared.test.rep` beside them. They are stored in one
  `c()`, so a character value coerced the pair to character and the later
  `excess_overlap <= "0.01"` became a string comparison, producing an incorrect
  `geometry_acceptable`; an NA value was instead propagated into
  `joint_support`.
* The frontier-crossing warning compared against a hard-coded 1% while
  `geometry_acceptable` used `geometry.max_overlap`, so widening the tolerance
  produced both a warning and a TRUE for the same row.

## New: `interval_collapsed`

A reported interval of zero did not always mean that the two thresholds
coincide. A comment in the source stated that every scale is a monotone
transformation of X, so that the magnitude is preserved; this holds for the
affine scales but not for `percentile`, which counts observations. Two distinct
thresholds with no observation between them share a rank, so the reported
width collapsed to zero while the actual gap remained positive, and the row
printed `threshold_status = "admissible interval"` beside an interval of 0.

The status is unchanged, since it is decided in actual units, which reflect
the geometry. The new logical column marks the rows in which the scale, not
the data, produced the zero, and a warning names them. An empty selection is
also now refused by `nsca_plot()` rather than passed downstream.

## Fixed: a misspelled `relevance` name could reverse the two components

`.nsca_validate_relevance()` accepted the named form only when both names
were correct:

    all(c("necessity", "sufficiency") %in% names(relevance))

One incorrect name made that `FALSE`, and the function fell back to positional
order without a warning. Thus

    relevance = c(sufficiency = 0.20, nec = 0.10)

assigned 0.20 to necessity and 0.10 to sufficiency, the reverse of what was
intended. A related path led to the same defect: in
`relevance = c(necessity = 0.10)`, `as.numeric()` removed the name, the length
of one triggered the "apply to both" branch, and a threshold stated for one
component was also applied to the other without a warning.

`relevance` determines `necessity_supported`, `sufficiency_supported` and
`joint_support`, so a typographical error could reverse a reported conclusion
without any visible sign. Names are now checked with `setequal()` before
`as.numeric()` removes them, and any names other than exactly the two
component names produce an error.

The existing test could not detect this: it used
`c(necessity = 0, sufficiency = 0.99)`, the canonical order, in which the named
and positional paths give the same answer. Testing the branch requires
reversed order and distinct values.

## Fixed: the permutation denominator no longer decreases without notice

A replication whose engine fits failed was left as `NA` and then dropped by
`.nsca_permutation_p()` from both the numerator and the denominator. A run
that requested 1000 replications could therefore report a p-value computed on
600 without indicating this.

Failed replications are now redrawn, up to ten times the requested count, so
`shared.test.rep = 1000` means a thousand complete replications. A replication
counts only when every requested cell is finite; a partial one is redrawn. If
the cap is reached, the sequence warns with the number obtained and the number
of attempts, because repeated failure is often not random, and a null built
from a biased subset of shuffles may not be the intended null distribution.

## Fixed: `.nsca_with_seed()` restores the stream in a fresh session

The restore was attached only when `.Random.seed` already existed. In a fresh
session it does not exist until the first draw, so `set.seed()` left its state
in the global environment, which the helper is documented as avoiding. When no
seed existed beforehand, the variable is now removed afterward.

## Fixed: an engine parameter that disappears no longer becomes NA

`.nsca_param()` wrapped every lookup in `tryCatch(error = NA_real_)`, which
also caught genuine interface errors along with the ordinary "this technique
has no slope" case. SCAtools now raises a classed condition when the engine no
longer reports a name it requests, and NSCA re-throws it instead of suppressing
it. The necessity side, which reads NCA directly, gained the same contract
check.

## Fixed: an upper-case Y in a condition name altered two columns

The cause was the same two-stage `sub(..., fixed = TRUE)` substitution that
SCAtools contained, copied into `NSCA/R/corners.R`; see the SCAtools 0.5.0
notes. `nsca_corners()` and `nsca_direction_map()` were not affected: they
return the generic wording and do not substitute variable names.

## Breaking

* `boundary` is removed from `nsca_thresholds()`; `inequality` replaces it.
* `nsca_extract()` no longer accepts retired column names. They remain
  available through `legacy = TRUE`, which copies the current column rather
  than translating a name, and therefore gives the same values.
* `Imports: SCAtools (>= 0.5.0)`.
* `DESCRIPTION` still described the region between the two frontiers as "the
  data zone", a term retired in 0.4.0.

# NSCA 0.4.5

## A missing condition value no longer moves that condition's scope

Each condition is fitted by the engines on its own complete pairs, but the
outcome bounds against which NSCA measured the admissible region were taken
from every observed outcome. A row with an outcome and no condition value
therefore widened the box on one side of the identity and not on the other. In
one such dataset, `admissible_region_share` fell from 0.167 to 0.056, the area
identity was off by 0.11, and `geometry_acceptable`, and with it
`joint_support`, became `FALSE`. With two conditions, the effect was not
confined to edge cases: missing values in one condition stretched its bounds to
the outcome range of the other, and the identity was off by 0.65. Outcome
bounds are now taken per condition, from the same pairs and the same scope
segment the engines use. Complete data are unaffected.

## Other fixes

- `relevance` given as a named character vector, such as
  `c(necessity = "0.1", sufficiency = "0.2")`, remained character, and the
  support rule then compared effect sizes with thresholds as strings. Named
  input is now reordered after it has been converted to numbers.
- `shared.seed` restored the caller's random stream only when one already
  existed. In a fresh session it left behind the stream it had created, so every
  later draw followed from that seed. It now removes it, as documented.
- `nsca_extract()` stated that every retired name had been renamed in 0.3.0;
  five of them were renamed in 0.4.0. The warning now names the correct
  release.
- `nsca_extract()` refuses a frontier that was not estimated on every route.
  The `nec:` and `suf:` routes previously returned `NA` for it, which resembles
  a missing measurement; the joint route now gives the same "was not
  estimated" message as `nsca_thresholds()` and `nsca_plot()`.

## Examples

CRAN review requested that `\dontrun{}` be reserved for code that cannot be
run. The one example that used it, `nsca_extract()`, relied on a `fit` object
it did not create; it now builds its own model and runs without a wrapper.

`nsca()`, `nsca_table()`, `nsca_results()` with its `print()` and `summary()`
methods, `nsca_thresholds()` and `nsca_plot()` with its `plot()` method had no
examples and now have runnable ones. The `nsca()` example keeps the default
1000 permutations, which take longer than CRAN's five-second limit, so it is
wrapped in `\donttest{}`; it still runs under `example()` and under
`R CMD check --as-cran`.

`nsca_reference()`'s example and its list of returned columns had been added
to the `.Rd` file by hand and were missing from the roxygen source, so the next
`roxygen2::roxygenise()` would have dropped them. They are now in the source,
and every help page is regenerated from it.

## Documentation

- The area identity in `?nsca_table` still used the name retired in 0.4.0,
  `data_zone_share`; it now reads `admissible_region_share`.
- `legacy` in `nsca_table()` and `nsca_thresholds()` was documented as
  appending the names of one release each; both append those of 0.3.0 and
  0.4.0.
- Both READMEs installed `NSCA_0.4.3.tar.gz` and described SCAtools as a
  source-only package. SCAtools is on CRAN, so installation is now
  `install.packages("NSCA")`, with the source-tarball route kept below it.

# NSCA 0.4.4

`URL`, `BugReports` and both READMEs pointed to
`youngchanresearcher/NSCA_condition_analysis`. They now point to the package's
repository, <https://github.com/youngchanresearcher/NSCA>. Nothing else changed;
the entry for 0.4.3 below had been relabeled 0.4.4 in that release and is
restored to its own version here.

# NSCA 0.4.3

## Authorship and maintainer

`Authors@R` named a placeholder, `NSCA Contributors <maintainer@example.org>`,
which CRAN does not accept: the maintainer is required to be a real person with
a working address. Both are now Young Chan <youngchanresearcher@gmail.com>.

The package had no `inst/CITATION`, so `citation("NSCA")` fell back to the
entry R generates from `DESCRIPTION`. That entry was correct but did not
identify the sources of the reported quantities. A CITATION file now names
NSCA, SCAtools and the NCA method article together, and reads the version from
`DESCRIPTION`.

## `region` was retired but not reachable as a legacy name

0.3.0 renamed the threshold-table column `region` to `threshold_status`, and
that rename is printed in the table below, in both READMEs and in the package
guide. The entry was not added to the legacy map, so `nsca_legacy_names()`
omitted it, `nsca_thresholds(legacy = TRUE)` did not append a `region` column,
and `nsca_extract("region")` reported it as an unknown parameter rather than a
retired name, contrary to the retirement policy. It is now in the map, and the
test now checks the whole documented mapping rather than a sample of it, so
that a missing entry would be detected.

## Documentation

Both READMEs gained an installation section. NSCA requires NCA and SCAtools to
be installed first, and the order matters, so this is stated at the start of
each README.

`DESCRIPTION` gained `URL` and `BugReports`, and both READMEs now link to the
repository. SCAtools is released separately, so its own link is given as well:

    URL: https://github.com/youngchanresearcher/NSCA_condition_analysis
    BugReports: https://github.com/youngchanresearcher/NSCA_condition_analysis/issues

## The interval and the overlap now follow `scale` (breaking)

`nsca_thresholds()` reported `necessity_threshold` and `sufficiency_threshold`
on the requested `scale`, but `necessity_sufficiency_interval` and
`overlap_width` were left in actual units on every scale. The three columns
print side by side, so a reader who subtracted the two thresholds did not
obtain the interval: on the default `"percentage.range"` the two thresholds
were percentages of the range and the width between them was in the units of
X. The behavior was documented (the help page stated that the interval is
"expressed in the units of X" and advised the reader to normalize), but the
package did not produce the normalized number, and the table did not indicate
that the row mixed two kinds of column.

Both columns are now reported on `scale`, so a row can be checked by
subtraction. `necessity_sufficiency_interval_actual` and `overlap_width_actual`
are new and hold the previous values, regardless of the requested `scale`; the
signed `threshold_gap_actual` is unchanged.

Only the reported magnitude changed. The sign is still determined in actual
units and carried over, because under `convention = "directional"` a low-X
direction is mirrored: taking the sign from the reported difference would swap
`"admissible interval"` and `"frontiers overlap"` in precisely the views where
the mirroring applies. `threshold_status`, `thresholds_compatible` and every
component and joint index are unchanged.

Under `scale = "percentile"` the column is a difference of ranks, so it is
read as the share of cases between the two thresholds rather than as a width.
The help page now states this and refers to `"percentage.range"` or `"sd"` for
a normalized width.

The previous values are available in the `_actual` columns or with
`scale = "actual"`.

# NSCA 0.4.2

This release contains three fixes identified by the first full
`R CMD check --as-cran` of the 0.4.x line. All three predate 0.4.1. No code,
argument, column, default or reported quantity changed.

## Corrected: malformed `nsca_table()` documentation (WARNING)

`checkRd` rejected `nsca_table.Rd`: `\emph` is not valid inside `\code`. The
formula for `balanced_joint_effect` had been written to the `.Rd` as
`\code{2 \emph{ sqrt(d_nec } d_suf)}`; a markdown converter had read the two
multiplication signs in `2 * sqrt(d_nec * d_suf)` as emphasis markers and
turned the text between them into `\emph`. The `.Rd` now contains the intended
formula. The roxygen source was already correct, so nothing in `R/results.R`
needed changing.

## Corrected: a test passed an argument its expectation does not accept (ERROR)

`test-invariants.R` passed `info = name` to `expect_gte()` and `expect_lte()`.
The ordered comparisons take `label` and `expected.label`, not `info`, so
current testthat raises "unused argument", which aborts the whole
`test_that()` block rather than failing one expectation. This explains why the
suite reported 1 error against 1032 passes: the assertion itself was correct
but had not been run. Both calls now pass `label`.

The neighboring `expect_true(..., info = ...)` calls are correct and remain
unchanged; `expect_true()` does take `info`.

## `README_zh-TW.md` moves into `inst/docs/` (NOTE)

`R CMD check` counts any unrecognized top-level file as a NOTE. From `inst/`,
the translation is still included in the package and can now be opened from R:
`file.show(system.file("docs", "README_zh-TW.md", package = "NSCA"))`.

## Not addressed: a NOTE unrelated to the package

The "unable to verify current time" NOTE arises when `R CMD check` cannot
reach a time server. It does not concern the package, and no change to the
package can remove it. Setting `_R_CHECK_SYSTEM_CLOCK_ = 0` before checking
skips it.

# NSCA 0.4.1

## "biconditional" is prose; `biconditional` stays a column name

The article calls the joint claim *necessary-and-sufficient*, and uses
"biconditional" only in contrast: section 2.1 says each component claim "is a
conditional, not a biconditional". Using the same word for the joint claim
would therefore be inconsistent with the article on which the package is
based. All prose now uses necessary-and-sufficient: the DESCRIPTION, the
roxygen documentation, the code comments, the printed inference notes, and
both READMEs.

No identifiers were renamed. The `biconditional` column of
`nsca_direction_map()` keeps its name: it holds the joint statement written
out, and the short logical word is the more usable column name.
`necessary_and_sufficient_statement` would be three times as long for the same
content, in a data frame whose other columns are one or two words. Every field
name, list element and status value is unchanged, including the
`necessary_and_sufficient` element of `nsca_results()`.

## OLS is available as a reference line

`"ols"` is an ordinary least-squares regression of the outcome on the
condition, fitted to all the data. It has always been refused as a component,
and still is: it summarizes central tendency, not an empty space, so it can be
neither the necessity nor the sufficiency side. Since 0.4.1 it can be drawn
alongside them. Section 6.1 of the framework presents average-effect analysis
and condition analysis side by side rather than merging them into one number,
and cites the combined importance-performance map as a precedent for showing
two logics in one figure without combining them.

- `nsca_analysis()` and `nsca()` gain `reference`. `reference = "ols"` draws
  the regression line on `nsca_plot()` as a dotted gray line, outside the
  component legend, with the caption identifying it.
- `ceilings` still refuses `"ols"`; the warning now indicates where the line
  can be requested instead.
- New `nsca_reference()` reports its intercept, slope, R-squared and the
  observations it was fitted to. It is available whether or not the line was
  drawn.
- `nsca_terms()` gains a "central tendency" row.
- No joint index, threshold, p value or support decision is affected. The
  reference line is reported separately from `nsca_table()` for that reason.

Requires SCAtools >= 0.4.1.

# NSCA 0.4.0

This release adopts the vocabulary of *condition analysis in degree*, the
framework that places SCA and NSCA beside NCA. Terminology follows the article;
the package's own names are retained only where the article does not provide
a term. Nothing was removed: every retired name still works.

## The admissible region and the necessity-sufficiency interval

The framework names the two quantities 0.3.0 had called the data zone.

- The part of the scope compatible with both components is the **admissible
  region**, so `data_zone_share` becomes `admissible_region_share`.
- For a fixed outcome target, the distance on the X scale between the two
  thresholds is the **necessity-sufficiency interval**, so `data_zone_width`
  becomes `necessity_sufficiency_interval`.

Two `threshold_status` values change with them: `"determinate"` becomes
`"exact correspondence"`, the limiting case in which one X value is both the
minimum required and the sufficient level, and `"data zone"` becomes
`"admissible interval"`.

The documentation now states three points that the framework addresses and
that 0.3.0 did not. The interval is expressed in the units of X, so it is
scale-dependent and comparable across studies only after normalization. It is
not a confidence interval and carries no statistical uncertainty. Finally, an
estimated interval of zero is consistent with exact correspondence but does
not by itself establish it; a small non-zero interval should be called near
correspondence only when a substantively justified tolerance has been defined.

## Joint support

Evidence for both components is **joint support**, which is a claim about the
content of the joint hypothesis rather than a requirement that the components
have identical effects or thresholds.

- `nsca_supported` becomes `joint_support`.
- `conjunction` becomes `joint_support_status`.

## The disjointness check

Section 5.2 of the framework states a consistency check that the package did
not previously implement. Whenever the two frontiers do not cross, the floor
lies at or below the ceiling throughout the scope, so no point belongs to both
expected empty spaces; because each component effect size is that space's
share of the same scope, the two shares cannot sum to more than one. Component
effect sizes summing to more than one therefore indicate that the frontiers
have crossed and that the result requires diagnosis, not that it is unusually
strong.

- `component_sum` reports `d_nec + d_suf`, so the framework's own rule can be
  applied directly by any reader.
- `frontiers_crossed` answers the same question from the measurement rather
  than from the sum.

The sum serves as a proxy when a crossing cannot otherwise be measured. Here it
can be measured: `overlap_share` is the area on which the two empty spaces
overlap, integrated from the reconstructed frontiers. Deciding from the
measurement avoids a limitation of the proxy: the sum comes from the engines
and the comparison threshold from the geometry, so the two can disagree in the
last few digits in precisely the configuration that the check is intended to
accept.

That configuration is the perfect diagonal. A step frontier overlaps itself by
construction: between two consecutive observations the upper staircase still
holds the earlier value while the lower one has already moved, so the two
spaces overlap on every tread, exactly `1 / (n - 1)` of the scope summed over
the `n - 1` treads, and a perfect diagonal attains it. That allowance is
subtracted before the comparison. For a smooth frontier the allowance is zero,
and any overlap is a genuine crossing.

`frontiers_crossed` asks whether the frontiers crossed at all;
`geometry_acceptable` asks whether they crossed by more than
`geometry.max_overlap`. A row can be flagged as crossed and still have
acceptable geometry. The distinction is useful: the first indicates that a
crossing occurred, and the second whether it exceeds the tolerance.

## Ceiling, floor, and the empty space

The framework distinguishes a **boundary**, the theoretical line separating an
expected empty space from the compatible region, from a **frontier**, that
boundary estimated from data, and names the direction in which Y is bounded a
**ceiling** or a **floor**. The package estimated the appropriate quantities
but did not state which was which, because the argument is called `ceilings`
in every direction.

`nsca_table()`, `nsca_results()`, `nsca_corners()` and `nsca_direction_map()`
now report, in the wording of Table 1:

- `relationship` -- the direction as a claim about the theorized X-Y
  relationship, with the user's own variable names, for example
  `"Higher training for higher score"`. This is a claim about the relationship,
  not about the shape of the estimated frontier, which may be a step function
  or a straight line under any of the four types.
- `necessity_empty_space` and `sufficiency_empty_space` -- what each expected
  empty space would contain.
- `necessity_boundary` and `sufficiency_boundary` -- `"ceiling"` or `"floor"`.
  The two are always one of each, because the corners are diagonally opposite.

`nsca_thresholds()` gains `relationship` and records the `inequality` it used.

## Boundary is no longer an argument name

`nsca_thresholds()` renames its `boundary` argument to `inequality`, matching
SCAtools 0.4.0. A boundary is the theoretical line by which an empty space is
separated; it is unrelated to whether a reported rule uses `>` or `>=`. The old
argument still works, with a warning.

## Retired names still work

- `nsca_legacy_names()` now covers the 0.4.0 retirements as well as the 0.3.0
  ones.
- `nsca_thresholds()` gains `legacy = TRUE`, which `nsca_table()` already had.
  Each table appends only the retired names whose replacement it carries.
- When asked for a dual-threshold column, `nsca_extract()` now explains that
  such a column is defined per outcome level and refers to `nsca_thresholds()`,
  instead of reporting a recognized name as unknown.

## Other changes

- `nsca_terms()` is a new machine-readable glossary linking each term of the
  framework to the name this package uses for it.
- Requires SCAtools >= 0.4.0, whose `sufficiency_threshold` columns this
  release reads. The two packages now use one name for one quantity.

# NSCA 0.3.0

This release adopts the unified conceptual specification: the reported
quantities now follow one geometric identity instead of being defined
separately, the vocabulary names the geometry rather than an interpretation,
and a statistic of both components can be tested appropriately for the first
time.

## The identity on which all area quantities are based

The two empty zones and the region between them tile the analytic scope, so by
inclusion and exclusion

    data_zone_share = 1 - d_nec - d_suf + overlap_share

Every other area quantity is derived from this identity rather than defined
separately. Two consequences are new:

- The bounds are exact and tight. `data_zone_share >= 0` gives
  `d_nec + d_suf <= 1 + O`, hence `weakest_effect <= 0.5 + O/2` and
  `balanced_joint_effect <= 1 + O`. Both are attained exactly by a perfect
  diagonal, where `O` is `1/(n - 1)`: at n = 40 the observed values are 0.5128
  and 1.0256 against bounds of 0.5128 and 1.0256. Earlier versions stated only
  that the excess was "of order 1/n"; it is now an equality, so a reader can
  subtract the discretization term instead of approximating it.
- Only the coverage is bounded unconditionally. `data_zone_share` is a share
  of the scope, so it and `joint_empty_zone_coverage` always lie in `[0, 1]`;
  the other two indices are bounded only up to the overlap. This provides a
  technical reason to report the coverage next to them, in addition to the fact
  that it answers a different question.

## New: `reconstruction_error`

The identity has two computational routes. `data_zone_share` and
`overlap_share` are integrated from the frontiers this package reconstructs;
`d_nec` and `d_suf` are reported by the engines from their own internal
frontiers. The residual between the routes is now reported per row.

It provides a direct check that the reconstruction matches what the engine
fitted, and it detects a failure mode that a plot cannot reveal: an incorrectly
rebuilt frontier still looks like a frontier and still yields plausible areas,
but the identity does not hold. This was previously checked on fixtures in the
test suite; it is now computed on the user's own data.

## Renamed columns

Old names remain available: `nsca_table(legacy = TRUE)` appends them as
duplicates, `nsca_extract()` accepts them with a warning, and
`nsca_legacy_names()` returns the mapping.

| retired | current |
|---|---|
| `d_nsca` | `weakest_effect` |
| `j_nsca` | `balanced_joint_effect` |
| `determinacy` | `joint_empty_zone_coverage` |
| `undetermined_share` | `data_zone_share` |
| `weaker_side` | `weaker_component` |
| `p_nsca` | `p_nsca_iut` |
| `undetermined_width` (thresholds) | `data_zone_width` |
| `region` (thresholds) | `threshold_status` |

`determinacy` was retired because it named an interpretation that the number
does not support: a constant outcome scores 1 on it. The region between the
frontiers is the **data zone**, the region in which observations are still
admitted; this does not imply that any were observed there, and the
documentation now states this explicitly. The `threshold_status` value
`"undetermined"` is correspondingly now `"data zone"`.

## New: `p_weakest_perm` and the shared permutation sequence

`shared.test.rep` replaces the engines' two separate permutation runs with one
sequence driving both. Each replication shuffles the outcome once, refits both
engines on that same shuffled frame, and records `d_nec`, `d_suf` and their
minimum together.

This is what makes a p-value for `min(d_nec, d_suf)` meaningful. An
intersection-union rule combines two p-values and requires no assumption about
how they covary, but the minimum is a statistic of both components, and its
null distribution depends on that covariation. Two independent permutation
runs would discard the constraint that the two empty zones are disjoint subsets
of one scope, and would produce a null distribution that no dataset could have
generated.

When the shared sequence is used, `p_nec` and `p_suf` come from it too, so all
four p-values are mutually consistent; `p_source` reports which regime produced
them. The p-values carry the Phipson-Smyth correction, so none can be zero.

`p_weakest_perm` does not replace `p_nsca_iut`, and the documentation now
gives a more precise reason. Its null is random pairing alone, which is a
proper subset of the union null "not necessary or not sufficient": that union
contains states, such as necessary-but-not-sufficient, that do not correspond
to random pairing. Calibrating against random pairing therefore provides no
level-alpha guarantee over the union. It is a sensitivity analysis; the
intersection-union rule remains the decision rule.

## New: the support decision

`necessity_supported`, `sufficiency_supported`, `geometry_acceptable` and
`nsca_supported` are reported separately, so a failure says which of the three
requirements failed rather than collapsing into one boolean.

- `relevance` takes pre-specified practical-relevance thresholds for the
  component effect sizes. Support then needs an effect at least that large as
  well as a significant test. Without it, support rests on significance alone,
  and `print()` states this explicitly.
- `geometry_acceptable` requires a non-degenerate axis, a small
  `reconstruction_error`, and frontier overlap no larger than the `1/(n - 1)`
  a step frontier produces by construction. A degenerate row therefore cannot
  be reported as supported.

## New: `balance`

`balanced_joint_effect` penalizes asymmetry but does not measure it: `(0.50,
0.125)` and `(0.25, 0.25)` both score 0.50. `balance` is the geometric mean
over the arithmetic mean, which is 1 exactly when the two components are equal.
The name of the index implied a property that only this column reports.

## Documentation is now generated from one source

Earlier versions carried hand-written `.Rd` files alongside roxygen comments.
Running `devtools::document()` replaced the `.Rd` files with roxygen output
without warning, and their content was lost; without
`Roxygen: list(markdown = TRUE)`, the result also rendered literal backticks
and `[function()]` links on every help page. The roxygen blocks are now the
single source, the shipped `.Rd` files are generated from them and carry the
roxygen2 header, and markdown is enabled.

## Other changes

- `print()` separates description from decision: section 3 reports the three
  joint summaries side by side, section 4 the p-values and support flags.
- The `nsca_analysis()` documentation gains sections on matched estimation,
  degenerate scopes and the shared sequence.
- Tests: `test-decision.R` and `test-shared.R` are new. The four tests
  distinguishing the joint indices by curvature, low-end inflation and the
  constant-outcome blind spot are restored.

# NSCA 0.2.0

- Adds `j_nsca = 2 * sqrt(d_nec * d_suf)` to `nsca_table()`, the balanced joint
  strength. It penalizes asymmetry rather than only tracking the weaker side:
  at an equal component sum of 0.80, `(0.40, 0.40)` gives 0.80 while
  `(0.70, 0.10)` gives 0.529.

  The factor of two is a normalization, not a new assumption. The two empty
  zones are disjoint subsets of one scope box, so `d_nec + d_suf <= 1`, and by
  AM-GM `sqrt(d_nec * d_suf) <= (d_nec + d_suf)/2 <= 0.5`. The index therefore
  spans `[0, 1]` and reaches 1 exactly when both components equal 0.5, which is
  the single-line case where the two frontiers coincide and the band closes. It
  is left unclamped so that the step frontier's `1/n` self-overlap remains
  visible.

- Adds `nsca_joint()`, which exposes the whole power-mean family that `min`,
  the geometric mean and the sum belong to, at `p = -Inf`, `p = 0` and `p = 1`.

  This revises how earlier versions framed the choice. `d_nsca` was described
  as though the minimum were the uniquely correct joint effect size and the
  alternatives were errors. The alternatives are members of the same family at
  different degrees of compensation, and the geometry of two empty zones does
  not determine the degree. The documentation now names the parameter instead
  of leaving the choice implicit in a formula.

- The geometric mean is documented as partially compensatory, not as a
  non-compensatory logical AND. A larger component does raise `j_nsca`, at a
  diminishing rate. `min` is the fully non-compensatory end of the family and
  the sum the fully compensatory end; `j_nsca` lies between them, and
  describing it otherwise would overstate the properties of a geometric mean.

- No magnitude benchmarks are asserted for `d_nsca`, `j_nsca` or
  `determinacy`, and the documentation now explains why. Conventions for a
  single NCA effect size do not transfer: the scales differ, and under
  independence both components carry a positive finite-sample bias that
  `j_nsca` amplifies rather than cancels, by an amount depending on `n`, the
  frontier technique and the scope. Calibration would require simulation.

- Adds a `degenerate` column and a warning at analysis time for an axis that
  does not vary inside the declared scope. This addresses a limitation shared
  by all three area summaries: a constant outcome at the center of an imposed
  scope leaves the upper and lower halves both empty, so
  `d_nec = d_suf = 0.5`, the band closes, and every index reports a perfect
  biconditional for data that contain no information. `j_nsca` is the most
  affected of the three, since `d_nsca` at least follows the weaker side when
  the constant lies off center. A second, milder warning is issued when the
  observations span less than five percent of the declared scope on either
  axis.

  The other case described as a horizontal line, a flat frontier over a filled
  scope, requires no special rule and has none. There is no empty zone on
  either side, so every index is correctly zero.

- `print()` and `summary()` now show all three summaries side by side with a
  one-line statement of what each measures, and state that no benchmarks exist.

# NSCA 0.1.2

- Adds a `determinacy` column to `nsca_table()`, equal to
  `1 - undetermined_share`: the share of the scope that the two claims together
  rule out. It is the complement of the band rather than the raw sum
  `d_nec + d_suf`, which measures the same thing but inherits the step
  frontier's own overlap and can exceed 1.

  It was added because `d_nsca` and the band answer different questions, and
  neither dominates the other. `d_nsca` is reduced by curvature: a perfect
  bijection such as `Y = X^3` scores 0.25 because the sufficiency zone it
  leaves is genuinely small. `determinacy` is curvature-neutral and reaches 1
  for `Y = X` and `Y = X^3` alike, but it is inflated at the low end, where a
  loose scatter can already score about 0.2.

  `determinacy` is not an effect size, and the documentation states this. A
  constant outcome closes the band completely and scores 1 while satisfying
  neither component. It is the permutation screen, not the index, that rules
  this case out.

# NSCA 0.1.1

- `nsca_results()` and the print/summary methods now expose all three requested
  layers explicitly: necessity, sufficiency, and the joint NSCA table.
- The joint table now reports `p_nsca = max(p_nec, p_suf)` and matching
  component/joint significance flags. The documentation now states the scope
  of this conjunction: the underlying NCA p-values test a
  random-pairing/independence null and do not directly establish
  biconditionality.
- Low-X (`LH` and `LL`) threshold gaps now use the reversed condition-axis
  ordering. Negative gaps are reported as frontier overlap instead of being set
  to zero without notice.
- `.nsca_band()` now returns the documented `crossed` flag, which resolves a
  failure in the package's own geometry test.
- The multi-condition regression test now uses an independent `X2`. The
  earlier visual test used `sort(runif(n))`, which made `X2` monotonically
  related to the already ordered outcome and invalidated the expected wide-band
  check.
- Declaring the ggplot2 data-mask variables removes the R CMD check NOTE.

# NSCA 0.1.0

First release.

- `nsca_analysis()` estimates both components of a biconditional statement on
  diagonally opposite corners of the same scatter plot: necessity through
  `NCA`, sufficiency through `SCAtools`.
- `nsca_table()` reports both component effect sizes, the descriptive joint index,
  the undetermined share of the scope, and an intersection-union verdict.
- `nsca_thresholds()` gives dual thresholds and the undetermined region for
  every outcome level.
- `nsca_plot()` draws both frontiers and shades the undetermined band.
- `nsca_corners()` and `nsca_direction_map()` expose the corner mapping.

## Design decisions

- The descriptive joint index is the weaker component.
  `min(d_nec, d_suf)` falls to zero as soon as either component does and stays
  on the components' own scale. The undetermined band remains separate because
  it answers how much of the scope neither claim rules out.
- The 0.5 ceiling is approximate for step frontiers. The two empty zones are
  disjoint only where the frontiers do not invert. A step frontier overlaps
  itself by about one tread between consecutive observations, inflating both
  components by a term of order `1/n`: a perfect diagonal through 40 points
  gives 0.513. `overlap_share` reports the overlap, and the frontier-inversion
  warning is issued only for smooth frontiers, where an inversion is
  informative.
- A constant outcome is refused, not scored. It has no scope of its own, so
  every empty area would be divided by zero height. With a theoretical scope
  imposed, the analysis runs, but the result depends on where the constant lies
  within that scope, and the error message states this.
- Both sides use the same frontier technique. The descriptive joint index is a
  comparison, and envelopment and regression frontiers do not measure the same
  thing. Mixing them across sides would allow the technique to determine which
  side appears weaker. Several techniques may be requested; each produces its
  own internally matched row.
- Plots are exempt from this rule. A chart makes no numerical comparison, so
  `nsca_plot()` draws every requested technique on both sides.
- The conjunction requires both directional screens. Their maximum p-value
  passes only if both reject. Because the underlying NCA permutation null is
  random pairing/independence, this does not directly establish logical or
  causal biconditionality.
- The engine's purity metrics are off by default. The engine computes them only
  for a corner with neither axis flipped. Of the two corners an NSCA model
  uses, at most one qualifies, and which one depends on the direction: `HH`
  gets them on the necessity side, `LL` on the sufficiency side, `LH` and `HL`
  on neither. They are also expensive to compute when many observations lie on
  the frontier, which is precisely the near-deterministic case that NSCA
  addresses. No NSCA statistic uses them. They can be computed with
  `purity = TRUE`.
- `ols` is refused. It estimates central tendency, not an empty space, so it
  yields no component of either claim.
