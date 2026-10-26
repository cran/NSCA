#' @section The method:
#' A condition is necessary for an outcome when the outcome does not occur
#' without it, and sufficient when the condition guarantees the outcome. In a
#' scatter plot of condition against outcome, each claim rules out one corner.
#' Necessary Condition Analysis (NCA; Dul, 2016) estimates the empty corner of
#' a necessity claim, and Sufficient Condition Analysis, implemented in the
#' SCAtools package, estimates the empty corner of a sufficiency claim. A
#' necessary-and-sufficient claim rules out both corners, which lie diagonally
#' opposite each other, so NSCA estimates both on the same observations, scope
#' and frontier technique and assesses them together. The band between the two
#' frontiers is the admissible region: the narrower it is, the closer the
#' relationship comes to a strict necessary-and-sufficient one. The framework,
#' and the vocabulary the tables use, follow Chan (in press), *From Necessity to
#' Sufficiency and Necessity-and-Sufficiency: Extending Condition Analysis in
#' Degree*.
#'
#' @section What it is for:
#' * Testing the claim that a condition is both required and sufficient for an
#'   outcome; the joint claim is supported only when both components are.
#' * Estimating the thresholds: for each outcome level, the condition level
#'   the outcome requires, the level that is sufficient for it, and the
#'   interval between them ([nsca_thresholds()]).
#' * Describing how closely the relationship approaches a
#'   necessary-and-sufficient one, through both component effect sizes and
#'   three summaries that combine them ([nsca_table()]).
#' * Comparing the result with regression ([nsca_compare()]) and planning the
#'   sample size of a study ([nsca_power()]).
#'
#' @section How to use it:
#' 1. The claim is fitted with [nsca_analysis()], which takes the condition,
#'    the outcome, the direction and, if theory provides one, the scope.
#'    `test.rep` sets the number of permutations for the two permutation
#'    tests; [nsca()] runs 1000 permutations and draws the plot in one step.
#' 2. The results are reported by [nsca_results()], [nsca_table()],
#'    [nsca_thresholds()] and [nsca_plot()]; [nsca_output()] prints and draws
#'    them together.
#' 3. [nsca_compare()] and [nsca_locationshift()] set the result beside
#'    regression.
#' 4. [nsca_test_plot()], [nsca_scope()], [nsca_scope_sensitivity()],
#'    [nsca_resample()] and [nsca_outliers()] provide diagnostic checks.
#'
#' @section Joint support:
#' The joint claim is supported only when both components are: both
#' directional permutation tests reject, combined by the intersection-union
#' rule (Berger, 1982), `p_nsca_iut = max(p_nec, p_suf)`; any pre-specified
#' relevance threshold is met; and the two estimated frontiers form an
#' acceptable geometry. The two tests can be based on one shared permutation
#' sequence, which also allows statistics involving both components to be
#' tested: the smaller of the two empty areas (`p_weakest_perm`) and the two
#' together (`p_sum_perm`). Each p-value addresses a different question, and
#' we suggest choosing the one that corresponds to the research question
#' before examining the results. Each component effect size can be read
#' against the conventional bands of Dul (2016). The joint summaries are
#' constructions specific to this package; they may be calibrated for the
#' design at hand with [nsca_power()] and [nsca_resample()].
#'
#' @section Reading the result:
#' A fitted model describes the empty-space pattern in the sample under the
#' stated scope and frontier technique, tested against random pairing of the
#' condition and the outcome, and its thresholds are empirical rules for the
#' population and design studied. Claims beyond the sample -- that the
#' condition causes the outcome, or that future cases will follow the
#' frontiers -- depend on the research design. Each direction predicts one
#' empty corner on each side, which assumes that the relationship is
#' monotone; it is therefore advisable to plot the data before choosing a
#' direction.
#'
#' @section What the permutation tests compare:
#' Each p-value in this package compares the observed empty areas with the
#' areas produced by randomly re-pairing the same values of the condition and
#' the outcome, the type of test that NCA also uses. A small p-value therefore
#' indicates that a corner is emptier than random pairing would leave it; on
#' its own it does not establish that the corner is empty in the population.
#'
#' This is most relevant when one corner is largely empty and the other is
#' not. A strong pattern on one side changes how the values of the outcome
#' are distributed; random pairing then fills the opposite corner more than
#' the data do, so that corner may appear emptier than under random pairing
#' even where it is not empty, and its test may reject more often than the
#' nominal level. A simple correlation can have the same effect on both
#' corners. In our simulations this tendency was noticeable, and its size
#' depended on the strength of the other side, the sample size and the
#' distributions involved; the scripts and their results are installed with
#' the package, in `system.file("simulations", package = "NSCA")`. The
#' intersection-union test does not remove it, since it combines the two
#' component tests as they are, and the two supplementary tests are more
#' susceptible to it (see [nsca_table()]).
#'
#' Researchers may wish to keep this in mind in three ways:
#' * reading each p-value together with the effect sizes, the scatter plot
#'   and theory, especially when one side is much stronger than the other;
#' * checking a planned design by simulating it with one side absent, for
#'   example `nsca_power(d_nec = 0.25, d_suf = 0)`, where every rejection of
#'   the joint claim is a false positive (see [nsca_power()]);
#' * stating in a report that the tests compare the data with random pairing.
#'
#' These are offered as considerations rather than rules; how much weight
#' they deserve depends on the design and the claim.
#'
#' @section Regression as a reference:
#' An ordinary least-squares or median regression line can be drawn beside the
#' two frontiers as a central-tendency reference. It summarizes the average
#' effect and is not a component of either claim. Chan (in press), *Mean
#' Effects Are Not Sufficiency Effects: How Sufficient Condition Analysis
#' Complements Regression Analysis*, explains why the two kinds of effect are
#' reported side by side.
#'
#' @seealso [nsca_analysis()], [nsca_table()], [nsca_thresholds()],
#'   [nsca_plot()], [nsca_terms()]
#' @references
#' Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance
#' sampling. *Technometrics*, 24(4), 295-300.
#'
#' Chan, Y. (in press). *From Necessity to Sufficiency and
#' Necessity-and-Sufficiency: Extending Condition Analysis in Degree*.
#' University of Leeds.
#'
#' Chan, Y. (in press). *Mean Effects Are Not Sufficiency Effects: How
#' Sufficient Condition Analysis Complements Regression Analysis*.
#' University of Leeds.
#'
#' Dul, J. (2016). Necessary Condition Analysis (NCA): Logic and methodology
#' of "necessary but not sufficient" causality. *Organizational Research
#' Methods*, 19(1), 10-52. \doi{10.1177/1094428115584005}
"_PACKAGE"
