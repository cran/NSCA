# Display helpers: one call that prints what a reported NSCA needs, and the
# one plot the package could not draw.

#' Plot the shared permutation null
#'
#' Draws the null distribution of each component effect size, of their
#' minimum and of their sum, from the single shuffle sequence
#' `shared.test.rep` produces. Each panel marks the observed value (solid
#' line) and the critical value (dashed line), that is, the `1 - alpha`
#' quantile that the observed value needs to exceed for the test to reject.
#'
#' @section Why this plot needs the shared sequence:
#' The engines' own permutation runs shuffle the outcome twice, independently,
#' so the pairs they produce are drawn from a product measure rather than from
#' the joint law. This does not affect `p_nec` and `p_suf`, whose marginals are
#' unchanged, but it invalidates any statistic that is a function of both. The
#' two empty spaces are disjoint subsets of one scope box, so `d_nec + d_suf`
#' is bounded and the two cannot both be large at the same time. Pairing a
#' large draw of one with a large draw of the other, taken from different
#' shuffles, yields combinations that no dataset could generate.
#'
#' `weakest_effect` and `component_sum` are functions of both components.
#' Their panels are therefore drawn only from
#' `attr(model, "nsca")$shared`, and this function
#' stops with an error when that sequence is absent rather than assembling a
#' substitute from the engines' separate draws. Such a substitute may appear
#' plausible, but it would place its upper tail too far to the right.
#'
#' @param model An object returned by [nsca_analysis()] with
#'   `shared.test.rep > 0`.
#' @param x A condition name or position. Defaults to the first condition.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param components Which panels to draw: any of `"necessity"`,
#'   `"sufficiency"`, `"weakest"` and `"sum"`.
#' @param bins Number of histogram bins.
#' @return A `ggplot` object.
#' @seealso [nsca_analysis()], [nsca_table()], [SCAtools::sca_test_plot()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' # A report would use 1000 or more replications.
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
#'                      shared.test.rep = 49, shared.seed = 1)
#' nsca_shared_plot(fit)
nsca_shared_plot <- function(model, x = NULL, ceiling = NULL,
                             components = c("necessity", "sufficiency",
                                            "weakest", "sum"),
                             bins = 30) {
  metadata <- .nsca_metadata(model)
  shared <- metadata$shared
  if (is.null(shared)) {
    stop(
      paste0(
        "No shared permutation sequence is stored in this model, so the null ",
        "distribution of the weakest effect is not available. It can be ",
        "obtained by refitting with shared.test.rep > 0. The engines' own ",
        "draws, produced by test.rep, come from two independent shuffles and ",
        "cannot be combined into a null distribution for a statistic of both ",
        "components."
      ),
      call. = FALSE
    )
  }

  components <- match.arg(components, several.ok = TRUE)
  # .nsca_select_conditions() returns every condition for NULL, which is right
  # for nsca_plot() and wrong here: one panel set shows one condition, so NULL
  # has to mean the first rather than all of them.
  if (is.null(x)) {
    x <- metadata$conditions[[1L]]
  }
  condition <- .nsca_select_conditions(model, x)
  if (length(condition) != 1L) {
    stop("A permutation-test plot requires exactly one condition.",
         call. = FALSE)
  }
  condition <- condition[[1L]]
  if (is.null(ceiling)) {
    ceiling <- metadata$ceilings[[1L]]
  }
  if (length(ceiling) != 1L || !(ceiling %in% metadata$ceilings)) {
    stop(
      sprintf(
        "Frontier '%s' was not estimated. Available: %s.",
        paste(ceiling, collapse = ", "),
        paste(metadata$ceilings, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  joint <- nsca_table(model)
  row <- joint[joint$condition == condition & joint$frontier == ceiling, ,
               drop = FALSE]
  if (nrow(row) == 0L) {
    stop(
      sprintf("No result for condition '%s' with frontier '%s'.",
              condition, ceiling),
      call. = FALSE
    )
  }
  row <- row[1L, ]

  labels <- c(
    necessity = "necessity  d_nec",
    sufficiency = "sufficiency  d_suf",
    weakest = "weakest  min(d_nec, d_suf)",
    sum = "sum  d_nec + d_suf"
  )
  observed <- c(
    necessity = row$d_nec,
    sufficiency = row$d_suf,
    weakest = row$weakest_effect,
    sum = row$component_sum
  )
  critical <- c(
    necessity = row$d_nec_critical,
    sufficiency = row$d_suf_critical,
    weakest = row$weakest_effect_critical,
    sum = row$component_sum_critical
  )
  p_values <- c(
    necessity = row$p_nec,
    sufficiency = row$p_suf,
    weakest = row$p_weakest_perm,
    sum = row$p_sum_perm
  )

  draws <- lapply(components, function(component) {
    values <- as.numeric(shared[[component]][, condition, ceiling])
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
      return(NULL)
    }
    data.frame(
      effect_size = values,
      component = factor(labels[[component]], levels = unname(labels)),
      stringsAsFactors = FALSE
    )
  })
  draws <- do.call(rbind, draws[!vapply(draws, is.null, logical(1L))])
  if (is.null(draws) || nrow(draws) == 0L) {
    stop("The shared sequence contains no usable draws for this selection.",
         call. = FALSE)
  }

  marks <- data.frame(
    component = factor(labels[components], levels = unname(labels)),
    observed = unname(observed[components]),
    critical = unname(critical[components]),
    stringsAsFactors = FALSE
  )
  bars <- marks[is.finite(marks$critical), , drop = FALSE]
  marks <- marks[is.finite(marks$observed), , drop = FALSE]

  replications <- dim(shared$weakest)[[1L]]
  requested <- metadata$shared.test.rep
  shortfall <- if (is.finite(requested) && replications < requested) {
    sprintf(" (%d of the %d requested succeeded)", replications, requested)
  } else {
    ""
  }

  chart <- ggplot2::ggplot(draws, ggplot2::aes(x = effect_size)) +
    ggplot2::geom_histogram(bins = bins, fill = "white", color = "gray25") +
    ggplot2::facet_wrap(~ component, ncol = 1L, scales = "free_y")
  if (nrow(bars) > 0L) {
    chart <- chart + ggplot2::geom_vline(
      data = bars,
      mapping = ggplot2::aes(xintercept = critical),
      color = "gray25", linetype = "dashed", linewidth = 0.6,
      inherit.aes = FALSE
    )
  }
  if (nrow(marks) > 0L) {
    chart <- chart + ggplot2::geom_vline(
      data = marks,
      mapping = ggplot2::aes(xintercept = observed),
      color = "#B22222", linewidth = 0.8, inherit.aes = FALSE
    )
  }

  caption <- if (isTRUE(metadata$test.rep > 0)) {
    paste0(
      "Solid line: observed value. Dashed line: critical value at alpha = ",
      format(metadata$alpha), ". A single shuffle sequence is used for both ",
      "components. The engines' separate runs, requested through test.rep, ",
      "are also stored in this model and may give slightly different ",
      "component p-values; nsca_table() reports the shared values, as ",
      "recorded in p_source."
    )
  } else {
    paste0(
      "Solid line: observed value. Dashed line: critical value at alpha = ",
      format(metadata$alpha), ". Each replication uses a single shuffle of ",
      "the outcome for both components."
    )
  }

  chart +
    ggplot2::labs(
      x = "Effect size under permutation",
      y = "Frequency",
      # Title, not caption: wraps narrower, because it is set larger and bold
      # and the panel is the same width for both.
      title = .nsca_wrap_title(
        paste0("NSCA shared permutation null: ",
               metadata$statement[[condition]])
      ),
      subtitle = .nsca_wrap(sprintf(
        "%s, %s frontier; %d replications%s. Observed: d_nec = %.3f (p = %.3f), d_suf = %.3f (p = %.3f), weakest = %.3f (p = %.3f), sum = %.3f (p = %.3f).",
        condition, ceiling, replications, shortfall,
        observed[["necessity"]], p_values[["necessity"]],
        observed[["sufficiency"]], p_values[["sufficiency"]],
        observed[["weakest"]], p_values[["weakest"]],
        observed[["sum"]], p_values[["sum"]]
      )),
      caption = .nsca_wrap(caption)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      strip.text = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

# ggplot2 wraps neither a title nor a caption, and a caption is right-aligned
# by default, so a long line loses its end on one and its beginning on the
# other with no warning. Wrapping here rather than asking the caller to guess
# an export width.
.nsca_wrap <- function(text, width = 95L) {
  if (!is.character(text) || length(text) != 1L || is.na(text)) {
    return(text)
  }
  paste(strwrap(text, width = width), collapse = "\n")
}

#' Display selected NSCA output
#'
#' Prints the tables and draws the plots that are typically needed to report a
#' necessary-and-sufficient analysis, in a single call. This is the
#' counterpart of [SCAtools::sca_output()].
#'
#' @param model An object returned by [nsca_analysis()].
#' @param summaries Whether to print the three result layers.
#' @param thresholds Whether to print the dual-threshold table.
#' @param plots Whether to draw the frontier plots.
#' @param tests Whether to draw the permutation-test plots. Both components are
#'   drawn in one figure by [nsca_test_plot()], [NCA::nca_output()] adds the
#'   engine's text on the p-accuracy, and the shared sequence, when one was
#'   run, is drawn by [nsca_shared_plot()].
#' @param selection Optional condition names or positions.
#' @param ceiling Frontier technique for the threshold and test output.
#'   Defaults to the first technique estimated.
#' @param scale,outcome_scale,convention Optional reporting overrides passed to
#'   [nsca_thresholds()].
#' @param inequality `"strict"` or `"inclusive"`.
#' @return Invisibly, a list of the tables and plots produced.
#' @seealso [nsca_results()], [nsca_plot()], [nsca_test_plot()],
#'   [nsca_shared_plot()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh")
#' nsca_output(fit, thresholds = TRUE)
nsca_output <- function(
    model,
    summaries = TRUE,
    thresholds = FALSE,
    plots = TRUE,
    tests = FALSE,
    selection = NULL,
    ceiling = NULL,
    scale = NULL,
    outcome_scale = NULL,
    convention = NULL,
    inequality = c("strict", "inclusive")) {
  metadata <- .nsca_metadata(model)
  inequality <- match.arg(inequality)
  conditions <- .nsca_select_conditions(model, selection)
  if (length(conditions) == 0L) {
    stop("No valid conditions were selected.", call. = FALSE)
  }
  if (is.null(ceiling)) {
    ceiling <- metadata$ceilings[[1L]]
  }

  result <- list()

  if (summaries) {
    layers <- nsca_results(model)
    keep <- function(table) {
      table[table$condition %in% conditions, , drop = FALSE]
    }
    result$necessity <- keep(layers$necessity)
    result$sufficiency <- keep(layers$sufficiency)
    result$necessary_and_sufficient <- keep(layers$necessary_and_sufficient)

    cat("\n1. Necessary Condition Analysis\n")
    print(result$necessity[c("condition", "frontier", "effect_size",
                             "p_value", "frontier_accuracy", "fit")],
          row.names = FALSE)
    cat("\n2. Sufficient Condition Analysis\n")
    print(result$sufficiency[c("condition", "frontier", "effect_size",
                               "p_value", "frontier_accuracy", "fit")],
          row.names = FALSE)
    cat("\n3. Necessary and Sufficient Condition Analysis\n")
    joint_columns <- c(
      "condition", "frontier", "d_nec", "d_suf", "weakest_effect",
      "balanced_joint_effect", "admissible_region_share", "p_nsca_iut"
    )
    if (!is.null(metadata$shared)) {
      joint_columns <- c(joint_columns, "p_weakest_perm", "p_sum_perm")
    }
    print(result$necessary_and_sufficient[c(joint_columns, "joint_support")],
          row.names = FALSE)
  }

  if (thresholds) {
    result$thresholds <- nsca_thresholds(
      model, ceiling = ceiling, x = conditions, scale = scale,
      outcome_scale = outcome_scale, convention = convention,
      inequality = inequality
    )
    cat("\n4. Dual thresholds\n")
    print(
      result$thresholds[c(
        "condition", "outcome_level", "necessity_threshold",
        "sufficiency_threshold", "necessity_sufficiency_interval",
        "threshold_status"
      )],
      row.names = FALSE
    )
  }

  if (plots) {
    result$plots <- nsca_plot(model, x = conditions)
    if (inherits(result$plots, "ggplot")) {
      print(result$plots)
    } else {
      invisible(lapply(result$plots, print))
    }
  }

  if (tests) {
    # Both components in one figure. Until 0.5.0 this drew the sufficiency
    # null through SCAtools and left the necessity null undrawn, printing the
    # engine's text for it instead -- so the two halves of a two-sided
    # analysis were reported in two different media, and the reader could not
    # see the one comparison the decision rests on.
    result$component_tests <- lapply(conditions, function(condition) {
      tryCatch(
        nsca_test_plot(model, x = condition, ceiling = ceiling),
        error = function(e) NULL
      )
    })
    names(result$component_tests) <- conditions
    invisible(lapply(Filter(Negate(is.null), result$component_tests), print))

    # The engine's own text as well, for the p-accuracy the plot has no room
    # for: a p-value of 0.001 from 1000 replications has a margin of error
    # that decides whether it is separable from 0.01.
    if (isTRUE(metadata$test.rep > 0)) {
      NCA::nca_output(model$necessity, plots = FALSE, summaries = FALSE,
                      test = TRUE)
    }

    # The one distribution neither engine can draw.
    if (!is.null(metadata$shared)) {
      result$shared_plots <- lapply(conditions, function(condition) {
        tryCatch(
          nsca_shared_plot(model, x = condition, ceiling = ceiling),
          error = function(e) NULL
        )
      })
      names(result$shared_plots) <- conditions
      invisible(lapply(Filter(Negate(is.null), result$shared_plots), print))
    }
  }

  # This function assembles a report, and a report has to say what it was run
  # with -- here including the two geometry tolerances that decide
  # joint_support. print() ends with the same block.
  .nsca_settings_block(metadata)
  invisible(result)
}
