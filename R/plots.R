# Plotting.
#
# Several frontier techniques may be drawn on one chart, because comparing them
# by eye is exactly what a chart is for. The restriction to a single technique
# applies to the computation only: there, the two sides are compared
# numerically, and a comparison between differently measured quantities is not
# interpretable. Here nothing is compared, so every requested technique is
# drawn on both sides and labeled.

.nsca_select_conditions <- function(model, x) {
  available <- .nsca_metadata(model)$conditions
  if (is.null(x)) {
    return(available)
  }
  if (is.numeric(x)) {
    x <- available[x]
  }
  unknown <- setdiff(x, available)
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown condition(s): %s.", paste(unknown, collapse = ", ")),
         call. = FALSE)
  }
  # NULL already means every condition, so an empty selection can only be a
  # filter that matched nothing. Saying so beats carrying it downstream.
  if (length(x) == 0L) {
    stop(
      sprintf(
        "'x' selected no conditions. NULL selects all of them; the available conditions are: %s.",
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  x
}

.nsca_frontier_frame <- function(metadata, model, condition, ceilings) {
  raw <- metadata$raw$conditions[[condition]]
  grid <- seq(raw$x_bounds[[1L]], raw$x_bounds[[2L]],
              length.out = .nsca_grid_size)
  nec_corner <- unname(metadata$necessity_corner[[condition]])
  suf_corner <- unname(metadata$sufficiency_corner[[condition]])

  frames <- list()
  index <- 0L
  for (ceiling in ceilings) {
    sides <- list(
      necessity = .nsca_frontier(
        raw$x, raw$y, nec_corner, ceiling, grid,
        slope = .nsca_nca_param(model$necessity, condition, ceiling, "Slope"),
        intercept = .nsca_nca_param(model$necessity, condition, ceiling, "Intercept"),
        peers = model$necessity$peers[[ceiling]][[condition]]
      ),
      sufficiency = .nsca_frontier(
        raw$x, raw$y, suf_corner, ceiling, grid,
        slope = .nsca_sca_param(model$sufficiency, condition, ceiling, "Slope"),
        intercept = .nsca_sca_param(model$sufficiency, condition, ceiling, "Intercept"),
        peers = model$sufficiency$peers[[ceiling]][[condition]]
      )
    )
    for (side in names(sides)) {
      # Clipped to the scope, because the area is. `.nsca_band()` integrates
      # the frontier against y_bounds, so a regression frontier that leaves
      # the box contributes nothing above the top edge -- while the drawn
      # line carried on and crossed the dashed scope boundary, which reads as
      # an error to anyone checking the figure against the number. Clipping
      # the line makes it the boundary of the region actually measured.
      values <- pmin(pmax(sides[[side]], raw$y_bounds[[1L]]),
                     raw$y_bounds[[2L]])
      if (!any(is.finite(values))) {
        next
      }
      index <- index + 1L
      frames[[index]] <- data.frame(
        x = grid, y = values, side = side, ceiling = ceiling,
        series = paste(side, ceiling, sep = " / "),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(frames) == 0L) {
    return(data.frame(x = numeric(), y = numeric(), side = character(),
                      ceiling = character(), series = character()))
  }
  out <- do.call(rbind, frames)
  out[is.finite(out$y), , drop = FALSE]
}

.nsca_plot_one <- function(model, condition, ceilings, shade, annotate) {
  metadata <- .nsca_metadata(model)
  raw <- metadata$raw$conditions[[condition]]
  points <- data.frame(x = raw$x, y = raw$y)
  lines <- .nsca_frontier_frame(metadata, model, condition, ceilings)
  xb <- raw$x_bounds
  yb <- raw$y_bounds

  chart <- ggplot2::ggplot(points, ggplot2::aes(x = x, y = y))

  if (shade && nrow(lines) > 0L) {
    # Shade the admissible region under the first technique: the region
    # neither component rules out. Its area is admissible_region_share.
    first <- lines[lines$ceiling == ceilings[[1L]], , drop = FALSE]
    nec <- first[first$side == "necessity", c("x", "y"), drop = FALSE]
    suf <- first[first$side == "sufficiency", c("x", "y"), drop = FALSE]
    if (nrow(nec) > 0L && nrow(suf) > 0L) {
      upper_is_nec <- identical(
        .nsca_corner_side(unname(metadata$necessity_corner[[condition]])), "upper"
      )
      upper <- if (upper_is_nec) nec else suf
      lower <- if (upper_is_nec) suf else nec
      # Join on the grid rather than assuming the two sides kept the same rows.
      # A frontier that is undefined somewhere loses those rows, and pairing by
      # position would then shift one frontier against the other silently.
      paired <- merge(upper, lower, by = "x", suffixes = c(".up", ".lo"))
      paired <- paired[order(paired$x), , drop = FALSE]
      ribbon <- data.frame(
        x = paired$x,
        ymax = pmin(pmax(paired$y.up, yb[[1L]]), yb[[2L]]),
        ymin = pmin(pmax(paired$y.lo, yb[[1L]]), yb[[2L]])
      )
      ribbon$ymin <- pmin(ribbon$ymin, ribbon$ymax)
      if (nrow(ribbon) > 1L) {
        chart <- chart + ggplot2::geom_ribbon(
          data = ribbon,
          mapping = ggplot2::aes(x = x, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE, fill = "#EAF2FB", color = NA
        )
      }
    }
  }

  chart <- chart +
    ggplot2::geom_point(shape = 21, color = "#174A7E", fill = "white",
                        size = 2) +
    ggplot2::geom_vline(xintercept = xb, linetype = "dashed",
                        color = "gray65") +
    ggplot2::geom_hline(yintercept = yb, linetype = "dashed",
                        color = "gray65")

  if (nrow(lines) > 0L) {
    # The linetype key is mapped only when more than one technique is drawn.
    # With a single technique it is a legend with one entry that distinguishes
    # nothing, and ggplot draws its sample line in black -- beside a color
    # legend in which every color means one of the two components, a black
    # key reads as a third component that does not exist.
    mapping <- if (length(ceilings) > 1L) {
      ggplot2::aes(x = x, y = y, color = side, linetype = ceiling)
    } else {
      ggplot2::aes(x = x, y = y, color = side)
    }
    chart <- chart + ggplot2::geom_line(
      data = lines, mapping = mapping, inherit.aes = FALSE, linewidth = 0.8
    )
  }

  # The central-tendency line, when one was requested. It stays out of the
  # component legend on purpose -- coloring it as a component would put an
  # average effect and an empty-space frontier in one key -- but it is
  # labeled *on the line*, because until 0.5.0 it had no key at all and a
  # reader saw a gray line through the middle of the cloud with nothing but a
  # caption sentence to identify it.
  if (length(metadata$reference) > 0L) {
    reference <- nsca_reference(model, x = condition,
                                line = metadata$reference)
    for (i in seq_len(nrow(reference))) {
      if (!is.finite(reference$slope[[i]]) ||
          !is.finite(reference$intercept[[i]])) {
        next
      }
      chart <- chart + ggplot2::geom_abline(
        intercept = reference$intercept[[i]],
        slope = reference$slope[[i]],
        color = "gray40",
        linetype = if (identical(reference$line[[i]], "ols")) "dotted"
                   else "longdash",
        linewidth = 0.7
      )
      chart <- chart + .nsca_reference_label(
        reference[i, ], xb, yb, metadata$reference.tau, i
      )
    }
  }

  if (annotate) {
    nec_corner <- unname(metadata$necessity_corner[[condition]])
    suf_corner <- unname(metadata$sufficiency_corner[[condition]])
    dx <- diff(xb)
    dy <- diff(yb)
    place <- function(corner, label) {
      ggplot2::annotate(
        "label",
        x = if (corner %in% c(1L, 3L)) xb[[1L]] + 0.03 * dx else xb[[2L]] - 0.03 * dx,
        y = if (corner %in% c(1L, 2L)) yb[[2L]] - 0.03 * dy else yb[[1L]] + 0.03 * dy,
        label = label,
        hjust = if (corner %in% c(1L, 3L)) 0 else 1,
        vjust = if (corner %in% c(1L, 2L)) 1 else 0,
        size = 3, fill = "white",
        # Each label sits inside a corner the claim says is empty, which is
        # exactly where a counterexample would show. An opaque box there hides
        # the observation a reader most needs to see. annotate = FALSE removes
        # both labels.
        alpha = 0.85
      )
    }
    chart <- chart +
      place(nec_corner, sprintf("empty if necessary\ncorner %d", nec_corner)) +
      place(suf_corner, sprintf("empty if sufficient\ncorner %d", suf_corner))
  }

  chart +
    ggplot2::coord_cartesian(xlim = xb, ylim = yb) +
    ggplot2::labs(
      x = condition,
      y = metadata$outcome,
      title = .nsca_wrap_title(paste0("NSCA: ", metadata$statement[[condition]])),
      subtitle = .nsca_wrap(sprintf(
        "Direction %s; necessity corner %d, sufficiency corner %d. The shaded band is the region that neither claim excludes.",
        unname(metadata$direction[[condition]]),
        unname(metadata$necessity_corner[[condition]]),
        unname(metadata$sufficiency_corner[[condition]])
      )),
      color = "Component",
      linetype = "Frontier",
      caption = .nsca_wrap(paste0(
        .nsca_reference_caption(metadata),
        "Empty-space patterns are empirical; a causal ",
        "necessary-and-sufficient claim requires a research design that ",
        "supports it."
      ))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "gray35", fill = NA),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9),
      # Left, not right. A wrapped caption aligned right has a ragged left
      # edge, and the caption is the last sentence of the figure.
      plot.caption = ggplot2::element_text(hjust = 0, size = 8)
    )
}

#' Plot both frontiers
#'
#' Draws the scatter plot with the necessity and sufficiency frontiers and
#' shades the admissible region between them. Several frontier techniques may
#' be drawn at once. Unlike the computation, a chart involves no numerical
#' comparison, so combining techniques here is appropriate and may be
#' informative.
#'
#' @param model An object returned by [nsca_analysis()].
#' @param x Optional condition names or positions.
#' @param ceilings Frontier techniques to draw. Defaults to all that were
#'   estimated.
#' @param shade Whether to shade the admissible region under the first
#'   technique.
#' @param annotate Whether to label the two corners that the claim requires to
#'   be empty. `FALSE` removes both labels, which may be preferable for small
#'   multiples or for a figure with its own caption.
#' @param ... Passed on to `nsca_plot()` by the `plot()` method.
#' @return A `ggplot` for one condition, or a named list of plots.
#' @seealso [nsca_analysis()], [nsca_table()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = c("ce_fdh", "cr_fdh"),
#'                      reference = "ols")
#' nsca_plot(fit)
#' plot(fit)
nsca_plot <- function(model, x = NULL, ceilings = NULL, shade = TRUE,
                      annotate = TRUE) {
  metadata <- .nsca_metadata(model)
  conditions <- .nsca_select_conditions(model, x)
  if (is.null(ceilings)) {
    ceilings <- metadata$ceilings
  }
  unknown <- setdiff(ceilings, metadata$ceilings)
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Frontier(s) not estimated: %s. Available: %s.",
        paste(unknown, collapse = ", "), paste(metadata$ceilings, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  plots <- lapply(
    conditions,
    function(condition) .nsca_plot_one(model, condition, ceilings, shade, annotate)
  )
  names(plots) <- conditions
  if (length(plots) == 1L) plots[[1L]] else plots
}

#' @rdname nsca_plot
#' @export
plot.nsca_result <- function(x, ...) {
  plots <- nsca_plot(x, ...)
  if (inherits(plots, "ggplot")) {
    print(plots)
  } else {
    invisible(lapply(plots, print))
  }
  invisible(plots)
}

# One component's permutation draws, from whichever engine produced them.
#
# The two engines store their tests in the same shape -- NCA writes
# `list(data, observed, p_value, threshold.value, names)` per condition per
# ceiling, and SCAtools passes `test.rep` through to the same engine -- so one
# reader serves both sides. Returns NULL rather than raising when a side has
# no test, because the caller has to distinguish "no test was run at all",
# which is a mistake worth a message, from "this frontier was not tested",
# which is not.
.nsca_component_test <- function(engine, condition, ceiling) {
  tests <- engine$tests
  if (is.null(tests) || is.null(tests[[condition]])) {
    return(NULL)
  }
  test <- tests[[condition]][[ceiling]]
  if (is.null(test) || is.null(test$data)) {
    return(NULL)
  }
  values <- as.numeric(test$data)
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(NULL)
  }
  list(
    values = values,
    observed = suppressWarnings(as.numeric(test$observed)),
    p_value = suppressWarnings(as.numeric(test$p_value))
  )
}

#' Plot the two component permutation tests
#'
#' Draws the null distribution of each component effect size against the
#' observed value: necessity in the upper panel, sufficiency in the lower. This
#' is the counterpart of [SCAtools::sca_test_plot()], which draws one panel
#' because a sufficiency analysis has one component.
#'
#' Each panel marks the observed effect size (solid line) and the critical
#' effect size (dashed line), the `1 - alpha` quantile of that component's
#' null distribution.
#'
#' These are the engines' own tests, requested through `test.rep`. Each
#' provides a valid null distribution for its own component, and the
#' intersection-union p-value `max(p_nec, p_suf)` reported in the subtitle is
#' valid whether or not the two were drawn together, because it rejects only
#' when both components reject and each component's test controls its own
#' error rate. These two panels cannot, however, provide a null distribution
#' for a statistic that depends on both components at once, such as
#' `min(d_nec, d_suf)`: the two histograms come from two independent shuffles,
#' so the pairs behind them were not drawn from one joint law.
#' [nsca_shared_plot()] serves that purpose and requires
#' `shared.test.rep > 0`.
#'
#' @param model An object returned by [nsca_analysis()] with `test.rep > 0`.
#' @param x A condition name or position. Defaults to the first condition.
#' @param ceiling Frontier technique. Defaults to the first technique
#'   estimated.
#' @param bins Number of histogram bins.
#' @return A `ggplot` object.
#' @seealso [nsca_shared_plot()], [nsca_plot()], [SCAtools::sca_test_plot()]
#' @export
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' dat <- data.frame(X = x, Y = pmin(pmax(x + rnorm(60, 0, 0.1), 0), 1))
#' # A report would use 1000 or more permutations.
#' fit <- nsca_analysis(dat, "X", "Y", ceilings = "ce_fdh",
#'                      test.rep = 99, test.seed = 1)
#' nsca_test_plot(fit)
nsca_test_plot <- function(model, x = NULL, ceiling = NULL, bins = 30) {
  metadata <- .nsca_metadata(model)

  # NULL means every condition to nsca_plot() and the first one here: one
  # panel set shows one condition.
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

  necessity <- .nsca_component_test(model$necessity, condition, ceiling)
  sufficiency <- .nsca_component_test(model$sufficiency, condition, ceiling)
  if (is.null(necessity) && is.null(sufficiency)) {
    stop(
      sprintf(
        paste0(
          "No permutation draws are stored for '%s' with frontier '%s'. ",
          "This model was fitted with test.rep = %s; a refit with ",
          "test.rep > 0 is required."
        ),
        condition, ceiling, format(metadata$test.rep)
      ),
      call. = FALSE
    )
  }

  labels <- c(necessity = "necessity  d_nec",
              sufficiency = "sufficiency  d_suf")
  parts <- list(necessity = necessity, sufficiency = sufficiency)
  present <- names(parts)[!vapply(parts, is.null, logical(1L))]

  draws <- do.call(rbind, lapply(present, function(side) {
    data.frame(
      effect_size = parts[[side]]$values,
      component = factor(labels[[side]], levels = unname(labels)),
      stringsAsFactors = FALSE
    )
  }))

  marks <- data.frame(
    component = factor(labels[present], levels = unname(labels)),
    observed = vapply(present, function(side) parts[[side]]$observed,
                      numeric(1L)),
    critical = vapply(present, function(side) {
      .nsca_critical(parts[[side]]$values, metadata$alpha)
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )
  bars <- marks[is.finite(marks$critical), , drop = FALSE]
  marks <- marks[is.finite(marks$observed), , drop = FALSE]

  p_nec <- if (is.null(necessity)) NA_real_ else necessity$p_value
  p_suf <- if (is.null(sufficiency)) NA_real_ else sufficiency$p_value
  p_iut <- suppressWarnings(max(p_nec, p_suf))

  # One panel when only one side was tested: a facet strip that names the only
  # component present, over an empty second panel, says less than the title.
  chart <- ggplot2::ggplot(draws, ggplot2::aes(x = effect_size)) +
    ggplot2::geom_histogram(bins = bins, fill = "white", color = "gray25") +
    ggplot2::facet_wrap(~ component, ncol = 1L, scales = "free_y",
                        drop = TRUE)
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

  subtitle <- sprintf(
    "%s, %s frontier. d_nec = %s (p = %s), d_suf = %s (p = %s); p_nsca_iut = max(p_nec, p_suf) = %s.",
    condition, ceiling,
    .nsca_format_number(if (is.null(necessity)) NA_real_ else necessity$observed),
    .nsca_format_number(p_nec),
    .nsca_format_number(if (is.null(sufficiency)) NA_real_ else sufficiency$observed),
    .nsca_format_number(p_suf),
    .nsca_format_number(p_iut)
  )

  chart +
    ggplot2::labs(
      x = "Effect size under permutation",
      y = "Frequency",
      title = .nsca_wrap_title(
        paste0("NSCA component permutation tests: ",
               metadata$statement[[condition]])
      ),
      subtitle = .nsca_wrap(subtitle),
      caption = .nsca_wrap(paste0(
        "Solid line: observed effect size. Dashed line: critical effect size ",
        "at alpha = ", format(metadata$alpha), ". The panels come from two ",
        "independent shuffles, one per component: each panel is the null ",
        "distribution for its own component, and together they give the ",
        "intersection-union p-value. nsca_shared_plot() draws the null ",
        "distribution for statistics of both components."
      ))
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

# %.3f prints "NA" as "NA" only for NA_real_ passed to sprintf in some
# locales and errors on NULL, and an untested component legitimately has no
# number at all. One formatter, so a missing side reads as "not tested"
# rather than as a number that failed to print.
.nsca_format_number <- function(value) {
  if (length(value) != 1L || !is.finite(value)) {
    return("not tested")
  }
  sprintf("%.3f", value)
}

# Name the unkeyed reference lines, since they are kept out of the component
# legend on purpose: they belong to neither side, and giving them a color in
# that key would put an average effect and an empty-space frontier in the same
# legend. A median line and a 95th-percentile line look identical on paper, so
# the quantile has to be stated.
.nsca_reference_caption <- function(metadata) {
  reference <- metadata$reference
  if (length(reference) == 0L) {
    return("")
  }
  tau <- metadata$reference.tau
  if (is.null(tau) || !is.finite(tau)) {
    tau <- 0.5
  }
  described <- vapply(reference, function(line) {
    if (identical(line, "ols")) {
      "dotted gray: OLS central tendency"
    } else {
      sprintf("long-dashed gray: quantile regression at tau = %g", tau)
    }
  }, character(1L))
  paste0(paste(described, collapse = "; "), ", drawn for comparison only. ")
}

# The rates nsca_powerplot() can draw, by the name `show` uses.
.nsca_power_labels <- c(
  necessity   = "necessity",
  sufficiency = "sufficiency",
  joint       = "joint (both reject)",
  weakest     = "weakest (p_weakest_perm)",
  sum         = "sum (p_sum_perm)"
)
.nsca_power_columns <- c(
  necessity   = "power_necessity",
  sufficiency = "power_sufficiency",
  joint       = "power_joint",
  weakest     = "power_weakest_perm",
  sum         = "power_sum_perm"
)

# The long form nsca_powerplot() draws from.
#
# `nsca_power()` returns one row per sample size with three power columns, not
# one column with a grouping variable, because the three are different
# quantities rather than three levels of one factor. A plot needs them stacked,
# and the stacking is what decides the legend order, so it lives here rather
# than in the plotting code.
.nsca_power_long <- function(df_power, show) {
  labels <- .nsca_power_labels
  columns <- .nsca_power_columns
  rows <- lapply(show, function(component) {
    data.frame(
      n = df_power$n,
      power = df_power[[columns[[component]]]],
      component = labels[[component]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # Joint FIRST, so it is drawn underneath and the two component curves show
  # through where all three coincide -- which is the ordinary case at high
  # power. Until 0.5.0 the joint line was drawn last and painted over them, so
  # a figure with three legend entries showed one curve.
  #
  # It leads the legend for the same reason it is drawn first: it is the
  # curve a study is powered on, and the caption says so.
  ordered <- c(labels[["joint"]],
               unname(labels[setdiff(show, "joint")]))
  out$component <- factor(out$component,
                          levels = intersect(ordered, unname(labels[show])))
  out[is.finite(out$power), , drop = FALSE]
}

#' Plot a necessary-and-sufficient power analysis
#'
#' Draws the three rejection rates that [nsca_power()] reports against sample
#' size. This is the counterpart of [SCAtools::sca_powerplot()], which draws
#' one curve because a sufficiency analysis has one test.
#'
#' Of the three curves, the joint curve, which is the lowest, is the relevant
#' one for planning. `joint` is the rate at which the intersection-union test
#' rejects, which requires both components to reject, so a study is adequately
#' powered when this curve reaches the target, not when either component curve
#' does. It is drawn with the heaviest line for that reason.
#'
#' `power_if_independent` is drawn as a dashed reference when `reference` is
#' `TRUE`. It is the joint power the two tests would have if they rejected
#' independently, which they do not, so the vertical distance between it and
#' the joint curve shows their dependence. When the joint curve lies above it,
#' the components tend to reject together and a study sized on the product
#' would be larger than necessary; when it lies below, they tend to trade off
#' and the product is optimistic.
#'
#' The curves are subject to the caveats of [nsca_power()]: the simulated
#' design is a simplification, and the realized effect sizes tend to overshoot
#' the population values at small n. The curves are therefore best read
#' alongside `d_nec_realized`, `d_suf_realized` and the declared effect sizes
#' rather than in isolation.
#'
#' `"weakest"` and `"sum"` draw the power of the two supplementary tests,
#' which [nsca_power()] reports with `shared = TRUE`. They address questions
#' other than the joint claim and are drawn only when requested.
#'
#' @param df_power A data frame from [nsca_power()].
#' @param show Which rates to draw: any of `"necessity"`, `"sufficiency"`,
#'   `"joint"`, `"weakest"` and `"sum"`.
#' @param reference Whether to draw `power_if_independent` as a dashed
#'   reference line.
#' @param target Optional horizontal line at a conventional target, `0.80` by
#'   default. `NULL` draws none.
#' @return A `ggplot` object.
#' @seealso [nsca_power()], [SCAtools::sca_powerplot()]
#' @export
#' @examples
#' power <- nsca_power(n = c(30, 60), rep = 2, test.rep = 39, seed = 1)
#' nsca_powerplot(power)
nsca_powerplot <- function(df_power,
                           show = c("necessity", "sufficiency", "joint"),
                           reference = TRUE,
                           target = 0.80) {
  show <- match.arg(show, names(.nsca_power_labels), several.ok = TRUE)
  required <- c("n", "power_necessity", "power_sufficiency", "power_joint")
  missing <- setdiff(required, colnames(df_power))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        paste0(
          "Missing power-result column(s): %s. The function expects the ",
          "full data frame returned by nsca_power(), not a subset of it."
        ),
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (nrow(df_power) == 0L) {
    stop("The power table has no rows to draw.", call. = FALSE)
  }
  supplementary <- intersect(show, c("weakest", "sum"))
  absent <- supplementary[vapply(supplementary, function(rate) {
    column <- .nsca_power_columns[[rate]]
    !column %in% colnames(df_power) || !any(is.finite(df_power[[column]]))
  }, logical(1L))]
  if (length(absent) > 0L) {
    stop(
      sprintf(
        paste0(
          "No %s power to draw: %s is reported only by ",
          "nsca_power(shared = TRUE)."
        ),
        paste(absent, collapse = " or "),
        paste(.nsca_power_columns[absent], collapse = " and ")
      ),
      call. = FALSE
    )
  }

  long <- .nsca_power_long(df_power, show)
  if (nrow(long) == 0L) {
    stop(
      paste0(
        "No finite power values to draw. Every replication failed at every ",
        "sample size, as recorded by nsca_power() in `shortfall`."
      ),
      call. = FALSE
    )
  }

  chart <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = n, y = power, color = component, group = component)
  )

  if (isTRUE(reference) && "power_if_independent" %in% colnames(df_power)) {
    benchmark <- df_power[is.finite(df_power$power_if_independent), , drop = FALSE]
    if (nrow(benchmark) > 1L) {
      chart <- chart + ggplot2::geom_line(
        data = benchmark,
        mapping = ggplot2::aes(x = n, y = power_if_independent),
        inherit.aes = FALSE,
        color = "gray45", linetype = "dashed", linewidth = 0.6
      )
    }
  }

  if (!is.null(target) && is.finite(target)) {
    chart <- chart + ggplot2::geom_hline(
      yintercept = target, color = "gray70", linetype = "dotted",
      linewidth = 0.6
    )
  }

  # The joint curve heaviest: it is the one a study is powered on.
  chart +
    ggplot2::geom_line(
      ggplot2::aes(linewidth = component == "joint (both reject)",
                   linetype = component)
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = component), size = 2.1) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.7, `TRUE` = 1.3),
                                    guide = "none") +
    # Three curves that coincide are the ordinary case at high power, and
    # a solid line painted over a solid line leaves a legend naming three
    # series of which one can be seen.
    ggplot2::scale_linetype_manual(
      values = c(`necessity` = "longdash", `sufficiency` = "dotted",
                 `joint (both reject)` = "solid",
                 `weakest (p_weakest_perm)` = "dotdash",
                 `sum (p_sum_perm)` = "twodash"),
      name = "Test"
    ) +
    ggplot2::scale_shape_manual(
      values = c(`necessity` = 16, `sufficiency` = 17,
                 `joint (both reject)` = 15,
                 `weakest (p_weakest_perm)` = 18,
                 `sum (p_sum_perm)` = 4),
      name = "Test"
    ) +
    # Break at the sample sizes that were actually simulated. A power
    # curve is read by picking an n off the axis, and the default breaks
    # landed on sizes the study never ran.
    ggplot2::scale_x_continuous(breaks = sort(unique(long$n))) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Sample size",
      y = "Rejection rate",
      color = "Test",
      title = "Power for a necessary-and-sufficient claim",
      caption = .nsca_wrap(paste0(
        "The joint curve requires both components to reject, so it is the ",
        "relevant curve for planning a study. The dashed line is the joint ",
        "power the two tests would have if they rejected independently, which ",
        "they do not; the gap between it and the joint curve reflects their ",
        "dependence. The simulated design is a simplification, and the ",
        "realized effect sizes tend to overshoot the population values at ",
        "small n, so these curves are best read alongside d_nec_realized and ",
        "d_suf_realized.",
        if (length(supplementary) > 0L) {
          paste0(
            " The weakest and sum curves show the supplementary tests on the ",
            "same shuffles; they address questions other than the joint claim."
          )
        } else {
          ""
        }
      ))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

# Titles wrap narrower than captions, because they are set in a larger,
# bolder face and the panel is the same width for both. 95 characters of
# caption fit; 95 of title overflowed and were clipped mid-word. The number
# below was found by rendering, not by arithmetic on point sizes.
.nsca_wrap_title <- function(text) {
  .nsca_wrap(text, width = 62L)
}

# Name the reference line where it is drawn.
#
# ggplot2 allows one scale per aesthetic, and both are taken: color is the
# component and linetype is the technique. A third key would mean bending one
# of those, which is the conflation the caller comment objects to. A label on
# the line names it without claiming it is either.
#
# Placed three quarters of the way along x, above the line, which keeps it
# clear of the two corner boxes at top-left and bottom-right. Clamped into the
# scope so a steep line cannot push it off the panel.
#
# Staggered by index, because an OLS line and a median line very nearly
# coincide on symmetric data -- and two labels at the same place is one
# label with something hidden under it.
.nsca_reference_label <- function(row, xb, yb, tau, index = 1L) {
  at <- xb[[1L]] + max(0.74 - 0.19 * (index - 1L), 0.15) * diff(xb)
  y <- row$intercept[[1L]] + row$slope[[1L]] * at
  lift <- 0.035 * diff(yb)
  y <- min(max(y + lift, yb[[1L]] + lift), yb[[2L]] - lift)
  label <- if (identical(row$line[[1L]], "ols")) {
    "OLS"
  } else {
    sprintf("QR tau=%g", if (is.null(tau) || !is.finite(tau)) 0.5 else tau)
  }
  ggplot2::annotate(
    "label", x = at, y = y, label = label,
    size = 2.7, color = "gray30", fill = "white",
    label.size = 0, label.padding = ggplot2::unit(0.12, "lines"),
    hjust = 0.5, vjust = 0.5
  )
}
