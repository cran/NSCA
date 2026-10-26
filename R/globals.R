# Variables captured by ggplot2's data-masked aesthetics.
utils::globalVariables(c(
  "side", "x", "y", "ymax", "ymin",
  # nsca_shared_plot() and nsca_test_plot()
  "effect_size", "component", "observed", "critical",
  # nsca_powerplot()
  "n", "power", "power_if_independent"
))
