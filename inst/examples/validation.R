library(NSCA)

# A visual check of all four directions. X2 is intentionally left unsorted,
# because sorting both X2 and an outcome already ordered by X1 would create a
# strong monotone relationship, which would not be a valid check of
# independence.
set.seed(11)
n <- 70L
x <- sort(runif(n))
base_y <- pmin(pmax(x + rnorm(n, 0, 0.11), 0), 1)

# In all four logical directions, the two components use diagonally opposite
# physical corners.
fits <- lapply(
  c("HH", "LH", "HL", "LL"),
  function(direction) {
    y <- if (direction %in% c("HH", "LL")) base_y else 1 - base_y
    nsca_analysis(
      data.frame(X = x, Y = y), "X", "Y",
      direction = direction, ceilings = "ce_fdh"
    )
  }
)
names(fits) <- c("HH", "LH", "HL", "LL")
stopifnot(vapply(fits, function(fit) inherits(nsca_plot(fit), "ggplot"), TRUE))

# X2 is unrelated to Y and is therefore expected to leave a much wider band.
two <- nsca_analysis(
  data.frame(X1 = x, X2 = runif(n), Y = base_y),
  c("X1", "X2"), "Y", ceilings = "ce_fdh"
)
joint <- nsca_table(two)
band_x1 <- joint$admissible_region_share[joint$condition == "X1"]
band_x2 <- joint$admissible_region_share[joint$condition == "X2"]
stopifnot(band_x2 > 0.4, band_x2 > band_x1 + 0.2)

# Each fitted model provides NCA, SCA and joint NSCA tables.
stopifnot(identical(
  names(nsca_results(fits$HH)),
  c("necessity", "sufficiency", "necessary_and_sufficient")
))
