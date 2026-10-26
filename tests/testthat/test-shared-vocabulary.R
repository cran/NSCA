# The vocabulary of condition analysis in degree now lives in SCAtools and NSCA
# delegates to it. These tests pin the interface between the two packages: the
# one lookup NSCA still keeps for itself must agree with SCAtools' copy, and
# the delegating helpers must return what they used to.

test_that("the sufficiency corner agrees with SCAtools", {
  # NSCA keeps this map so that the diagonal relationship with the necessity
  # map stays visible in one place. That is the only reason it is duplicated,
  # so the duplication has to be checked rather than trusted.
  directions <- c("HH", "LH", "HL", "LL")
  expect_equal(
    unname(NSCA:::.nsca_sufficiency_corner[directions]),
    SCAtools::sca_corner(directions)
  )
})

test_that("the necessity corner is the diagonal opposite of the sufficiency one", {
  directions <- c("HH", "LH", "HL", "LL")
  necessity <- unname(NSCA:::.nsca_necessity_corner[directions])
  sufficiency <- unname(NSCA:::.nsca_sufficiency_corner[directions])
  # Corners 1 and 4 are one diagonal, 2 and 3 the other, so opposite corners
  # sum to 5 in every direction.
  expect_equal(necessity + sufficiency, rep(5L, length(directions)))
})

test_that("the delegating helpers still answer for all four corners", {
  expect_equal(
    NSCA:::.nsca_boundary_type(1:4),
    c("ceiling", "ceiling", "floor", "floor")
  )
  expect_equal(
    NSCA:::.nsca_empty_space_text(4L),
    "High X with low Y"
  )
  expect_equal(
    NSCA:::.nsca_relationship_text("HH"),
    "Higher X for higher Y"
  )
  expect_equal(NSCA:::.nsca_validate_direction("hh ", n = 2L), c("HH", "HH"))
})

test_that("a variable name containing an upper-case Y survives substitution", {
  # The fault this merge was partly for: two successive sub(..., fixed = TRUE)
  # calls filled X first and then hit the Y inside the name just substituted.
  # SYSTEM, Yield, YEARS and QUALITY all triggered it.
  expect_equal(
    NSCA:::.nsca_empty_space_text(4L, "SYSTEM", "score"),
    "High SYSTEM with low score"
  )
  expect_equal(
    NSCA:::.nsca_relationship_text("HH", "Yield", "GDP"),
    "Higher Yield for higher GDP"
  )
  expect_equal(
    SCAtools:::.sca_empty_space_text(4L, "QUALITY", "output"),
    "High QUALITY with low output"
  )
})

test_that("the vocabulary functions refuse a corner that does not exist", {
  expect_error(SCAtools::sca_boundary_type(5), "1 through 4")
  expect_error(SCAtools::sca_corner_location(0), "1 through 4")
  expect_error(SCAtools::sca_empty_space(NA), "1 through 4")
})
