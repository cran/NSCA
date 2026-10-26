# Direction handling for necessary-and-sufficient statements.
#
# A necessary-and-sufficient statement fixes two empty corners at once. Writing the
# sufficiency claim as X -> Y and the necessity claim as Y -> X, and applying
# contraposition to each, puts the two empty regions on opposite corners of the
# same scatter plot. The necessity corner is therefore always the diagonal
# opposite of the sufficiency corner, in every direction.

.nsca_directions <- c("HH", "LH", "HL", "LL")
.nsca_direction_number <- stats::setNames(seq_along(.nsca_directions),
                                          .nsca_directions)
.nsca_equivalence_class <- c(
  HH = "HH/LL", LH = "LH/HL", HL = "LH/HL", LL = "HH/LL"
)

# Corner numbering follows SCAtools: 1 upper-left, 2 upper-right,
# 3 lower-left, 4 lower-right.
#
# The two are kept side by side deliberately, because the diagonal relationship
# between them is the whole content of this file and it is legible only when
# both are visible. The sufficiency row must agree with
# SCAtools::sca_corner(); a test pins that, rather than a second lookup
# derived at build time from a package that is loaded lazily.
.nsca_sufficiency_corner <- c(HH = 4L, LH = 3L, HL = 2L, LL = 1L)
.nsca_necessity_corner <- c(HH = 1L, LH = 2L, HL = 3L, LL = 4L)

# An empty corner at the top means the frontier bounds the data from above.
.nsca_corner_side <- function(corner) {
  ifelse(corner %in% c(1L, 2L), "upper", "lower")
}

# Vocabulary of condition analysis in degree, section 3.2 and Table 1.
#
# The corner-to-wording maps used to live here as well, character for character
# identical to SCAtools' copies. Two copies of one lookup table is two places
# to fix, and they had already drifted: the upper-case-Y substitution fault
# existed twice over, and the warning on the rounded-text fallback existed in
# one package and not the other. Since 0.5.0 SCAtools owns the vocabulary and
# exports it, and the four helpers below delegate.
#
# The delegation is by :: rather than ::: on purpose. Reaching into another
# package's internals is flagged by R CMD check and would make this file depend
# on names SCAtools is free to change; the exported functions are a stated
# interface between the two packages.
#
# Table 1 read through these corners: necessity uses a ceiling exactly when the
# outcome level is high (HH, LH) and a floor when it is low (HL, LL);
# sufficiency, sitting on the diagonally opposite corner, is the reverse.

.nsca_boundary_type <- function(corner) {
  SCAtools::sca_boundary_type(corner)
}

# Fill the variable names into the Table 1 wording.
.nsca_empty_space_text <- function(corner, x = "X", y = "Y") {
  SCAtools::sca_empty_space(corner, x, y)
}

.nsca_relationship_text <- function(direction, x = "X", y = "Y") {
  SCAtools::sca_relationship(direction, x, y)
}

.nsca_validate_direction <- function(direction, n = length(direction)) {
  SCAtools::sca_check_direction(direction, n)
}

.nsca_level <- function(letter) ifelse(letter == "H", "High", "Low")

.nsca_statement <- function(direction, x = "X", y = "Y") {
  direction <- .nsca_validate_direction(direction)
  sprintf(
    "%s %s is necessary and sufficient for %s %s",
    .nsca_level(substr(direction, 1L, 1L)), x,
    .nsca_level(substr(direction, 2L, 2L)), y
  )
}

.nsca_necessity_statement <- function(direction, x = "X", y = "Y") {
  direction <- .nsca_validate_direction(direction)
  sprintf(
    "%s %s is necessary for %s %s",
    .nsca_level(substr(direction, 1L, 1L)), x,
    .nsca_level(substr(direction, 2L, 2L)), y
  )
}

.nsca_sufficiency_statement <- function(direction, x = "X", y = "Y") {
  direction <- .nsca_validate_direction(direction)
  sprintf(
    "%s %s is sufficient for %s %s",
    .nsca_level(substr(direction, 1L, 1L)), x,
    .nsca_level(substr(direction, 2L, 2L)), y
  )
}

#' Empty corners of a necessary-and-sufficient statement
#'
#' Returns the two scatter-plot corners that a necessary-and-sufficient
#' statement implies are empty. For each direction, the necessity corner is
#' diagonally opposite the sufficiency corner.
#'
#' @param direction One or more of `"HH"`, `"LH"`, `"HL"`, `"LL"`. The first
#'   letter refers to the level of the condition and the second to the level of
#'   the outcome.
#' @return A data frame with the necessity and sufficiency corners and the side
#'   of the point cloud each frontier bounds.
#' @export
#' @examples
#' nsca_corners("HH")
#' nsca_corners(c("HH", "LL"))
nsca_corners <- function(direction = "HH") {
  direction <- .nsca_validate_direction(direction)
  necessity <- unname(.nsca_necessity_corner[direction])
  sufficiency <- unname(.nsca_sufficiency_corner[direction])
  data.frame(
    direction = direction,
    nsca_number = unname(.nsca_direction_number[direction]),
    relationship = .nsca_relationship_text(direction),
    necessity_corner = necessity,
    necessity_location = SCAtools::sca_corner_location(necessity),
    necessity_empty_space = .nsca_empty_space_text(necessity),
    necessity_boundary = .nsca_boundary_type(necessity),
    necessity_bound = .nsca_corner_side(necessity),
    sufficiency_corner = sufficiency,
    sufficiency_location = SCAtools::sca_corner_location(sufficiency),
    sufficiency_empty_space = .nsca_empty_space_text(sufficiency),
    sufficiency_boundary = .nsca_boundary_type(sufficiency),
    sufficiency_bound = .nsca_corner_side(sufficiency),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Complete mapping between directions, corners and component statements
#'
#' Reproduces Table 1 of the condition analysis in degree framework, which sets
#' out the main directional types of the X-Y relationship together with the
#' expected empty space and boundary each implies for necessity and for
#' sufficiency, and adds the corner numbering the frontier engines use.
#'
#' The direction is a claim about the theorized relationship, not about the
#' shape of the estimated frontier, which may be a step function or a straight
#' line under any of the four types. Theory may instead posit a curvilinear or
#' U-shaped relationship, for which the expected empty space is no longer a
#' single corner; such relationships are outside the scope of this package.
#'
#' @return A data frame describing all four necessary-and-sufficient
#'   directions. The `biconditional` column retains this name because it
#'   contains the joint statement written out in full, and the shorter logical
#'   term is more convenient as a column name.
#' @references
#' Chan, Y. (in press). *From Necessity to Sufficiency and
#' Necessity-and-Sufficiency: Extending Condition Analysis in Degree*.
#' University of Leeds.
#' @seealso [nsca_corners()], [nsca_terms()]
#' @export
#' @examples
#' nsca_direction_map()
nsca_direction_map <- function() {
  direction <- .nsca_directions
  corners <- nsca_corners(direction)
  data.frame(
    nsca_number = unname(.nsca_direction_number[direction]),
    direction = direction,
    relationship = .nsca_relationship_text(direction),
    equivalence_class = unname(.nsca_equivalence_class[direction]),
    # Column name kept: `biconditional` is the concise logical term for the
    # joint statement. The prose everywhere else says necessary-and-sufficient.
    biconditional = .nsca_statement(direction),
    sufficiency_component = .nsca_sufficiency_statement(direction),
    necessity_component = .nsca_necessity_statement(direction),
    necessity_empty_space = corners$necessity_empty_space,
    necessity_boundary = corners$necessity_boundary,
    sufficiency_empty_space = corners$sufficiency_empty_space,
    sufficiency_boundary = corners$sufficiency_boundary,
    necessity_corner = corners$necessity_corner,
    sufficiency_corner = corners$sufficiency_corner,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
