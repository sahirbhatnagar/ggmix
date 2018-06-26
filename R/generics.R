#' Extract Random Effects
#'
#' @description Generic function for extracting the random effects. This is the same generic
#' (and same name) defined in the nlme package.
#'
#' @return a numeric vector of length n of subject-specific random effects
#'
#' @seealso \code{\link{gic.ggmix}}
#' @param object any fitted model object from which random effects estimates can be
#'   extracted. Currently supports "gic.ggmix" object
#' @param ... some methods for this generic function require additional
#'   arguments.
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname ranef
#' @export
random.effects <- function(object, ...) {
  UseMethod("ranef")
}

#' @rdname ranef
#' @method random.effects default
#' @export
random.effects.default <- function(object) {
  "Unknown class"
}

#' @rdname ranef
#' @export
ranef.default <- function(object) {
  "Unknown class"
}
