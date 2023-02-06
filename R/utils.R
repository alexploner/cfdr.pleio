## Generic helper and utility functions for cfdr.pleio

#' Extract locations of non-zero rows in the column of sparse matrix
#'
#' Extract the row-indices of non-zero entries for a given column of a
#' sparse matrix in compressed column format (i.e. the R default)
#'
#' @param x A sparse matrix
#' @param j The column index of interest
#'
#' @details Functionally the same as `which(x[,j])`, but more than 1E3 times
#' faster
#' @returns A numerical vector no longer than the number of rows of x
#' @examples
#' mm <- cbind( c(1, 0, 0, 0), c(0, 1, 0, 1), c(1, 0, 1, 0) )
#' mm <- Matrix::Matrix(mm, sparse = TRUE)
#'
#' which_col(mm, 1)
#' which_col(mm, 2)
#' which_col(mm, 3)
#' @export
which_col <- function(x, j)
{
  x@i[(x@p[j]+1):x@p[j+1]] + 1
}


#' Exact binomial confidence intervals for a series of experiments
#'
#' Given a series of binomial experiments, summarized as a vector of
#' successes and a vector of trials, this calculates exact confidence
#' intervals slightly faster than a wrapper for base-R `binom.test`
#'
#' @param x A vector that contains the number of successes for each experiment
#' @param n A vector that contains the number of trials per experiment
#' @param conf.level The desired confidence level; default 0.95
#'
#' @details A simple pass-through wrapper for function `binom.confint` in
#'          package `binom`.
#' @returns A matrix with one row per experiment and two columns, holding the
#'          lower and upper confidence limit
#' @seealso \code{\link[binom]{binom.confint}}
#' @examples
#' x <- c(0, 10, 240)
#' n <- c(5, 100, 500)
#' prop_ci(x, n)
#' @export
prop_ci = function(x, n, conf.level = 0.95) {
  ret <- binom::binom.confint(x = x , n = n, conf.level = conf.level,
                              method = "exact")[, c("lower", "upper")]
  as.matrix(ret)
}

#' Logit and inverse logit as matrix transformations
#'
#' These functions calculate the logit and inverse logit (sigmoid function)
#' for each element of a matrix of proportions, preserving dimension names if
#' defined.
#'
#' @param x Designed for a numerical matrix of proportions, but should work
#'          for vectors and general arrays, too.
#' @param d The dimension of the result; default is the dimension of the input `x`.
#'        Set to `NULL` to return a vector.
#' @param dn A list of dimension names for the result; default are the dimension
#'        names of the input `x`. Set to `NULL` to suppress dimension names.
#'
#' @details No checks of parameters `d` and `dn` are performed.
#' @returns A vector, matrix or array, depending on the input format and the
#'          specification for argument `d`.
logit <- function(x, d = dim(x), dn = dimnames(x)) {
  ret <- log( x / (1 - x ))
  if (!is.null(d)) dim(ret) <- d
  if (!is.null(dn)) dimnames(ret) <- dn
  ret
}

#' @rdname logit
#' @export
sigmoid <- function(x, d = dim(x), dn = dimnames(x)) {
  ret <- 1 / ( 1 + exp(-x) )
  if (!is.null(d)) dim(ret) <- d
  if (!is.null(dn)) dimnames(ret) <- dn
  ret
}

#' Convert between p-values and test statistics
#'
#' Helper functions for converting between two-sided p-values on the -log10-scale
#' and positive test statistics
#'
#' @param p A vector of p-values on the -log10-scale (i.e. zero to infinity)
#' @param z A vector of test statistics from a standard normal null distribution
#'
#' @returns A vector of test statistics or p-values, respectively
#' @examples
#' p2z( -log10(0.05) )
#' z2p( 1.96 )
#' @export
p2z <- function(p) abs( qnorm(0.5*10^(-p)) )

#' @rdname p2z
#' @export
z2p <- function(z) -log10( 2*pnorm(-abs(z)) )
