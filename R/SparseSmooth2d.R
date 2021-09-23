## SparseSmooth2d.R
##
## This is (almost) carbon copy of the matalb function SparseSmooth2d that is part
## of pleiofdr https://github.com/precimed/pleiofdr
##
## This function implements smoothing on a 2d grid, using differencing matrices
## and unrolling of long vectors which I find suggestive, but somehow I can't
## really connect it. So I'll just re-implement it in R.
##
##  *Almost* a carbon copy, as the original code has some really obvious
##  inefficiences (carrying around a vector initialized to zero and never
##  changed) that suggest to me (together with the absence of any meaningful
##  comments) that the original author(s) have cripped this code from
##  somewhere else...
##
## Note that this piece of code, as the original pleiofdr, is made available
## under the Gnu Public Licence GPL-3.

#' Smoothing of values on a 2D grid
#'
#' This implements an algorithm that accepts a 2d-grid of values (a matrix of
#' proportions for the package-specific use case) and smoothes them based on
#' some facny differencing.
#'
#' @param values  an n x p matrix of values on a grid
#' @param weights a matrix of weights, with the same dimensions
#' @param smooth_param a vector of length 2, which probably holds the a smoothing
#'              parameter for rows / columns
#' @param weight_thresh minimal value for weights, smaller values are replaced by
#'              this threshold (default 0.001)
#'
#' @returns A matrix of smoothed values
#' @examples
#' ## Small example, reasonably easy to check
#' vm <- rbind( c(1, 2, 3, 2, 5), c(2, 2, 3, 1, 4),  c(1, 4, 2, 1, 3), c(6, 3, 1, 2, 2) )
#' wm <- matrix(1, nrow = 4, ncol = 5)
#' sparseSmooth2d(vm, wm)
#'
#' ## A bigger example (close to use-case size)
#' nr <- 30
#' nc <- 100
#' set.seed(1)
#' vm <- matrix(rnorm(nr*nc), nrow = nr)
#' wm <- matrix(runif(nr*nc), nrow = nr)
#' sm <- sparseSmooth2d(vm, wm)
#' @export
sparseSmooth2d <- function(values, weights, smooth_param = c(10, 10), weight_thres = 0.001)
{
  nr <- nrow(values)
  nc <- ncol(values)

  ## Set up the basic matrix for 2nd differences along columns
  ## This is a (nr-2) x nr matrix
  D2_along_col <- colDiffs(colDiffs(diag(nr)))
  D2_along_col <- Matrix(D2_along_col, sparse = TRUE)

  ## Expand the basic matrix so that it works for an unrolled, stacked
  ## grid matrix (in natural order, i.e. by column)
  ## Turns out his is just a kinda block matrix, with the diagonal blocks
  ## being the basic 2nd diff operator
  ## This is a (nr-2)*nc x nr*nc matrix
  L1 <- .bdiag( lapply(1:nc, function(x) D2_along_col) )

  ## Now the basic matrix for second differences along rows
  ## This is a (nc-2)  x nc matrix
  D2_along_row <- colDiffs(colDiffs(diag(nc)))
  D2_along_row <- Matrix(D2_along_row, sparse = TRUE)

  ## Expand the basic matrix so that works for the unrolled, stacked grid
  ## matrix.
  ## Because the differences are to be taken along rows, and rows elements
  ## are not contiguous, we need to dilute the basic matrix with zeros
  ## to eliminate the rest of the elements between the row neighbors
  ## We do this once for the basic matrix...
  ## ... transpose to turn into by-row order
  ## FIXME: why do I have to invoke the transpose-method for sparse matrices
  ## explicitly if the package actually depends on Matrix?!
  tmp <- as.vector( Matrix::t(D2_along_row) )
  ## ... replace each element by a vector consisting of itself and the
  ##     (nr-1) zeros, forming a complete column of the original matrix
  tmp_infl <- unlist( lapply( tmp, function(x) c(x, rep(0, nr-1))) )
  ## ... fill in the inflated vector into a matrix with the same rows as the
  ##     basic operator, but do it by row to reverse the first step
  ##     This is a matrix of size (nr-2) x nr*nc
  tmp_basic <- matrix(tmp_infl, nrow = nc-2, byrow = TRUE)
  ## Now build up up the inflated matrix, by adding the the same block shifted
  ## to the right at the bottom, until we have a full size matrix
  ## .... which is (nc-2)*nr x nr*nc FIXME: NOT EFFICIENT
  L2 <- tmp_basic
  for (i in 2:nr) {
    tmp_basic <- cbind(rep(0, nc-2), tmp_basic[, -nr*nc])
    L2 <- rbind(L2, tmp_basic)
  }
  L2 <- Matrix(L2, sparse = TRUE)

  ## Stack the expanded / inflated second difference operators, using
  ## the specified scaling factors for each part
  ## This is a (nr-2)*nc + (nc-2)*nr x nr*nx matrix
  L <- rbind( sqrt(smooth_param[1]) * L1, sqrt(smooth_param[2]) * L2)

  ## And this is the full matrix of weighted, squared and mixed second order
  ## difference operators
  ## This is a nr*nc x nr*nc matrix, so sparseness is not too bad
  ## FIXME: why explicit method call necessary?!
  LtL <- Matrix::t(L) %*% L

  ## Now we start with the actual data; we unwrap & stack the value- and
  ## weight matrices; so these are now vectors of length nr*nc
  values  <- as.vector( values )
  weights <- as.vector( weights )

  ## Set up the linear equation system whose solution contains the smoothed values
  ## This is just the squared/mixed operator matrix with the (truncated)
  ## weights added along the diagonal
  ## Still a matrix of size nr*nc x nr*nc
  H <- LtL + bandSparse( n = nr*nc, k = 0, diagonals = list(pmax(weights, weight_thres)) )
  ## RHS are just the weighted values
  g <- values * weights
  ## Solve it
  values_sm <- solve(H, g, sparse = TRUE)

  ## Return the smoothed matrix, again nr x nc
  matrix(values_sm, nrow = nr, ncol = nc)
}
