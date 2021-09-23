## Code for handling reference data

#' Specify the location of the genetic reference data
#'
#' @param refdir Character; name of the directory that holds the reference data
#' @param vartab Character; mame of the file that contains the main table listing
#'                 all variants in the reference data; if missing, the pre-computed
#'                 default `all_chr_perVariant.rds`
#' @param chrmat Character vector; names of the files that hold the local LD structure
#'                 as a sparse matrix, with one file per chromosome; if not specified,
#'                 the pre-computed defaults `chrXX_LDpairs.rds` are used, where
#'                 `XX` is the two-digit chromosome number.
#' @param chr_order Numeric vector of same length as `chrmat`; the chromosome
#'                    numbers (1-22) corresponding to the entries in `chrmat`;
#'                    if missing, it is assumed that the files cover all 22
#'                    non-sex chromosomes in the given order.
#'
#' @return An object of class `refdata_location`. This is a list with three
#'         entries:
#'         * `refdir`: character, the directory that holds the reference data
#'         * `vartab`: character, the file that holds the list of all variants in
#'                     the reference data
#'         * `chrmat`: a named character vector with the files that hold the
#'                     per-chromosome LD structure; the vector names are the
#'                     chromosome numbers
#'
#' @details More about the structure of the reference data in vignette FIXME.
#'
#' A pre-computed set of reference data can be downloaded from FIXME; see the
#' instructions in the vignette or the included README file.
#'
#' @examples
#' \dontrun{
#' refdata_location("~/ExtraData/refdata")
#' }
#' @export
refdata_location <- function(refdir, vartab, chrmat, chr_order)
{
  ## Missing values
  if ( missing(vartab) )  {
    vartab <- "all_chr_perVariant.rds"
  }
  if ( missing(chr_order) ) {
    chr_order <- 1:22
  }
  if ( missing(chrmat) ) {
    chrmat <- paste0( "chr", (formatC(chr_order, flag = "0", width = 2)), "_LDpairs.rds" )
  }

  ## Build return object
  ret <- list()
  ret$refdir <- refdir
  ret$vartab <- vartab
  names(chrmat) <- as.character(chr_order)
  ret$chrmat <- chrmat
  class(ret) <- "refdata_location"

  ret
}

#' Load genetic reference files
#'
#' Given an object of class `refdata_location`, load the specified reference data
#' files (variant list, LD-structure matrices) from disk into R.
#'
#' @param x An object of class `refdata_location`
#' @param chr Integer or character, specifying the chromosome for which to load
#'            the LD-structure matrix
#'
#' @returns A character string with the path of the specified file
#'
#' @examples
#' \dontrun{
#' exref <- refdata_location("~/ExtraData/refdata")
#' get_variants(exref)
#' get_chrmat(exref, chr = 1)
#' }
#' @export
get_variants <- function(x)
{
  file.path(x$refdir, x$vartab)
}

#' @rdname get_variants
#' @export
get_chrmat <- function(x, chr)
{
  file.path(x$refdir, x$chrmat[ as.character(chr) ])
}

#' Get information on the chromosomes in a genetic reference data set
#'
#' Given an object of class `refdata_location`, report either the total number
#' of chromosomes covered by the reference data, or the vector of the individual
#' chromosome numbers.
#'
#' @param x An object of class `refdata_location`
#'
#' @returns A numeric vector with either the number of chromosomes (length one)
#'          or the individual chromosome numbers
#'
#' @examples
#' \dontrun{
#' exref <- refdata_location("~/ExtraData/refdata")
#' get_chr_num(exref)
#' get_chr_set(exref)
#' }
#' @export
get_chr_num <- function(x)
{
  length( get_chr_set(x) )
}

#' @rdname get_chr_num
#' @export
get_chr_set <- function(x)
{
  as.numeric( names( x$chrmat ))
}
