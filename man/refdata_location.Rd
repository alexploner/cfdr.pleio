% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refdata.R
\name{refdata_location}
\alias{refdata_location}
\title{Specify the location of the genetic reference data}
\usage{
refdata_location(refdir, vartab, chrmat, chr_order)
}
\arguments{
\item{refdir}{Character; name of the directory that holds the reference data}

\item{vartab}{Character; mame of the file that contains the main table listing
all variants in the reference data; if missing, the pre-computed
default \code{all_chr_perVariant.rds}}

\item{chrmat}{Character vector; names of the files that hold the local LD structure
as a sparse matrix, with one file per chromosome; if not specified,
the pre-computed defaults \code{chrXX_LDpairs.rds} are used, where
\code{XX} is the two-digit chromosome number.}

\item{chr_order}{Numeric vector of same length as \code{chrmat}; the chromosome
numbers (1-22) corresponding to the entries in \code{chrmat};
if missing, it is assumed that the files cover all 22
non-sex chromosomes in the given order.}
}
\value{
An object of class \code{refdata_location}. This is a list with three
entries:
* \code{refdir}: character, the directory that holds the reference data
* \code{vartab}: character, the file that holds the list of all variants in
the reference data
* \code{chrmat}: a named character vector with the files that hold the
per-chromosome LD structure; the vector names are the
chromosome numbers
}
\description{
Specify the location of the genetic reference data
}
\details{
More about the structure of the reference data in vignette FIXME.

A pre-computed set of reference data can be downloaded from FIXME; see the
instructions in the vignette or the included README file.
}
\examples{
\dontrun{
refdata_location("~/ExtraData/refdata")
}
}
