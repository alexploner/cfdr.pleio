## The master file for package cfdr.pleio

#' A class for calculating conditional & conjunctional fdr for pleiotropy analysis
#' @description A class that captures all data for calculating the conditional
#'   and conjuncational fdr for pleiotropy analysis, and manages the necessary
#'   analysis flow
#' 
#' @examples
#' cfdr_pleio$new()
#' @export
cfdr_pleio <- R6::R6Class("cfdr_pleio",
public = list(
  #' @field trait1_data Genetic data for the first trait
  trait1_data = NULL,
                        
  #' @field trait2_data Genetic data for the second trait
  trait2_data = NULL,
                        
  #' @field traitnames Character vector of length two, the names of the traits
  traitnames = NULL,
                        
  #' @description Initialize an empty 
  #' @param trait1 Genetic data for the first trait
  #' @param trait2 Genetic data for the second trait
  #' @param traitnames Character vector of length two, the names of the traits
  initialize = function(trait1 = NA, trait2 = NA, traitnames) {
    if (missing(traitnames)) {
      traitnames <- c("(Trait1)", "(Trait2)")
    }
    self$trait1_data <- trait1
    self$trait2_data <- trait2
    self$traitnames <- traitnames
  }
) )
