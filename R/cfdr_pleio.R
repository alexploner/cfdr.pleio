## The master file for package cfdr.pleio

#' @import R6
#' @import Matrix
#' @import data.table
#' @import matrixStats
#' @importFrom stats pnorm qnorm
#' @importFrom interp bilinear
NULL

#' A class for calculating conditional & conjunctional fdr for pleiotropy analysis
#'
#' @description A class that captures all data for calculating the conditional
#'   and conjuncational fdr for pleiotropy analysis, and manages the necessary
#'   analysis flow
#'
#' @examples
#' cfdr_pleio$new()
#' @export
#' @import R6
cfdr_pleio <- R6::R6Class("cfdr_pleio", public = list(

  #' @field trait_data Genetic data for both traits, aligned with each other and
  #' the reference data
  #' @field trait_names Character vector of length two, the names of the traits
  #' @field trait_columns Named character vector of length three, listing the
  #'        names of the required columns in the input trait data. These are
  #'        colums are:
  #'
  #'        * `id`: unique variant identifiers; default `"SNP"` (character)
  #'        * `beta`: variant-specific effect sizes (log-odds ratios or similar);
  #'                  default `"BETA"` (numeric)
  #'        * `pval`: p-value corresponding to `beta`; default `"PVAL"` (numeric).
  trait_data    = NULL,
  trait_names   = NULL,
  trait_columns = c(id = "SNP", beta = "BETA", pval = "PVAL"),

  #' @field refdat_orig An object of class `refdata_location` that points at the
  #'        original (uncompressed) genetic reference data
  #' @field refdat_local An object of class `refdata_location` that points at the
  #'        pre-processed (compressed) genetic reference data used for the current
  #'        analysis
  #' @field ref_columns Named character vector of length seven, listing the
  #'        names of the required columns in the input trait data. These columns
  #'        are:
  #'
  #'        * `id`: unique variant identifier; default `"SNP"` (character)
  #'        * `chr`: variant chromosome number; default `"CHR"` (numeric)
  #'        * `bp`: variant base pair position; default "`BP`" (numeric)
  #'        * `a1`, `a2`:  major and minor alleles at position; default `"A1"`
  #'                       and `"A2"` (character)
  #'        * `maf`: minor allele frequency; default `"MAF"` (numeric)
  #'        * `interg`: flag indicating whether the variant is intergenic (1)
  #'                    or not (0); default `"INTERGENIC"` (numeric 0/1)
  refdat_orig  = NULL,
  refdat_local = NULL,
  ref_columns  = c(id = "SNP", chr = "CHR", bp = "BP", a1 = "A1", a2 = "A2",
                   maf = "MAF", interg = "INTERGENIC"),

  #' @field index A logical matrix of size (number of variants) x (number of iterations)
  #'        This is the container that holds the multiple randomly pruned subsets
  #'        of the full data from which the cfdr is estimated. Default `NULL`
  #'        (i.e. uninitialized)
  #' @field n_iter Numerical; indicates the number of random prunings to be used
  #'        for estimating the cFDR. Default `NULL`(i.e. uninitialized)
  #' @field seed Numerical; suitable seed for setting the random number generator
  #'        before starting random pruning of data. Default `NULL`(i.e. uninitialized)
  #' @field gc_correct List; either NULL, i.e. not initialized, or with two
  #'                   elements:
  #'
  #'                      * `corrfac`: a numerical vector of length two, the actual
  #'                      correction factors applied
  #'                      * `corrfac_all`: a numerical matrix with two columns and
  #'                      `n_iter` rows, which holds the per-pruning-index correction
  #'                      factors
  index  = NULL,
  n_iter = NULL,
  seed   = NULL,
  gc_correct = list(),

  #' @field fdr_grid_par Named character vector of length 5, listing the
  #'        default values for setting up the grid over which the conditional
  #'        fdr is calculated
  #'
  #'        * `fdr_max`: maximum value for the fdr-trait axis of the grid -
  #'           larger values will be aggregated in a catch-all category; default 10
  #'        * `fdr_nbrk`: number of breaks along the fdr-trait axis, from zero
  #'           to `fdr_trait_max`; default 1001
  #'        * `fdr_thinfac`: factor controlling the thinning out of the bins
  #'           along the fdr-trait axis for smoothing and interpolation; default 10
  #'        * `cond_max`: maximum value for the conditioning trait axis of
  #'        the grid - larger values will be aggregated in a catch-all category;
  #'        default 3
  #'        * `cond_nbrk`: number of breaks along the conditioning trait
  #'        axis, from zero to `cond_trait_max`; defaults to 31.
  #' @field cfdr12 Estimated fdr for trait1 conditional on trait2; a list.
  #' @field cfdr21 Estimated fdr for trait2 conditional on trait1; a list.
  fdr_grid_par = c(fdr_max = 10, fdr_nbrk = 1001, fdr_thinfac = 10,
                   cond_max = 3, cond_nbrk = 31),
  cfdr12 = list(),
  cfdr21 = list(),


  #' @description Initialize an empty `cfdr_pleio` object
  #'
  #' @param trait1 data.table object, genetic data for the first trait
  #' @param trait2 data.table, genetic data for the second trait
  #' @param trait_names Character vector of length two, the names of the traits
  #' @param trait_columns Character vector overriding the default required
  #'        column names in the trait data; see field description.
  #' @param refdat An object of class `refdata_location` that contains
  #'        the directory and file names with the genetic reference data
  #' @param ref_columns Character vector overriding the default required
  #'        column names in the reference data `refdat`; see field description.
  #' @param local_refdat_path Character, the directory for storing the reduced
  #'        and compressed genetic reference data for the current analysis. This
  #'        directory will be created if it does not exist. If it exists and
  #'        is not empty, the content will be overwritten. FIXME: some kind of
  #'        failsafe?
  #' @param correct_GC logical; correct trait p-values for genetic correlation?
  #'        Default `TRUE`
  #' @param correct_SO logical; correct trait data for sample overlap? Default
  #'        `FALSE`, currently not implemented, stops with error if `TRUE`
  #' @param exclusion_range data.frame with three columns, specifying genomic
  #'        regions to be excluded from the analysis FIXME: detailed description
  #' @param filter_maf_min numeric, minimum require minor allele frequency;
  #'        variants below this threshold are exluded from the analysis. Default
  #'        0.005, set to `NA` to turn off frequency filtering
  #' @param filter_ambiguous logical, exlude ambiguous variants? Default `FALSE`;
  #'        currently not implemented, stops with error if `TRUE`
  #' @param filter_fisher logical; perform Fisher filtering? Default `FALSE`;
  #'        currently not implemented, stops with error if `TRUE`
  #' @param verbose logical; indicates whether to display progress; default `TRUE`
  #'
  #' @seealso \code{\link{refdata_location}}
  init_data = function(trait1, trait2, trait_names, trait_columns,
           refdat, ref_columns, local_refdat_path = "./local_refdat",
           correct_GC = TRUE, correct_SO = FALSE, exclusion_range,
           filter_maf_min = 0.005, filter_ambiguous = FALSE,
           filter_fisher = FALSE, verbose = TRUE
          ) {
    ## Check that arguments have right class
    if ( !inherits(trait1, "data.table") ) trait1 <- data.table( trait1 )
    if ( !inherits(trait2, "data.table") ) trait2 <- data.table( trait2 )
    stopifnot( inherits(refdat, "refdata_location") )

    ## Check format of exclusion ranges
    if (!missing(exclusion_range)) {
      stopifnot( is.data.frame( exclusion_range ) )
      stopifnot( ncol(exclusion_range) == 3 )
      stopifnot( all( colnames(exclusion_range) == c("chr", "from", "to") ) )
      stopifnot( all ( apply(exclusion_range, 2, is.numeric) ) )
    }

    ## Check for not yet implemented features
    stopifnot("correct_SO is not yet implemented" = !correct_SO)
    stopifnot("filter_ambiguous is not yet implemented" = !filter_ambiguous)
    stopifnot("filter_fisher is not yet implemented" = !filter_fisher)

    if (verbose) cat("Preparing trait data...\n")

    ## Set default trait names, if missing
    if (missing(trait_names)) {
      trait_names <- c("(Trait1)", "(Trait2)")
    }
    ## Currently, only default names are accepted
    stopifnot("Only default trait_columns are currently accepted" = missing(trait_columns))
    trait_columns <- self$trait_columns
    stopifnot("Only default ref_columns are currently accepted" = missing(ref_columns))
    ref_columns <- self$ref_columns

    if (verbose) cat("Preparing reference data...\n")

    ## Check the target directory for the local / compressed reference data
    if ( !dir.exists(local_refdat_path) ) {
      message("Creating directory for pre-processed reference data:", local_refdat_path)
      dir.create(local_refdat_path)
    }
    self$refdat_local <- refdata_location( local_refdat_path )

    ## Column names in trait data?
    stopifnot( all( trait_columns %in% colnames(trait1) ) )
    stopifnot( all( trait_columns %in% colnames(trait2) ) )
    ## Unique identifiers in trait data
    stopifnot( !any(duplicated( trait1[, ..trait_columns["id"] ])) )
    stopifnot( !any(duplicated( trait2[, ..trait_columns["id"] ])) )
    ## FIXME: remove (unique?) missing identifiers ??!?

    if (verbose) cat("Aligning trait and reference data...\n")

    ## FIXME: drop non-shared SNPs between traits with a warning

    ## Now, inner join; still wary of data.table notation
    self$trait_data <- merge(trait1[, ..trait_columns], trait2[, ..trait_columns],
                             by = "SNP", sort = FALSE, suffixes = c("1", "2"))
    self$trait_names <- trait_names
    ## Set standard names, convert p-values
    colnames(self$trait_data)  <- c("SNP", "BETA1", "LOG10PVAL1", "BETA2", "LOG10PVAL2")
    self$trait_data$LOG10PVAL1 <- -log10( self$trait_data$LOG10PVAL1 )
    self$trait_data$LOG10PVAL2 <- -log10( self$trait_data$LOG10PVAL2 )

    ## Set the reference data location
    self$refdat_orig <- refdat

    ## Load the reference data overview table
    reftab <- readRDS( get_variants(refdat)  )
    ## Column names in trait data?
    stopifnot( all( ref_columns %in% colnames(reftab) ) )
    ## Unique identifiers in trait data?
    stopifnot( !any(duplicated( reftab[, ..ref_columns["id"] ])) )

    ## Only keep the required columns (this will require renaming down the road)
    reftab <- reftab[ , ..ref_columns]

    ## Only keep the variants that are part of the trait data
    ## This may be horrible, but it seems to be competitive with more
    ## obscure data.table solutions
    ## FIXME: keep shared traits that are NOT in the reference (and be informative about it)
    setkey(self$trait_data, SNP)
    setkey(reftab, SNP)
    common_snps <- intersect(self$trait_data$SNP, reftab$SNP)
    self$trait_data <- self$trait_data[common_snps]
    reftab <- reftab[common_snps]

    ## Sort both trait data and reference data by location: CHR, BP (SNP as tie breaker)
    ## This is absolutely crucial, do NOT mess this up
    rnk <- frank(reftab, CHR, BP, SNP, ties.method = "first")
    self$trait_data <- self$trait_data[rnk]
    reftab <- reftab[rnk]

    if (verbose) cat("Filtering aligned data...\n")

    ## Start filtering: by default, we keep everybody
    ndx_keep <- TRUE
    ## Filter on MAF
    if ( !is.na(filter_maf_min) ) {
      ndx_keep <- ndx_keep & reftab$MAF > filter_maf_min
    }

    ## Specified exclusions here
    if ( !missing(exclusion_range) ) {
      for (i in 1:nrow(exclusion_range)) {
        ndx_keep <- ndx_keep & (
        ## Keep if....
          (reftab$CHR != exclusion_range$chr[i]) |  ## ... on wrong chromosome
          (reftab$BP < exclusion_range$from[i])  |  ## ... before the range
          (reftab$BP > exclusion_range$to[i])       ## ... after the range
        )
      }
    }

    ## Filter for ambiguous SNPs, based on reference allele information, if required
    ## (not default, but specified for Edu/Swb example)
    ##
    ## So ambiguous is when the SNP variants are the naturally complementary pairs,
    ## cos in this case, a simple strand-mixup could lead to switch that is not
    ## really a SNP (iow, exclude SNPs that are essentially just a flip between
    ## strands)... I think, though this may be more complicated
    ## https://www.snpedia.com/index.php/Ambiguous_flip
    ##
    ## Anyhoozle, this is actually the definition employed for the pre-cooked
    ## reference data, the confusing python code not withstanding; see
    ## read_matlab.R for details and verification, but this is exactly the
    ## definition in used in
    ## https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/ref9545380_1kgPhase3eur_LDr2p1.mat
    if ( !missing(filter_ambiguous) ) {
      ambi <- c("C:G", "G:C", "A:T", "T:A")
      AB   <- paste(reftab$A1, reftab$A2, sep = ":")
      ndx_keep <- ndx_keep & !( AB %in% ambi )
    }

    ## Apply the filtering to trait- and reference data
    self$trait_data <- self$trait_data[ndx_keep]
    reftab <- reftab[ndx_keep]

    ## Save the compressed local reference data
    if (verbose) cat("Saving aligned reference data to disk...\n")

    ## Break the list of variants into chromosome-size bites
    ## Note, this preserves the BP-order within chromosome
    snp_per_chr = split(reftab$SNP, reftab$CHR)

    ## Loop over the LD-type matrices, and cut them down to size
    ## Check: no SNPs in GWAS that do not appear in reference?
    ## There should not be, as we have done the alignment with the per-variant
    ## reference, but we check here anyway, cos why not
    ## FIXME: make the check of the reference data part of the definition
    ## somewhere around refdata_location
    if (verbose) pb <-txtProgressBar(min = 1, max = get_chr_num(self$refdat_orig), style = 3)
    for ( i in get_chr_set(self$refdat_orig) ) {
      if (verbose) setTxtProgressBar(pb, i)

      ## Load the sparse matrix from file
      tmp <- readRDS( get_chrmat(self$refdat_orig, i) )
      ## Match the SNPs on the current chromosome to the SNPs in the
      ## sparse matrix
      ndx <- match(snp_per_chr[[i]], colnames(tmp))
      ## Throw an error if we find an unmatched variant
      stopifnot( all(!is.na(ndx)) )
      ## ... and put the baby away
      tmp <- tmp[ndx, ndx]
      saveRDS(tmp, get_chrmat(self$refdat_local, i))

    }
    if (verbose) close(pb)

    ## Add the reduced variant table, and we're done
    saveRDS(reftab, get_variants(self$refdat_local))

    ## ... aaaand we're done
    invisible(self)
  },

  #' @description Set up the randomly pruned index matrix for cFDR estimation
  #' @param n_iter Number of random prunings
  #' @param seed Integer; random seed for pruning; if missing, the function will
  #'        randomly generate and store a seed to make the analysis reproducible
  #' @param force Logical; indicate whether to override an existing index.
  #'        Default `FALSE`, which stops with a warning if the index has been
  #'        defined.
  #' @param verbose logical; indicates whether to display progress; default `TRUE`
  initialize_pruning_index = function(n_iter, seed, force = FALSE, verbose = TRUE) {

    ## FIXME: has the object been initialized?! i.e. state check

    if (verbose) cat("Preparing random pruning...\n")

    ## Check: does index exist? Only kill if forced to
    if (is.matrix(self$index)) {
      if (!force) {
        stop("Pruning index already defined, use 'force=TRUE' to discard")
      } else {
        message("Discarding existing pruning index")
      }
    }

    ## Check & set number of iterations
    stopifnot(is.numeric(n_iter))
    stopifnot(n_iter > 0)
    stopifnot(floor(n_iter) == ceiling(n_iter))
    self$n_iter = n_iter

    ## FIXME: check that the local referece data still exists / is valid??
    ## Could be using cryptographic hash.

    ## Fill in a seed, if not specified
    if (missing(seed) ) {
      set.seed(NULL)
      seed <- as.integer(runif(1, -2*1E9, 2*1E9))
    }
    self$seed <- seed

    ## Set up the logical matrix that will hold the indices
    ## We do not add row names here, tempting as it is; depending
    ## on memory and messiness of code, I may want to revisit this
    ## decision
    self$index = matrix(TRUE, nrow = nrow( self$trait_data ), ncol = self$n_iter)

    ## Initialize counter for chromosome offset; loop over
    ## chromosomes: we fill in the index in row blocks, one per
    ## chromosome
    set.seed(self$seed)
    chr_set   <- get_chr_set(self$refdat_local)
    start_chr <- 1
    ## FIXME - check: is this order correct?!
    if (verbose) cat("Starting random pruning...\n")
    if (verbose) pb <-txtProgressBar(min = 1, max = length(chr_set), style = 3)
    for (i in chr_set) {
      if (verbose) setTxtProgressBar(pb, i)

      tmp <- readRDS( get_chrmat( self$refdat_local, i))
      chr_len <- nrow(tmp)

      ## Loop over iterations (columns in row block)
      for (j in 1:n_iter) {
        ## Generate a random permutation of the variants on the current
        ## chromosome
        ro <- sample(chr_len)

        ## Initialize the vector of keepers with true: everybody
        ## has the chance to play
        keep <- rep(TRUE, chr_len)

        ## Now loop over all variants, in random order: if the current
        ## variant is still a keeper, we set all variants in LD to
        ## not selected; if not, we continue
        for (k in 1:chr_len) {
          if (keep[ro[k]]) {
            keep[ which_col( tmp, ro[k] ) ] <- FALSE
          }
        }

        ## Put this into the full index matrix
        self$index[start_chr:(start_chr + chr_len -1), j] <- keep

      } ## end of loop across one chromosome

      ## Update the start of chromosome pointer
      start_chr <- start_chr + chr_len

    }  ## End of loop across chromosomes
    if (verbose) close(pb)

    ## Calculate a per-pruning-index correction factor for genetic correlation
    ##
    ## Note that this is based on genetic correlation across variants in the same
    ## pruning iteration, i.e. nominally not in LD (as used for the Edu/Swb
    ## example, see config.txt). Thefore, using the integenic variants seems somewhat
    ## of an overkill.
    ##
    ## Note that we use a more recent & corrected reference, namely
    ## Dadd, Weale & Lewis, Genetic Epidemiology 2009 at
    ## https://onlinelibrary.wiley.com/doi/pdf/10.1002/gepi.20379
    ##
    ## FIXME: inspect & possibly implement the "in-house" method, which takes the
    ## median over fairly non-extreme percentiles of the z-scores (instead of
    ## the z-scores) in order to calculate the correction factors

    if (verbose) cat("Calculating genomic correction factors...\n")

    ## Load the list of annotated variants for the intergenic-flag
    vartab <- readRDS( get_variants( self$refdat_local ))

    ## Calculate the test statistics corresponding to the p-values
    zscore  <-self$trait_data[, lapply(.SD, p2z), .SDcols = c("LOG10PVAL1","LOG10PVAL2")]
    gc_cf   <- matrix(0, nrow = self$n_iter, ncol = 2)
    for (i in 1:self$n_iter) {
      ## Denominator is not quite a quantile of the chi2-distribution
      ## with df=1, see Dadd et al.
      gc_cf[i, ] <- as.matrix( zscore[self$index[,i] & vartab$INTERGENIC, lapply(.SD, function(x) median(x)/sqrt(0.4549))] )
    }

    ## Summarize across iterations
    gc_comb <- apply(gc_cf, 2, median, na.rm = TRUE)
    ## Correct FIXME: nicer for data.table??
    zscore  <- sweep( zscore, 2, gc_comb, FUN = "/" )
    logpval <- apply( zscore, 2, z2p )
    ## Play back
    self$trait_data[, c("LOG10PVAL1corr", "LOG10PVAL2corr")] <- as.data.table(logpval)
    self$gc_correct = list(corrfac = gc_comb, corrfac_all = gc_cf)

    invisible(self)
  },


  #' @description This is the workhorse function that calculates the conditional
  #'              fdr estimate from the randomly pruned subsets of variants
  #' @param fdr_trait Integer indicating which for which trait the conditional
  #'                  fdr is to be calculated; must be 1 or 2.
  #' @param cond_trait Integer indicating which for which trait the conditional
  #'                  fdr is to be calculated; must be 1 or 2. Only one of
  #'                  `fdr_trait` or `cond_trait` needs to be specified.
  #' @param smooth Logical; whether to smooth the grid of crude fdr estimates;
  #'               default `TRUE`.
  #' @param adjust Logical; whether to correct the fdr estimates to enforce
  #'               monotonicity; default `TRUE`.
  #' @param correct_gc logical; indicates whether to apply a correction factor
  #'                   for genetic correlation to the test statistics underlying
  #'                   the trait p-values; default `TRUE`.
  #' @param fdr_max TBA
  #' @param fdr_nbrk TBA
  #' @param fdr_thinfac TBA
  #' @param cond_max TBA
  #' @param cond_nbrk TBA
  #' @param verbose logical; indicates whether to display progress; default `TRUE`
  #' @returns The updated `cfdr_pleio`-onject, invisibly.
  calculate_cond_fdr = function(fdr_trait, cond_trait, smooth = TRUE,
                                adjust = TRUE, correct_gc = TRUE,
                                fdr_max, fdr_nbrk, fdr_thinfac, cond_max,
                                cond_nbrk, verbose = TRUE
                              ) {

    if (verbose) cat("Preparing cfdr calculation...\n")

    ## FIXME: check inputs, setup etc.
    if ( missing(fdr_trait) ) {
      if ( missing(cond_trait) ) {
        stop("Specify either fdr_trait or cond_trait")
      } else {
        stopifnot( cond_trait %in% 1:2 )
        fdr_trait <- 3 - cond_trait
      }
    } else {
      if ( missing(cond_trait) ) {
        stopifnot( fdr_trait %in% 1:2 )
        cond_trait <- 3 - fdr_trait
      } else {
        stopifnot( fdr_trait %in% 1:2 )
        stopifnot( cond_trait %in% 1:2 )
        stopifnot( cond_trait + fdr_trait != 3 )
      }
    }
    ## Set the variable names
    fdr_trait_pname  <- paste0("LOG10PVAL", fdr_trait)
    cond_trait_pname <- paste0("LOG10PVAL", cond_trait)

    ## Switch to corrected version if required
    if (correct_gc) {
      fdr_trait_pname  <- paste0(fdr_trait_pname, "corr")
      cond_trait_pname <- paste0(cond_trait_pname, "corr")
    }

    ## FIXME: warning /check if already calculated fdr exists?

    ## Defaults; currently fixed
    stopifnot("fdr_max currently fixed at default" = missing(fdr_max))
    stopifnot("fdr_nbrk currently fixed at default" = missing(fdr_nbrk))
    stopifnot("fdr_thinfac currently fixed at default" = missing(fdr_thinfac))
    stopifnot("cond_max currently fixed at default" = missing(cond_max))
    stopifnot("cond_nbrk currently fixed at default" = missing(cond_nbrk))

    fdr_max <- self$fdr_grid_par["fdr_max"]
    fdr_nbrk <- self$fdr_grid_par["fdr_nbrk"]
    fdr_thinfac <- self$fdr_grid_par["fdr_thinfac"]
    cond_max <- self$fdr_grid_par["cond_max"]
    cond_nbrk <- self$fdr_grid_par["cond_nbrk"]

    ## Define the full breaks (with open-ended class to the right)
    fdr_breaks  <- c( seq(0, fdr_max, length = fdr_nbrk), Inf )
    cond_breaks <- c( seq(0, cond_max, length = cond_nbrk), Inf )

    #' Breaks for trait1, but not trait2, are shifted half a bin-width downward;
    #' comment in original function ind_look states "shift half-bin to imitate old codes"
    fdr_breaks = fdr_breaks - (fdr_breaks[2] - fdr_breaks[1]) / 2

    ## Set up a matrix to accumulate the fdr-tables
    ## Size is somewhat weird, due to unmotivated??? reductions in matrix size
    ## in the original code, see comment below
    fdr_table <- matrix(0, nrow = cond_nbrk - 1, ncol = round( (fdr_nbrk - 1)/ fdr_thinfac ) )
    ## We store intermediate results for diagnostics
    fdrtab_iter <- array(NA, c(nrow(fdr_table), ncol(fdr_table), self$n_iter))
    fdr_table_unadj <- fdr_table

    ## Loop over the random prunings
    if (verbose) cat("Starting cfdr calculation...\n")
    if (verbose) pb <-txtProgressBar(min = 1, max = self$n_iter, style = 3)
    for (i in (1:self$n_iter)) {
      if (verbose) setTxtProgressBar(pb, i)

      ## Extract the p-values - as vectors
      p_fdr_trait  <- ( self$trait_data[ self$index[,i], ..fdr_trait_pname ] )[[1]]
      p_cond_trait <- ( self$trait_data[ self$index[,i], ..cond_trait_pname ] )[[1]]

      ## Faster with .bincode (and same parameters), but no nice row / col category names
      fdr_bins  <- cut(p_fdr_trait, breaks = fdr_breaks, right = FALSE, include.lowest = TRUE)
      cond_bins <- cut(p_cond_trait, breaks = cond_breaks, right = FALSE, include.lowest = TRUE)

      ## Note: set factor levels explicitly if using .bincode
      tab <- table(cond_bins, fdr_bins)

      ## From the original code: cumulative counts per cell, conditional on the
      ## binned p-value for trait1, accumulating from most to least significant
      ## p-value (i.e. 0->1).
      ## The flipping of the order is because we have binned -log10(p), thereby
      ## reversing the direction from least to most significant.
      cumtab12 <- matrixStats::colCumsums( tab[nrow(tab):1, ] )[nrow(tab):1, ]

      ## Calculate the fully conditional, fully cumulative table of cell counts
      ## by repeating the row-steps from above for the columns: i.e. reversing
      ## the column order, calculating the cumulative sums, and re-reversing
      ## the columns.
      ## At the end, we have in each cell the total number of variants for which
      ## p1 and p2 are smaller than or equal to the left / top edge of the bin
      ## (taking reversed directions into account)
      cumtab <- matrixStats::rowCumsums(cumtab12[, ncol(cumtab12):1])[, ncol(cumtab12):1]

      ## Calculate the total number of variants for each bin of trait2 p-values
      ## (marginal sum), and repeat this vector for each column of the original
      ## matrix
      n <- matrix( matrixStats::rowSums2(cumtab12), nrow = nrow(cumtab12), ncol = ncol(cumtab12) )

      ## For compatibility with pleioFDR; comment in ind_look states,
      ## "trim to conform to legacy codes"
      ##
      ## FIXME: current hypothesis is that this trims the last, open-on-one side
      ## catch-all interval, as we could not assign a bin center for using it in
      ## estimation - VERIFY
      cumtab <- cumtab[ 1:(nrow(cumtab) - 1), 1:(ncol(cumtab) - 1)]
      n      <- n[ 1:(nrow(n) - 1), 1:(ncol(n) - 1)]

      ## Here we have per row the conditional cumulative distribution function
      ## F(p_1|p_2) as seen in the denominator of eq. 6 of Smeland et al,
      ## Human Genetics 2020
      F12 <- cumtab / n

      ## Thin out the column space of the distribution
      ## I guess to keep the memory down for the smoothing procedure
      fdr_ndx_thin <- seq(1, ncol(F12), by = fdr_thinfac)
      F12 <- F12[, fdr_ndx_thin]

      ## If required, smooth the empirical distribution function
      ## This is done on the log-odds scale, and uses the length of
      ## the confidence interval as weight
      ## Smoothing the individual iterations, instead of the average?
      if (smooth) {

        ## Convert cumulative proportions to log-odds
        log_odds <- logit(F12)

        ## Returns a matrix with two rows
        ## FIXME: do we need this for the dense or sparse grid??
        ## FIXME: really should F12 (thinned out version) instead of cumtab
        ##        also, fix weights matrix definition below
        ci <- prop_ci(cumtab, n)
        ## Calculate the weights as 1 / squared interval lengths
        weights <- as.vector( 1/matrixStats::rowDiffs(ci)^2 )
        weights <- matrix(weights, nrow = nrow(F12))[, fdr_ndx_thin]

        ## Should any infinite value have it made so far, we zero them out
        ndx2 <- is.infinite(log_odds)
        log_odds[ndx2] <- 0
        weights[ndx2]  <- 0

        ## Call the smoother
        log_odds <- sparseSmooth2d(log_odds, weights)

        ## Undo the log odds
        F12 <- sigmoid(log_odds)
      }

      ## Calculate the cdf F0 under the null hypothesis: this is just the
      ## uniform distribution on [0,1] on the original p-value scale, and the same
      ## regardless of whether we condition on p2 or not: F0(p1) = F0(p1|p2) = p1
      ## See again Smeland et al.
      ## So all we have to do is to transform the original breaks for -log10(p1)
      ## back to the p-value scale (thinned out)
      ## Note: these are the left edges of the p1-bins
      null_lp <- fdr_breaks[fdr_ndx_thin]
      F0 <- matrix(10^(-null_lp), nrow = nrow(F12), ncol = ncol(F12), byrow = TRUE)

      ## Here we go: per-iteration fdr table, suitable truncated
      fdr <- F0 / F12
      fdr[fdr < 0] <- 0
      fdr[fdr > 1] <- 1

      ## For inspection
      fdrtab_iter[,,i] <- fdr

      ## Accumulate the sums
      ## FIXME: also sums of squares, for (imprecise) variance estimate
      fdr_table <- fdr_table + fdr
    }
    if (verbose) close(pb)
    if (verbose) cat("Finalizing cfdr calculations...\n")

    ## Average the sums
    fdr_table <- fdr_table / self$n_iter

    ## Adjust to be monotonous, if required
    if (adjust) {

      ## Save for reference
      fdr_table_unadj <- fdr_table
      for (nc in 1:ncol(fdr_table)) {
        for (nr in 2:nrow(fdr_table)) {
          fdr_table[nr, nc] <- min(fdr_table[nr, nc], fdr_table[nr-1, nc])
        }
        if (nc > 1) {
          fdr_table[, nc] <- matrixStats::rowMins(fdr_table[, (nc-1):nc])
        }
      }
    }

    ## Make per-variant predictions by fitting the full list of
    ## p-values into the table and doing linear interpolation
    ## Note that we scale the p-values to the grid of the fdr table,
    ## which requires fiddling with bin widths and such
    ## FIXME: should be moved to interp::bilinear for handling indices etc.
    ##        careful: out-of-table statistics
    fdr_breaks_short <- fdr_breaks[fdr_ndx_thin]
    fdr_bin_width    <- fdr_breaks_short[2] - fdr_breaks_short[1]
    col_ndx <- self$trait_data[[fdr_trait_pname]] / fdr_bin_width
    col_ndx <- pmax(1, pmin(ncol(fdr_table), 1 + col_ndx))

    cond_bin_width <- cond_breaks[2] - cond_breaks[1]
    row_ndx <- self$trait_data[[cond_trait_pname]] / cond_bin_width
    row_ndx <- pmax(1, pmin(nrow(fdr_table), 1 + row_ndx))

    ## Interpolate & add to main data element
    fdr_variants <- interp::bilinear(1:nrow(fdr_table), 1:ncol(fdr_table), fdr_table, row_ndx, col_ndx)$z
    ##self$data <- self$data[, fdr12 := fdr_variants]

    ## Internal information: the final lookup-table, both adjusted and
    ## unadjusted, and the per-iteration lookup-tables
    ## Also, the relevant bin limits for heatmap plotting
    ## FIXME: what about the shifted trait1-bins here?! See above
    element_name <- paste0("cfdr", fdr_trait, cond_trait)
    assign( x = element_name,
            value = list(fdr_table = fdr_table,
                         fdr_table_unadj = fdr_table_unadj,
                         fdrtab_iter = fdrtab_iter,
                         fdr_coord = c( fdr_breaks_short, fdr_max ),
                         cond_coord = cond_breaks[1:(nrow(fdr_table)+1)],
                         fdr_variants = fdr_variants
                        ),
            envir = self
          )

    invisible( self )
  },

  #' @description This function returns the original trait data as a data.table
  #'              and adds any conditional fdrs that have been calculated as
  #'              extra columns; if both conditional fdrs have been calculated,
  #'              the function also adds the conjunctional fdr.
  #' @returns A data.table with 5-8 columns
  get_trait_results = function() {

    ret <- cbind(self$trait_data[, 1:5],
                 cfdr12 = self$cfdr12$fdr_variants,
                 cfdr21 = self$cfdr21$fdr_variants
                )
    if (ncol(ret) == 5) {
      warning("No conditional fdrs calculated - returning the original data")
    }
    if (ncol(ret) == 7) {
      ## Note the extra []: supposed to make the return value from this
      ## function print *directly*, not after the second eval
      ## See https://stackoverflow.com/questions/32988099/data-table-objects-assigned-with-from-within-function-not-printed
      ret[, conj_fdr := pmax(cfdr12, cfdr21)][]
    }
    ret
  }
) ) ## End class definition
