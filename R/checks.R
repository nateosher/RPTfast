######################checks for RPTfast inputs#########################

# Check that the fasta file is a dna bin object
#'
#' @inheritParams fsp_fast
#'
#' @noRd
#'
check_dna_bin <- function(fasta){
  if(all(class(fasta) != "DNAbin")){
    stop(paste("The fasta object must be of class DNAbin, you have supplied an object of class ",
               class(fasta)))
  }
}

#' Check locs input to fst and fsp functions
#'
#' @inheritParams fsp_fast
#'
#' @noRd
#'
check_locs <- function(locs){
  #check that the locs object is a named vector
  #if it is not a vector or has no names then not good (check both of these in tests)
  if(!(is.vector(locs)) || is.null(names(locs))){
    stop("The locs and pt objects must be a named list of locations named by sample IDs")
  }
  #check that the locs object has at least 2 items
  if(length(locs) < 2){
    stop(paste("You have only supplied locations or patient IDs for "), length(locs),
         " isolates. Please supply a named vector of locations or patient IDs for at least 2 isolates")
  }
}

#' Check that the names of the isolates in locs actually exist in the SNV matrix
#'
#' @inheritParams fsp_fast
#'
#' @noRd
#'
check_fasta_vs_locs <- function(fasta, locs){
  #check that there are less than or equal to the number of samples in the vector than in the dists matrix
  if(length(locs) > nrow(fasta)){
    warning(paste("You have supplied a list of more isolates (n = ", length(locs),
                  ") with locations than exist in your SNV distance matrix (n = ",
                  nrow(fasta),
                  "). Will subset"))
  }
  #check that there are at least 2 dists and locs in common
  if(length(intersect(rownames(fasta), names(locs))) < 2){
    stop(paste("You have not provided locations of at least 2 isolates in your SNV distance matrix (dists). Please provide locations for at least 2 isolates in your SNV distance matrix."))
  }
  #warn if they will be subsetting??
  if(!setequal(names(locs), rownames(fasta))){
    warning("You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset")
  }
}

#' Check that the names of the isolates in locs match pt
#'
#' @inheritParams fsp_fast
#'
#' @noRd
#'
check_locs_vs_pt <- function(locs, pt){
  #check that there are equal number of samples in locs and pt
  if(length(locs) != length(pt)){
    stop(paste("You have supplied a list of more isolates (n = ", length(locs),
                  ") with locations than exist in your patient isolate list (n = ",
                  length(pt),
                  ")"))
  }
  #check that the names match
  if(!all(names(locs) %in% names(pt))){
    stop(paste("You have supplied a list of different isolates with locations than your patient isolate list"))
  }
}

#' Check input to fsp_fast, fst_fast
#'
#' @inheritParams fsp_fast
#'
#' @noRd
#'
check_fsp_fst_input <- function(fasta, locs, pt){
  #check everything that is common (aka no pt) first
  #check that the dists object is the snv object returned by dist.dna
  check_dna_bin(fasta)
  #check that the locs object is a named vector
  check_locs(locs)
  #check that the pt object is a named vector
  check_locs(pt)
  #check that locs, pt have the same isolates
  check_locs_vs_pt(locs, pt)
  #check that the locs names exist in the dists dataframe
  check_fasta_vs_locs(fasta, locs)
}

