#' Calculate Fsp from Donker et al. 2017
#'
#' @inheritParams fst_fast
#'
#' @return facility x facility lower tri matrix of Fsp values
#' @export
#'
#' @examples
#' \dontrun{
#' locs <- regentrans::metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
#' pt <- regentrans::metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
#' fasta <- regentrans::aln
#' fsp <- fsp_fast(fasta = fasta, locs = locs, pt = pt)
#' }
fsp_fast = function(fasta, locs, pt){
  #source
  sourceCpp("src/bottleneck_functions.cpp")
  #check input
  check_fsp_fst_input(fasta, locs, pt)

  # If this isn't true, it will cause issues with the tibbles
  # This functionality has been moved to input checks
  #if(length(locs) != length(pt))
  #  stop("length(locs) != length(pt)")

  # Create a tibble of location data, patient ids, etc.
  data_tibble = full_join(
    tibble(
      id = names(locs),
      locs = unname(locs)
      ),
    tibble(
      id = names(pt),
      pt = unname(pt)
    ),
    by = "id"
  )

  # Which locations have multiple samples?
  locs_count = table(data_tibble$locs)
  loc_names = names(locs_count)
  names(locs_count) = NULL

  locs_with_multiple = loc_names[locs_count > 1]

  # Subset tibble to those locations
  data_to_analyze = data_tibble %>%
    filter(
      locs %in% locs_with_multiple
    ) %>%
    mutate(
      locnum = map_dbl(locs, ~ which(locs_with_multiple == .x))
    )

  # Subset fasta, and make it into a numeric matrix for usage with
  # Rcpp functions
  fasta_to_analyze = apply(fasta[data_to_analyze$id,], 2, as.numeric)

  # Get meta sequences, as well as resulting location/patient info tibble
  meta_and_locs = make_meta_seqs_fast(fasta_to_analyze,
                                       data_to_analyze)

  # Split up those objects
  fasta_sub_meta = meta_and_locs$ms_mat
  final_loc_tibble = meta_and_locs$samples_tib %>%
    mutate(
      locnum = map_dbl(locs, ~ which(locs_with_multiple == .x))
    )

  # Finally, get distance for each facility
  facil_dist = get_facil_dist(1:length(locs_with_multiple),
                              fasta_sub_meta,
                              final_loc_tibble$locnum)
  # Name resulting matrix according to actual names, not numeric
  # codes
  colnames(facil_dist) = locs_with_multiple
  rownames(facil_dist) = locs_with_multiple
  return(facil_dist)
}

make_meta_seqs_fast = function(fasta, ptl_tib){
  # Count samples from each patient
  patient_samples_counts = ptl_tib %>%
    group_by(pt, locs) %>%
    summarise(
      count = n(),
      ids = list(id),
      .groups = "keep"
    )

  # Make a reference sequence, but fast
  ref_seq = make_ref_cpp(fasta)

  # For each patient, if they have one sample, just return thatsequence;
  # if they have 2 or more, make a meta sequence
  # meta_seqs is a list of vectors of numbers corresponding to base pair
  # codes
  meta_seqs = pmap(patient_samples_counts, function(pt, locs, count, ids){
    if(count == 1){
      cur_seq = ptl_tib$id %in% ids
      return(as.numeric(fasta[cur_seq,]))
    }else{
      cur_subset = ptl_tib$id %in% ids
      return(find_major_alleles_cpp(fasta[cur_subset,], ref_seq))
    }
  })

  # Bind list of vectors into a matrix
  ms_mat = do.call(rbind, meta_seqs)

  # If patient samples have more than one id, just use
  # the first one (this is what the original code did)
  patient_samples_counts = patient_samples_counts %>%
    mutate(
      first_id = map_chr(ids, ~ .x[1])
    )

  # Not ideal, and don't need to return the whole tibble, but useful
  # for debugging
  rownames(ms_mat) = patient_samples_counts$first_id
  ret_obj = list(
    ms_mat = ms_mat,
    samples_tib = patient_samples_counts
  )
  return(ret_obj)
}
