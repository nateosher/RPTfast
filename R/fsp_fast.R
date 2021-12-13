library(Rcpp)
library(RcppArmadillo)
library(tibble)
library(dplyr)
library(purrr)
sourceCpp("src/bottleneck_functions.cpp")


fsp_fast = function(fasta, locs, pt){
  if(length(locs) != length(pt))
    stop("length(locs) != length(pt)")

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

  locs_count = table(data_tibble$locs)
  loc_names = names(locs_count)
  names(locs_count) = NULL

  locs_with_multiple = loc_names[locs_count > 1]

  data_to_analyze = data_tibble %>%
    filter(
      locs %in% locs_with_multiple
    ) %>%
    mutate(
      locnum = map_dbl(locs, ~ which(locs_with_multiple == .x))
    )

  fasta_to_analyze = apply(fasta[data_to_analyze$id,], 2, as.numeric)

  meta_and_locs = make_meta_seqs_fast(fasta_to_analyze,
                                       data_to_analyze)

  fasta_sub_meta = meta_and_locs$ms_mat
  final_loc_tibble = meta_and_locs$samples_tib %>%
    mutate(
      locnum = map_dbl(locs, ~ which(locs_with_multiple == .x))
    )

  facil_dist = get_facil_dist(1:length(locs_with_multiple),
                              fasta_sub_meta,
                              final_loc_tibble$locnum)
  colnames(facil_dist) = locs_with_multiple
  rownames(facil_dist) = locs_with_multiple
  return(facil_dist)
}

make_meta_seqs_fast = function(fasta, ptl_tib){

    patient_samples_counts = ptl_tib %>%
      group_by(pt, locs) %>%
      summarise(
        count = n(),
        ids = list(id),
        .groups = "keep"
      )

    ref_seq = make_ref_cpp(fasta)


    meta_seqs = pmap(patient_samples_counts, function(pt, locs, count, ids){
      if(count == 1){
        cur_seq = ptl_tib$id %in% ids
        return(as.numeric(fasta[cur_seq,]))
      }else{
        cur_subset = ptl_tib$id %in% ids
        return(find_major_alleles_cpp(fasta[cur_subset,], ref_seq))
      }
    })

    ms_mat = do.call(rbind, meta_seqs)

    patient_samples_counts = patient_samples_counts %>%
      mutate(
        first_id = map_chr(ids, ~ .x[1])
      )

    # Not ideal,
    rownames(ms_mat) = patient_samples_counts$first_id
    ret_obj = list(
      ms_mat = ms_mat,
      samples_tib = patient_samples_counts
    )
    return(ret_obj)
}
