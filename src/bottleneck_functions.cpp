#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double get_allele_freq_btwn_cpp(colvec alleles,
                                colvec locations,
                                int allele, int location){
  colvec ar = alleles.rows(find(alleles == allele && locations == location));
  double numer = ar.n_elem;
  colvec lr = locations.rows(find(locations == location));
  double denom = lr.n_elem;
  if(numer > denom){
    throw std::out_of_range("GAFB: numer > denom");
  }
  return  numer/denom;
}


// [[Rcpp::export]]
mat get_facil_dist(arma::colvec locs_unique,
                   arma::mat fasta_sub_meta,
                   arma::colvec sample_locs) {
  // Intentionally starting at 1, since diagonals are all 0
  // Further, note I am only filling in bottom diagonal

  // Final matrix for storing distances
  mat final_dists(locs_unique.n_elem, locs_unique.n_elem, fill::zeros);

  for(int i = 1; i < locs_unique.n_elem; i++){
    for(int j = 0; j < i; j++){
      // CONSTRUCTING SUB-MATRIX
      int row_loc = locs_unique(i);
      int col_loc = locs_unique(j);

      mat cur_submat = fasta_sub_meta.rows(find(sample_locs == row_loc ||
        sample_locs == col_loc));
      mat cur_rowloc_submat = fasta_sub_meta.rows(find(sample_locs == row_loc));
      mat cur_colloc_submat = fasta_sub_meta.rows(find(sample_locs == col_loc));
      colvec cur_submat_locs = sample_locs.rows(find(sample_locs == row_loc ||
        sample_locs == col_loc));

      // To store variance at each site
      vec between(cur_submat.n_cols, fill::zeros);
      vec within_l1(cur_submat.n_cols, fill::zeros);
      vec within_l2(cur_submat.n_cols, fill::zeros);

      for(int k = 0; k < cur_submat.n_cols; k++){
        // Get alleles for this column
        colvec cur_sm_col = cur_submat.col(k);
        colvec cur_rl_col = cur_rowloc_submat.col(k);
        colvec cur_cl_col = cur_colloc_submat.col(k);
        colvec unique_alleles = unique(cur_sm_col);
        colvec unique_locs = unique(cur_submat_locs);

        // If there are two alleles at the location:
        if(unique_alleles.n_rows == 2){
          double a1l1 = get_allele_freq_btwn_cpp(cur_sm_col,
                                                 cur_submat_locs,
                                                 unique_alleles(0),
                                                 unique_locs(0));
          double a1l2 = get_allele_freq_btwn_cpp(cur_sm_col,
                                                 cur_submat_locs,
                                                 unique_alleles(0),
                                                 unique_locs(1));
          double a2l1 = get_allele_freq_btwn_cpp(cur_sm_col,
                                                 cur_submat_locs,
                                                 unique_alleles(1),
                                                 unique_locs(0));
          double a2l2 = get_allele_freq_btwn_cpp(cur_sm_col,
                                                 cur_submat_locs,
                                                 unique_alleles(1),
                                                 unique_locs(1));

          between(k) = a1l1 * a1l2 * a2l1 * a2l2;
          within_l1(k) = (a1l1 * a2l1) * (a1l1 * a2l1);
          within_l2(k) = (a1l2 * a2l2) * (a1l2 * a2l2);
        }
      }

      double btwn_sum = sum(between);
      double within_l1_sum = sum(within_l1);
      double within_l2_sum = sum(within_l2);
      final_dists(i,j) = (((within_l1_sum + within_l2_sum) / 2) - btwn_sum) /
        ((within_l1_sum + within_l2_sum) / 2);
    }
  }
  return final_dists;
}

// [[Rcpp::export]]
NumericVector make_ref_cpp(NumericMatrix fasta){
  NumericVector final_ref(fasta.nrow());
  // Hashmap to make counting easier (i.e. can just use array)
  std::unordered_map<int, int> m = {
    {72,0},
    {40,1},
    {136,2},
    {24,3},
    {4,4},
    {240,5}
  };
  // Reverse map- just an array
  int rev_map[6] = {72, 40, 136, 24, 4, 240};

  for(int i = 0; i < final_ref.length(); i++){
    int cur_counts[6] = {0, 0, 0, 0, 0, 0};
    for(int j = 0; j < fasta.nrow(); j++){
      cur_counts[m[fasta(j,i)]] += 1;
    }

    // Are any nonzero besides the reference?
    int which_max = 0;
    int cur_max = cur_counts[0];
    for(int j = 1; j < 6; j++){
      if(cur_counts[j] > cur_max){
        which_max = j;
        cur_max = cur_counts[j];
      }
    }

    final_ref[i] = rev_map[which_max];
  }

  return final_ref;
}

// [[Rcpp::export]]
NumericVector find_major_alleles_cpp(NumericMatrix fasta, NumericVector ref){
  NumericVector final_ref(ref.length());
  // Hashmap to make counting easier (i.e. can just use array)
  std::unordered_map<int, int> m = {
    {72,0},
    {40,1},
    {136,2},
    {24,3},
    {4,4},
    {240,5}
  };
  // Reverse map- just an array
  int rev_map[6] = {72, 40, 136, 24, 4, 240};

  for(int i = 0; i < final_ref.length(); i++){
    int cur_counts[6] = {0, 0, 0, 0, 0, 0};
    int cur_ref_index = m[ref[i]];
    for(int j = 0; j < fasta.nrow(); j++){
      cur_counts[m[fasta(j,i)]] += 1;
    }

    // Are any nonzero besides the reference?
    int which_nz = -1;
    for(int j = 0; j < 6; j++){
      if(j != cur_ref_index && cur_counts[j] > 0){
        which_nz = j;
        break;
      }
    }

    if(which_nz >= 0){
      final_ref[i] = rev_map[which_nz];
    }else{
      // I think this should technically be whichever entry is nonzero,
      // but this works too
      final_ref[i] = ref[i];
    }
  }

  return final_ref;
}
