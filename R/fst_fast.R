sourceCpp("src/bottleneck_functions.cpp")
source("R/fsp_fast.R")

############################
fst_fast = function(fasta, locs, pt){
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

  naming=str_c(locs_with_multiple,".prob")

  SNP1 = get_facil_dist_fst(1:length(locs_with_multiple),
                            fasta_sub_meta,
                            final_loc_tibble$locnum)
  colnames(SNP1) = naming
  SNP1 = as_tibble(SNP1)

  #Finding S_{si} which represent the sum of all proportions of allele 1 in all population at SNP location i
  SNP1$psum=rowSums(SNP1[naming],na.rm=T)
  #Finding the number of populations that do not have this SNP location i tested
  SNP1$count=rowSums(!is.na(SNP1[naming]))
  #If SNP location i is tested in one population, it is excluded as mean in all other populations is non-existent
  SNP_reduced=SNP1 %>% dplyr::filter(count!=1) # %>% dplyr::select(-all)
  #Calculating n_k for all SNP locations i
  nk=SNP_reduced$count
  #Making a variable for S_{si}
  psum=SNP_reduced$psum
  #Making matrix needed to calculate FST
  SNP_reduced2=SNP_reduced%>%  dplyr::select(-count,-psum)%>% as.matrix()
  ph=SNP_reduced2
  ph_minus1=1-SNP_reduced2

  ps=apply(SNP_reduced2,MARGIN=2,function(x){
    (psum-x)/(nk-1)
  })
  ps_minus1=1-ps
  H_P=colSums(ph*ph_minus1*ps*ps_minus1,na.rm=TRUE)
  H_Sh=colSums((ph*ph_minus1)^2,na.rm=TRUE)
  H_Ss=colSums((ps*ps_minus1)^2,na.rm=TRUE)
  H_S=(H_Sh+H_Ss)/2
  fst=(H_S-H_P)/H_S
  names(fst)=str_c(locs_with_multiple,".fst")
  return(fst)
}
