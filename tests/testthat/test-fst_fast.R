locs <- regentrans::metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- regentrans::metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
fasta <- regentrans::aln

test_locs_4 <- locs[locs %in% c("A", "F", "H")]
test_pt_3 <- pt[names(pt) %in% names(test_locs_4)]
test_fasta_2 <- fasta[names(test_locs_4),]

fst <- fst_fast(fasta = test_fasta_2, locs = test_locs_4, pt = test_pt_3)

test_that("fst_fast works", {
  #class
  expect_true(class(fst) == "numeric")
  #length
  expect_true(length(fst) == 3)
  #names
  expect_true(all(names(fst) == c("A.fst", "F.fst", "H.fst")))
  #values between 0 and 1
  expect_true(all(fst <= 1) && all(fst >= 0))
})
