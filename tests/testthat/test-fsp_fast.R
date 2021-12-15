locs <- regentrans::metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- regentrans::metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
fasta <- regentrans::aln

test_locs_4 <- locs[locs %in% c("A", "F", "H")]
test_pt_3 <- pt[names(pt) %in% names(test_locs_4)]
test_fasta_2 <- as.DNAbin(fasta[names(test_locs_4),])

fsp <- fsp_fast(fasta = test_fasta_2, locs = test_locs_4, pt = test_pt_3)

test_that("fsp_fast works", {
  #type
  expect_true(any(class(fsp) == "matrix"))
  #number of cols
  expect_true(ncol(fsp) == 3)
  #number of rows
  expect_true(nrow(fsp) == 3)
  #colnames
  expect_true(all(colnames(fsp) == c("A", "F", "H")))
  #rownames
  expect_true(all(rownames(fsp) == c("A", "F", "H")))
  #all numeric between 0, 1
  expect_true(all(fsp <= 1) & all(fsp >= 0))
  #zeros on the diagonal
  expect_true(all(diag(fsp) == 0))
  #half of the board = 0s (one upper tri or lower tri == 0 )
  expect_true(all(fsp[upper.tri(fsp)] == 0) || all(fsp[lower.tri(fsp)] == 0))
})
