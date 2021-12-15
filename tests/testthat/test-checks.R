#create test data
locs <- regentrans::metadata %>% dplyr::select(isolate_id, facility) %>% tibble::deframe()
pt <- regentrans::metadata %>% dplyr::select(isolate_id, patient_id) %>% tibble::deframe()
fasta <- regentrans::aln

test_locs <- locs[1:4]
test_locs_2 <- locs[1:3]
test_locs_3 <- locs[2:5]
test_pt <- pt[1:4]
test_pt_2 <- pt[2:5]
test_fasta <- as.DNAbin(fasta[names(test_locs),])

test_that("check_fsp_fst_input works", {
  #one good test
  expect_true(is.null(check_fsp_fst_input(test_fasta, test_locs, test_pt)))
  #one test missing something
  expect_error(
    check_fsp_fst_input(fasta = test_fasta, locs = test_locs),
    'argument "pt" is missing, with no default',
    fixed = TRUE
  )
  #one test where fasta isn't a DNAbin
  expect_error(
    check_fsp_fst_input(fasta = "test_fasta", locs = test_locs, pt = test_pt),
    'The fasta object must be of class DNAbin, you have supplied an object of class  character',
    fixed = TRUE
  )
  #one test where locs isn't a named list
  expect_error(
    check_fsp_fst_input(fasta = test_fasta, locs = unname(test_locs), pt = test_pt),
    'The locs and pt objects must be a named list of locations named by sample IDs',
    fixed = TRUE
  )
  #one test where pt isn't a named list
  expect_error(
    check_fsp_fst_input(fasta = test_fasta, locs = test_locs, pt = "test_pt"),
    'The locs and pt objects must be a named list of locations named by sample IDs',
    fixed = TRUE
  )
  #one test where locs, pt have different # of isolates
  expect_error(
    check_fsp_fst_input(fasta = test_fasta, locs = test_locs_2, pt = test_pt),
    'You have supplied a list of more isolates (n =  3 ) with locations than exist in your patient isolate list (n =  4 )',
    fixed = TRUE
  )
  #one test where locs, pt have different isolates
  expect_error(
    check_fsp_fst_input(fasta = test_fasta, locs = test_locs_3, pt = test_pt),
    'You have supplied a list of different isolates with locations than your patient isolate list',
    fixed = TRUE
  )
  #one where locs names don't exist the fasta
  expect_warning(
    check_fsp_fst_input(fasta = test_fasta, locs = test_locs_3, pt = test_pt_2),
    'You have provided an isolate location vector of fewer isolates than are contained in your SNV distance matrix (dists). Will subset',
    fixed = TRUE
  )
})

