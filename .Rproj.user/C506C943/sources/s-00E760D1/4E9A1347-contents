test_that("number_of_outfiles_correct",{
  test_outdir=paste0(tempdir(), "/bedutilsr_test_out")
  valid_input_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_input.bed", package = "bedutilsr")
  valid_regions_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_regions.bed", package = "bedutilsr")
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")

  for (i in c(1, 12, 16, 50, 100)) {
    bed_randomly_shuffle_intervals(input = valid_input_bed,
                                   regions = valid_regions_bed,
                                   genome = genome_file,
                                   iterations = i,
                                   outdir = paste0(test_outdir,"/number_of_outfiles_correct_",i),
                                   force = TRUE,
                                   seed = 1)
    n_files_produced = length(dir(paste0(test_outdir, "/number_of_outfiles_correct_",i)))
    expect_equal(n_files_produced, i)

    unlink(test_outdir, recursive = TRUE)

    ##ADD test to make sure output chromosomes are all as expected given genomefile and regions bedfile
  }
})

testthat::test_that("outfiles are different to input",{
  #if(file.exists(test_outdir)) unlink(test_outdir)
  test_outdir=paste0(tempdir(), "/bedutilsr_test_out")
  valid_input_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_input.bed", package = "bedutilsr")
  valid_regions_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_regions.bed", package = "bedutilsr")
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")

  outdir_out_files_are_different_to_input <- paste0(test_outdir,"/out_files_are_different_to_input")
  bed_randomly_shuffle_intervals(input = valid_input_bed,
                                 regions = valid_regions_bed,
                                 genome = genome_file,
                                 iterations = 1,
                                 outdir = outdir_out_files_are_different_to_input,
                                 force = TRUE,
                                 seed = 1)

  output_file = dir(outdir_out_files_are_different_to_input, full.names = TRUE)
  utilitybelt::assert_that(length(output_file) == 1, msg = utilitybelt::fmterror("number of output files produced in test_that 'outfiles are different to input' is expected to be 1. Set iterations to 1"))
  outbed.df = read.csv(output_file, header = FALSE, sep = "\t")
  input.df = read.csv(valid_input_bed, header=FALSE, sep="\t")

  expect_false(all(paste(outbed.df[[1]], outbed.df[[2]], outbed.df[[3]]) == paste(input.df[[1]], input.df[[2]], input.df[[3]])))
})



test_that("Outfile contains only the expected chromosomes",{
  test_outdir=paste0(tempdir(), "/bedutilsr_test_out")
  valid_input_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_input.bed", package = "bedutilsr")
  valid_regions_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_regions.bed", package = "bedutilsr")
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")
  i=10
  chromosomes_expected_in_outfile = read.csv(valid_regions_bed, header = FALSE, sep = "\t")[[1]]


  bed_randomly_shuffle_intervals(input = valid_input_bed,
                             regions = valid_regions_bed,
                             genome = genome_file,
                             iterations = i,
                             outdir = paste0(test_outdir,"/number_of_outfiles_correct"),
                             force = TRUE,
                             seed = 1)


  chromosomes_in_outfile = dir(paste0(test_outdir, "/number_of_outfiles_correct"), full.names = TRUE) %>%
    purrr::map(function(outbed) {
      outbed_df <- read.csv(outbed, header = FALSE, sep = "\t")
      return(outbed_df[[1]])
    })

  all_chromosomes_are_as_expected <- chromosomes_in_outfile %>%
    unlist() %>%
    unique() %>%
    magrittr::is_in(chromosomes_expected_in_outfile) %>%
    all()

  expect_true(all_chromosomes_are_as_expected)

  unlink(test_outdir, recursive = TRUE)
})


testthat::test_that("outfile intervals all lie within regions",{
  #if(file.exists(test_outdir)) unlink(test_outdir)
  test_outdir=paste0(tempdir(), "/bedutilsr_test_out")
  valid_input_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_input.bed", package = "bedutilsr")
  valid_regions_bed=system.file("testfiles/test_input/randomly_shuffled_intervals/valid_regions.bed", package = "bedutilsr")
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")

  outdir_out_files_are_different_to_input <- paste0(test_outdir,"/out_intervals_all_lie_within_regions")
  bed_randomly_shuffle_intervals(valid_input_bed,
                                 regions = valid_regions_bed,
                                 genome = genome_file,
                                 iterations = 1,
                                 outdir = outdir_out_files_are_different_to_input,
                                 force = TRUE,
                                 seed = 1)

  output_file = dir(outdir_out_files_are_different_to_input, full.names = TRUE)
  utilitybelt::assert_that(length(output_file) == 1, msg = utilitybelt::fmterror("number of output files produced in test_that 'outfiles are different to input' is expected to be 1. Set iterations to 1"))
  expect_true(bed_intervals_all_occur_within_regions(input_bed = output_file, regions_bed = valid_regions_bed))

  unlink(test_outdir, recursive = TRUE)
})
