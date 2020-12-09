test_that("We can accurately detect sorted vs unsorted bedfiles", {
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")
  sorted_bed = system.file("testfiles/test_input/sorted.bed",package = "bedutilsr")
  unsorted_bed_chromosome_level = system.file("testfiles/test_input/unsorted.chrlevel.bed",package = "bedutilsr")
  unsorted_bed_position_level = system.file("testfiles/test_input/unsorted.position_level.bed",package = "bedutilsr")

  expect_true(bed_is_sorted(bedfile_path = sorted_bed, genome = genome_file, warn = FALSE))
  expect_false(bed_is_sorted(bedfile_path = unsorted_bed_chromosome_level, genome = genome_file, warn = FALSE))
  expect_false(bed_is_sorted(bedfile_path = unsorted_bed_position_level, genome = genome_file, warn = FALSE))


})

test_that("We can identify whether bedfiles contain chromosomes not in genome file", {
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")
  valid_chromosomes_sorted_bed = system.file("testfiles/test_input/sorted.bed",package = "bedutilsr")
  invalid_chromosomes_sorted_bed = system.file("testfiles/test_input/sorted.invalid_chromosomes.bed",package = "bedutilsr")
  invalid_chromosomes_sorted_2_bed = system.file("testfiles/test_input/sorted.invalid_chromosomes.2.bed",package = "bedutilsr")

  expect_true(bed_chromosomes_match_genome_file(bedfile_path = valid_chromosomes_sorted_bed, genome = genome_file, warn = FALSE))
  expect_false(bed_chromosomes_match_genome_file(bedfile_path = invalid_chromosomes_sorted_bed, genome = genome_file, warn = FALSE))
  expect_false(bed_chromosomes_match_genome_file(bedfile_path = invalid_chromosomes_sorted_2_bed, genome = genome_file, warn = FALSE))
})

test_that("bedtools binary packaged with tool can be found: ", {
  expect_true(assert_that_bedtools_binary_available(GLOBAL_BEDTOOLS_DIRECTORY, silent = TRUE))
  })


test_that("bed_intervals_overlap works", {
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")
  overlapping_intervals_bed = system.file("testfiles/test_input/overlapping_intervals.bed",package = "bedutilsr")
  overlapping_intervals_with_header_bed = system.file("testfiles/test_input/overlapping_intervals_with_header.bed",package = "bedutilsr")
  non_overlapping_intervals_bed = system.file("testfiles/test_input/non_overlapping_intervals.bed",package = "bedutilsr")
  non_overlapping_intervals_with_header_bed = system.file("testfiles/test_input/non_overlapping_intervals_with_header.bed",package = "bedutilsr")
  unsorted_overlapping_intervals_bed = system.file("testfiles/test_input/unsorted_overlapping_intervals.bed",package = "bedutilsr")

  expect_true(bed_intervals_overlap(bedfile_path = overlapping_intervals_bed, force_strandedness = FALSE, header = FALSE, warn = FALSE, verbose = FALSE, cleanup_temp = TRUE))
  expect_true(bed_intervals_overlap(bedfile_path = overlapping_intervals_with_header_bed, force_strandedness = FALSE, header = TRUE, warn = FALSE, verbose = FALSE, cleanup_temp = TRUE))
  expect_false( bed_intervals_overlap(bedfile_path = non_overlapping_intervals_bed, force_strandedness = FALSE, header = FALSE, warn = FALSE, verbose = FALSE, cleanup_temp = TRUE))
  expect_false( bed_intervals_overlap(bedfile_path = non_overlapping_intervals_with_header_bed, force_strandedness = FALSE, header = TRUE, warn = FALSE, verbose = FALSE, cleanup_temp = TRUE))
  expect_error(bed_intervals_overlap(bedfile_path = unsorted_overlapping_intervals_bed, force_strandedness = FALSE, header = FALSE, warn = FALSE, verbose = FALSE, cleanup_temp = TRUE))
  })

test_that("is_bed works", {
  genome_file = system.file("testfiles/test_input/hs37d5.genome",package = "bedutilsr")
  sorted_bed = system.file("testfiles/test_input/sorted.bed",package = "bedutilsr")
  unsorted_bed_chromosome_level = system.file("testfiles/test_input/unsorted.chrlevel.bed",package = "bedutilsr")
  unsorted_bed_position_level = system.file("testfiles/test_input/unsorted.position_level.bed",package = "bedutilsr")
  invalid_chromosomes_sorted_bed = system.file("testfiles/test_input/sorted.invalid_chromosomes.bed",package = "bedutilsr")
  invalid_chromosomes_sorted_2_bed = system.file("testfiles/test_input/sorted.invalid_chromosomes.2.bed",package = "bedutilsr")
  overlapping_intervals_bed = system.file("testfiles/test_input/overlapping_intervals.bed",package = "bedutilsr")
  overlapping_intervals_with_header_bed = system.file("testfiles/test_input/overlapping_intervals_with_header.bed",package = "bedutilsr")
  non_overlapping_intervals_bed = system.file("testfiles/test_input/non_overlapping_intervals.bed",package = "bedutilsr")
  non_overlapping_intervals_with_header_bed = system.file("testfiles/test_input/non_overlapping_intervals_with_header.bed",package = "bedutilsr")
  unsorted_overlapping_intervals_bed = system.file("testfiles/test_input/unsorted_overlapping_intervals.bed",package = "bedutilsr")

  expect_true(file_is_valid_bed(filepath = sorted_bed, genomefile = genome_file))
  expect_true(file_is_valid_bed(filepath = non_overlapping_intervals_bed, genomefile = genome_file))
  expect_true(file_is_valid_bed(filepath = non_overlapping_intervals_with_header_bed, genomefile = genome_file))

  expect_false(file_is_valid_bed(filepath = unsorted_bed_chromosome_level, genomefile = genome_file))
  expect_false(file_is_valid_bed(filepath = unsorted_bed_position_level, genomefile = genome_file))
  expect_false(file_is_valid_bed(filepath = overlapping_intervals_bed, genomefile = genome_file))
})
