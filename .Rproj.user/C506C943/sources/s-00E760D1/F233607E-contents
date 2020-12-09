test_that("outfile looks exactly as expected", {
  input_bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")
  expected_result = system.file("/testfiles/test_input/variants_per_interval/expected_result_input_x_variants.bed", package="bedutilsr")
  outfile = system.file("testfiles/test_input/variants_per_interval",package="bedutilsr") %>%
    paste0("/actual_result_input_x_variants.bed")

  unlink(outfile)

  bed_variants_per_interval(bedfile = input_bed, variantfile = variants.bed, genomefile = genome, outfile = outfile)
  #browser()
  expect_equivalent(tools::md5sum(outfile), tools::md5sum(expected_result))

})


test_that("outfile looks exactly as expected, this time using regions.bed", {
  input_bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")
  expected_result = system.file("/testfiles/test_input/variants_per_interval/expected_result_region_x_variants.bed", package="bedutilsr")
  outfile = system.file("testfiles/test_input/variants_per_interval",package="bedutilsr") %>%
    paste0("/actual_result_region_x_variants.bed")

  unlink(outfile)

  bed_variants_per_interval(bedfile = regions.bed, variantfile = variants.bed, genomefile = genome, outfile = outfile)
  #browser()
  expect_equivalent(tools::md5sum(outfile), tools::md5sum(expected_result))

})

test_that("variants per interval works the same when returning R dataframes mode", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")
  expected_result_input = system.file("/testfiles/test_input/variants_per_interval/expected_result_input_x_variants.bed", package="bedutilsr")
  expected_result_region = system.file("/testfiles/test_input/variants_per_interval/expected_result_region_x_variants.bed", package="bedutilsr")

  actual_regions_res = bed_variants_per_interval(bedfile = regions.bed, variantfile = variants.bed, genomefile = genome)
  expected_regions_res = data.table::fread(expected_result_region)
  expect_equivalent(actual_regions_res, expected_regions_res)

  actual_input_res = bed_variants_per_interval(bedfile = input.bed, variantfile = variants.bed, genomefile = genome)
  expected_input_res = data.table::fread(expected_result_input)
  expect_equivalent(actual_input_res, expected_input_res)

})



test_that("bed variant density estimation works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")

  expect_error(bed_variant_density(bedfile = input.bed, variantfile = variants.bed, genomefile = genome, scaling="asdas"))
  expect_equivalent(bed_variant_density(bedfile = input.bed, variantfile = variants.bed, genomefile = genome, scaling="b"), 0.008)
  expect_equivalent(bed_variant_density(bedfile = regions.bed, variantfile = variants.bed, genomefile = genome, scaling="b"), 0.001230769230769)
  expect_equivalent(bed_variant_density(bedfile = input.bed, variantfile = variants.bed, genomefile = genome, scaling="kb"), 0.008*1000)
  expect_equivalent(bed_variant_density(bedfile = regions.bed, variantfile = variants.bed, genomefile = genome, scaling="kb"), 0.001230769230769*1000)
  expect_equivalent(bed_variant_density(bedfile = input.bed, variantfile = variants.bed, genomefile = genome, scaling = "mb"), 0.008*1000000)
  expect_equivalent(bed_variant_density(bedfile = regions.bed, variantfile = variants.bed, genomefile = genome, scaling = "mb"), 0.001230769230769*1000000)

})

test_that("vectorized bed variant density estimation works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")

  expect_equivalent(bed_variant_densities(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling="b"), expected = c(0.008, 0.001230769))

})

test_that("Comparison of variant density for two bedfiles by subtracting second_bed density from the first works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")

  expect_equivalent(bed_variant_densities_compare_subtract(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "b"), expected = 0.006769231)

})

test_that("Comparison of variant density for two bedfiles by division works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")

  expect_equivalent(bed_variant_densities_compare_divide(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "kb"), expected = 6.5)

})


test_that("Comparison of variant density for two bedfiles by absolute difference works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")
  #browser()
  expect_equivalent(bed_variant_densities_compare_abs_difference(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "b"), expected = 0.006769231)

  expect_equivalent(bed_variant_densities_compare_abs_difference(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "b"), expected = 0.006769231)

})

test_that("Cbed_variant_densities_compare_directory_of_beds_to_single_ref works as expected", {
  input.bed = system.file("/testfiles/test_input/variants_per_interval/input.bed", package="bedutilsr")
  regions.bed = system.file("/testfiles/test_input/variants_per_interval/regions.bed", package="bedutilsr")
  variants.bed = system.file("/testfiles/test_input/variants_per_interval/variants.bed", package="bedutilsr")
  genome = system.file("/testfiles/test_input/hs37d5.genome", package="bedutilsr")
  dir_containing_beds = system.file("testfiles/test_input/directory_containing_beds", package = "bedutilsr")
  dir_containing_beds_and_readme = system.file("testfiles/test_input/directory_containing_beds_and_invalid", package = "bedutilsr")

  #browser()
  expect_equivalent(
    bed_variant_densities_compare_directory_of_beds_to_single_ref(
      directory_containing_bed = dir_containing_beds,
      reference_bedfile = regions.bed,
      variantfile = variants.bed,
      genomefile = genome,
      scaling = "b",
      comparator = magrittr::subtract,
      output_checking_function = is.numeric),
    c(0.006769231, 0)
    )

  expect_error(
    bed_variant_densities_compare_directory_of_beds_to_single_ref(
      directory_containing_bed = dir_containing_beds_and_readme,
      reference_bedfile = regions.bed,
      variantfile = variants.bed,
      genomefile = genome,
      scaling = "b",
      comparator = magrittr::subtract,
      output_checking_function = is.numeric)
  )
  #browser()
 # expect_equivalent(bed_variant_densities_compare_abs_difference(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "b"), expected = 0.006769231)

  #expect_equivalent(bed_variant_densities_compare_abs_difference(bedfiles = c(input.bed, regions.bed), variantfile = variants.bed, genomefile = genome, scaling = "b"), expected = 0.006769231)

})


