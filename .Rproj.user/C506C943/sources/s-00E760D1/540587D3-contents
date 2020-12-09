#' Do intervals within a bedfile overlap
#'
#' Checks whether intervals WITHIN a bedfile overlap. If so, it is usually a good idea to run bedtools merge
#'
#' @param bedfile_path path to bedfile (string)
#' @param force_strandedness should we consider it impossible for intervals on one strand to overlap with those on the other (logical)
#' @param header does my bedfile have a header (logical)
#' @param warn print warnings to stderr  (logical)
#' @param verbose verbose mode (logical)
#' @param cleanup_temp delete temporary files after use (logical)
#'
#' @return TRUE if redundancy is present, FALSE if each region is completely distinct (logical)
#' @export
#'
bed_intervals_overlap <- function(bedfile_path, force_strandedness = FALSE, header=FALSE,warn=TRUE, verbose=TRUE, cleanup_temp=TRUE) {
  #"QC"
  # if (is.null(bedtools_directory)
  #     bedtools_directory = GLOBAL_BEDTOOLS_DIRECTORY
  #
  #assert_that_bedtools_binary_available(directory_containing_binary_search = bedtools_directory, silent = !verbose)

  function_name = paste0("[",match.call()[[1]], "] ")
  utilitybelt::assert_that(file.exists(bedfile_path), msg = utilitybelt::fmterror(function_name, "could not find file: ", bedfile_path))
  utilitybelt::assert_that(is.logical(force_strandedness), msg = utilitybelt::fmterror(function_name, "force_strandedness argument [", force_strandedness, "] must be logical, not: ", class(force_strandedness)))
  utilitybelt::assert_that(is.logical(header), msg = utilitybelt::fmterror(function_name, "header argument [", header, "] must be logical, not: ", class(header)))

  outfile = tempfile("bedtest")

  #Convert arguments to options
  if(force_strandedness == FALSE) s=NULL else s=TRUE
  if(header == FALSE) h=NULL else h=TRUE

  #Merge bedfile
  bedtoolsr::bt.merge(i = bedfile_path, s = s, output = outfile, header = h, d = -1)
  #utilitybelt::assert_that(file.exists(outfile), msg = utilitybelt::fmterror(function_name, "bedtools merge has failed to produce an output file"))

  rows_in_infile=nrow(data.table::fread(file = bedfile_path, sep = "\t", select = 1, fill = TRUE))
  rows_in_outfile=nrow(data.table::fread(file = outfile, sep = "\t", select = 1, fill = TRUE))

  utilitybelt::assert_that(rows_in_outfile > 0, msg = utilitybelt::fmterror(function_name, "bedtools merge produced an empty file ... something has gone wrong"))

  redundancy_present.lgl = (rows_in_outfile!=rows_in_infile)

  if (warn & redundancy_present.lgl){
    message(utilitybelt::fmtwarning(function_name, "Overlapping intervals are present. I would HIGHLY reccomend running a bedtools merge BEFORE running any analysis. BED region redundancy can lead to unexpected results"))
  }

  if(cleanup_temp)
    unlink(outfile)
  else
    message(function_name, "temporary outfile can be found at: ", outfile, ". Use cleanup_temp=TRUE to automatically delete this file")

  return(redundancy_present.lgl)
}




# Randomly Shuffle Intervals ----------------------------------------------
#' Randomly Shuffle Intervals
#'
#' Takes an input bed and reshuffles intervals randomly throughout the genome (restricted to regions defined by \strong{regions} option).
#' The size of each interval is preserved, its just moved randomly to a new location.
#' This is function is a useful first step if looking to see if particular regions are enriched for some feature because it allow you to examine such features (e.g. variant density) in the original intervals vs randomly reshuffled intervals.
#'
#' @param input path to input bed containing intervals to re-shuffle (string)
#' @param regions path to bed file representing regions within which input intervals can be reshuffled (string)
#' @param genome path to genome file (string)
#' @param iterations number of times to independently reshuffle the bedfile
#' @param outdir directory to save output
#' @param force force overwrite if directory exists
#' @param seed seed for random number generation, used for reproducibility. When seed is supplied & iterations > 1, the nth iteration is seeded with \strong{seed + (n-1)} (numeric)
#' @param verbose verbose mode
#'
#' @return path to directory containing the input files (string)
#' @export
#'
bed_randomly_shuffle_intervals <- function(input, regions, genome, iterations=1L, outdir, force=FALSE, seed=NULL, verbose=FALSE){

  # Functions ---------------------------------------------------------------
  utilitybelt::library_quietload("dplyr", verbose)
  utilitybelt::library_quietload("assertthat", verbose)
  utilitybelt::library_quietload("crayon", verbose)
  utilitybelt::library_quietload("optparse", verbose)
  utilitybelt::library_quietload("bedtoolsr", verbose)

  input_bed.path = input
  input_genome.path = genome
  intervals_within_which_input_bed_intervals_can_be_placed.path = regions
  outdir = outdir
  prefix = utilitybelt::path_process(input_bed.path, extract_basename = TRUE, remove_extension = TRUE)
  iterations = iterations
  force_overwrite = force
  seed=seed

  # Assertions --------------------------------------------------------------
  utilitybelt::assert_that(is.null(seed) | is.numeric(seed), msg = utilitybelt::fmterror("Seed must be a number, or not supplied at all"))
  utilitybelt::assert_that(is.integer(as.integer(iterations)), iterations > 0, msg = utilitybelt::fmterror("Number of iterations must be an integer"))
  utilitybelt::assert_that(file.exists(input_bed.path), msg = utilitybelt::fmterror("File not found"))
  utilitybelt::assert_that(file.exists(input_genome.path), msg = utilitybelt::fmterror("File not found"))
  utilitybelt::assert_that(file.exists(intervals_within_which_input_bed_intervals_can_be_placed.path), msg = utilitybelt::fmterror("File not found"))
  utilitybelt::assert_that(assertthat::is.flag(force_overwrite), msg = utilitybelt::fmterror("force should be either TRUE or FALSE"))
  utilitybelt::assert_that(assertthat::is.string(outdir), msg = utilitybelt::fmterror("File not found"))


  # Create complement of bedfile --------------------------------------------
  # We need to generate the complement of our include regions and use it in the bt.shuffle call as otherwise reshuffled start positions will be kept inside include regions BUT end positions may overhang'
  intervals_within_which_input_bed_intervals_cannot_be_placed.path = paste0(tempdir(), "/regions_to_exclude.bed")
  if(file.exists(intervals_within_which_input_bed_intervals_cannot_be_placed.path)) unlink(intervals_within_which_input_bed_intervals_cannot_be_placed.path)
  bedtoolsr::bt.complement(i = intervals_within_which_input_bed_intervals_can_be_placed.path, g = input_genome.path, output = intervals_within_which_input_bed_intervals_cannot_be_placed.path)
  utilitybelt::assert_that(file.exists(intervals_within_which_input_bed_intervals_cannot_be_placed.path), msg = "bt.complement failed to produce an output file")

  # Create output directory -------------------------------------------------------
  utilitybelt::dir_create(outdir, force_overwrite = force_overwrite)

  # Randomly Generate Intervals ---------------------------------------------

  for (i in 1:iterations) {
    if (verbose){
      if(i%%50==0)
        message(i, "/", iterations, " iterations complete")
    }
    bedtools_call = paste0("bedtools shuffle",
                           " -i ", input_bed.path,
                           " -g ", input_genome.path,
                           " -incl ", intervals_within_which_input_bed_intervals_can_be_placed.path,
                           " -excl ", intervals_within_which_input_bed_intervals_cannot_be_placed.path,
                           " -maxTries ", "1000000",
                           ifelse(is.null(seed), yes = "", no= paste0(" -seed ", seed + i-1)),
                           " > ", outdir, "/", prefix, ifelse(iterations > 1, yes = paste0(".", i), no = ""),
                           ".bed"
                           )
    #message("Running: ", bedtools_call, "\n")
    system(bedtools_call)
    #
    #browser()
    #
    # bedtoolsr::bt.shuffle(
    #   i = input_bed.path,
    #   g = input_genome.path,
    #   incl = intervals_within_which_input_bed_intervals_can_be_placed.path,
    #   excl = intervals_within_which_input_bed_intervals_cannot_be_placed.path,
    #   allowBeyondChromEnd=FALSE,
    #   seed = seed,
    #   maxTries = "1000000",
    #   output = paste0(outdir, "/", prefix, ifelse(iterations > 1, yes = paste0(".", i), no = ""),".reshuffled.bed")
    # )
  }

  return(outdir)
}




#' Do All Intervals Lie Within Regions
#'
#' Takes two bedfiles and tests whether intervals in the first are ALL encapsulated by intervals in the second.
#' Useful function for checking that restriction of some feature to particular bed intervals worked as expected
#'
#' May also work with other bedtools compatible input formats (e.g. vcfs) although this functionality is untested. Use at your own risk.
#'
#' @param input_bed path to a bedfile whose intervals you want to test (string)
#' @param regions_bed path to a bedfile defining regions you expect will encapsulate input_bed intervals (string)
#'
#' @return TRUE if all input_bed intervals are encapsulated in regions_bed intervals. Otherwise returns FALSE. (logical)
#' @export
#'
bed_intervals_all_occur_within_regions <- function(input_bed, regions_bed){
  res=bedtoolsr::bt.subtract(a = input_bed, b = regions_bed, f = 1)
  return(nrow(res)==0)
}

#' bed_variants_per_interval
#'
#' @param bedfile path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param outfile file to write results to. If not supplied, will return an R dataframe instead(string)
#'
#' @return a file/dataframe (depends on whether or not outfile is supplied) with a bed format layout containing an extra column with thenumber of
#' @export
#'
bed_variants_per_interval <- function(bedfile, variantfile, genomefile, outfile=NULL) {

  utilitybelt::assert_that(assertthat::is.string(bedfile), assertthat::is.string(variantfile), assertthat::is.string(genomefile))

  if(assertthat::is.string(outfile)) {
    if(file.exists(outfile)) unlink (outfile)
  }
  res_df <- bedtoolsr::bt.coverage(a = bedfile, b = variantfile, counts = TRUE, g = genomefile, sorted = TRUE, output = outfile)


  if(!is.null(outfile))
   return(invisible(outfile))
  else
    return(res_df)
}

#' Density of variants in bed file
#'
#' @param bedfile path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @return variants per base / kilobase / megabase depending on 'scaling' (number)
#' @export
#'
bed_variant_density <- function(bedfile, variantfile, genomefile, scaling="kb" ){
  scaling_factors <- c(1, 1000, 1000000) %>% magrittr::set_names(c("b", "kb", "mb"))
  utilitybelt::assert_that(assertthat::is.string(scaling), scaling %in% names(scaling_factors))
  scaling_factor <- scaling_factors[scaling]
  coverage_count_df <- bed_variants_per_interval(bedfile = bedfile, variantfile =  variantfile, genomefile = genomefile)
  total_variants.n <- utilitybelt::df_sum_last_column(coverage_count_df)
  total_variants.n <- total_variants.n * scaling_factor
  total_region_size <- bed_total_region_size(coverage_count_df)
  variant_density = total_variants.n/total_region_size
  utilitybelt::assert_that(assertthat::is.number(variant_density))
  return(variant_density)
}


#' bed_variant_densities
#'
#' Takes n bed files and returns the average variant density (variants per bp) for each
#'
#' @param bedfiles path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @param bedfilepaths_to_names a function that takes a string (the bedfile path) and outputs a string that will be used as the 'name' of the returned value. Default uses the full filename (function)
#' @return named vector where values are variants per base, and names are derived from each bedfile path supplied (double)
#' @export
#' @examples
#' \dontrun{
#' bed_variant_densities(
#' bedfiles = c(input_bed,regions.bed),
#' variantfile = variants.bed,
#' genomefile = genome,
#' bedfilepaths_to_names = function(x) {
#' utilitybelt::path_process_vectorized(paths = x, extract_basename = TRUE, remove_extension = TRUE)
#' })
#' }
bed_variant_densities <- function(bedfiles, variantfile, genomefile, scaling="kb", bedfilepaths_to_names = function(x) { return(x) }){
  utilitybelt::assert_that(utilitybelt::files_exist_all(bedfiles), msg=utilitybelt::fmterror("The following bedfiles could not be found: \n", paste0(bedfiles[utilitybelt::files_exist_which(bedfiles)], collapse = "\n")))

  variant_densities.v <- purrr::map_dbl(.x = bedfiles, .f = purrr::partial(bed_variant_density, variantfile=variantfile, genomefile=genomefile, scaling=scaling))
  utilitybelt::assert_that(length(variant_densities.v) == length(bedfiles))

  variant_densities.v <- magrittr::set_names(variant_densities.v, bedfilepaths_to_names(bedfiles))
  return(variant_densities.v)
}


#' bed_variant_densities_compare
#'
#' Compares variant density of two bedfiles with a function of your choosing
#'
#' @param bedfiles path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @param bedfilepaths_to_names a function that takes a string (the bedfile path) and outputs a string that will be used as the 'name' of the returned value. Default uses the full filename (function)
#' @param comparator a function that takes two arguments. 1. the variant density of bedfile. 2. the variant density of bedfile2. The function should returns some metric of the difference. For example, magrittr::subtract would subtract the second argument from the first. (function)
#' @return result of comparator_function (variable -- usually numeric but could be boolean or anything else)
#' @export
bed_variant_densities_compare <- function(bedfiles, variantfile, genomefile, scaling="kb", bedfilepaths_to_names = function(x) { return(x) }, comparator){
  utilitybelt::assert_that(utilitybelt::files_exist_all(bedfiles), msg=utilitybelt::fmterror("The following bedfiles could not be found: \n", paste0(bedfiles[utilitybelt::files_exist_which(bedfiles)], collapse = "\n")))
  utilitybelt::assert_that(length(bedfiles)==2, msg=utilitybelt::fmterror("This function requires exactly two bedfiles, not ", length(bedfiles)))
  utilitybelt::assert_that(is.function(comparator))

  variant_densities.v <- bed_variant_densities(bedfiles = bedfiles, variantfile = variantfile, genomefile = genomefile, scaling = scaling, bedfilepaths_to_names = bedfilepaths_to_names)
  variant_difference = comparator(variant_densities.v[1], variant_densities.v[2])
  names(variant_difference) <- paste0(variant_densities.v, collapse = "_VS_")
  return(variant_difference)
}

#' bed_variant_densities_compare_divide
#'
#' Compares variant density of two bedfiles by dividing the density of the first by the second
#' Variant density is given in variants/base
#'
#' @param bedfiles path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @param bedfilepaths_to_names a function that takes a string (the bedfile path) and outputs a string that will be used as the 'name' of the returned value. Default uses the full filename (function)
#'
#' @return result of comparator_function (variable -- usually numeric but could be boolean or anything else)
#' @export
bed_variant_densities_compare_divide <- function(bedfiles, variantfile, genomefile, scaling = "kb", bedfilepaths_to_names = function(x) { return(x) }){
  comparator=function(density1, density2) { density1/density2 }
  bed_variant_densities_compare(bedfiles = bedfiles, variantfile = variantfile, genomefile = genomefile, scaling = scaling, bedfilepaths_to_names = bedfilepaths_to_names, comparator = comparator)
}

#' bed_variant_densities_compare_divide
#'
#' Compares variant density of two bedfiles by comparing the absolute difference between the two density values
#' Variant density is given in variants/base
#'
#' @param bedfiles path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @param bedfilepaths_to_names a function that takes a string (the bedfile path) and outputs a string that will be used as the 'name' of the returned value. Default uses the full filename (function)
#'
#' @return result of comparator_function (variable -- usually numeric but could be boolean or anything else)
#' @export
bed_variant_densities_compare_abs_difference <- function(bedfiles, variantfile, genomefile, scaling = "kb", bedfilepaths_to_names = function(x) { return(x) }){
  comparator=function(density1, density2) { abs(density1-density2) }
  bed_variant_densities_compare(bedfiles = bedfiles, variantfile = variantfile, genomefile = genomefile, scaling=scaling, bedfilepaths_to_names = bedfilepaths_to_names, comparator = comparator)
}

#' bed_variant_densities_compare_divide
#'
#' Compares variant density of two bedfiles by comparing the absolute difference between the two density values
#' Variant density is given in variants/base
#'
#' @param bedfiles path to bed file (string)
#' @param variantfile path to bedtools-compatible variant file, e.g. bed/vcf (string)
#' @param genomefile path to genome file (string)
#' @param scaling the unit of variant density. For example, the default 'kb' will force the function to return variants/kilobase. Options include 'b'  (base), 'kb' (kilobase) and 'mb' (megabase) (string)
#' @param bedfilepaths_to_names a function that takes a string (the bedfile path) and outputs a string that will be used as the 'name' of the returned value. Default uses the full filename (function)
#'
#' @return result of comparator_function (variable -- usually numeric but could be boolean or anything else)
#' @export
bed_variant_densities_compare_subtract <- function(bedfiles, variantfile, genomefile, scaling = "kb", bedfilepaths_to_names = function(x) { return(x) }){
  comparator=function(density1, density2) { density1-density2 }
  bed_variant_densities_compare(bedfiles = bedfiles, variantfile = variantfile, genomefile = genomefile, scaling=scaling, bedfilepaths_to_names = bedfilepaths_to_names, comparator = comparator)
}

#'
#' bed_variant_densities_compare_directory_of_beds_to_single_ref
#'
#' @param directory_containing_bed A directory containing exclusively bedfiles
#' @param output_checking_function a function that takes the output vector and returns TRUE if it is as expected and FALSE if it is not. Output of this function is fed into an assertion. (function w/ logi)
#' @inheritParams bed_variant_densities
#' @inheritParams bed_variant_densities_compare
#' @param reference_bedfile path to a bedfile whose variant density is going to be compared to each variant density of the bedfiles in 'directory_containing_bed' (string)
#'
#' @return vector of densities_compared_to_ref (usually numeric but could vary depending on comparator)
#' @export
bed_variant_densities_compare_directory_of_beds_to_single_ref <- function(directory_containing_bed, reference_bedfile, variantfile, genomefile, scaling = "kb", bedfilepaths_to_names = function(x) { return(x) }, comparator = magrittr::subtract, output_checking_function = is.numeric){
  utilitybelt::assert_that(assertthat::is.string(reference_bedfile))
  utilitybelt::assert_that(file.exists(reference_bedfile))
  utilitybelt::assert_that(length(dir(directory_containing_bed)) > 0 , msg = "Supplied directory is empty")
  #browser()
  utilitybelt::assert_that(all(file_is_valid_bed_vectorized(filepaths = dir(directory_containing_bed, full.names = TRUE), genomefile = genomefile)))

  ref_variant_density = bed_variant_density(bedfile = reference_bedfile, variantfile = variantfile, genomefile = genomefile, scaling = scaling)

  densities_compared_to_ref = bed_variant_densities(bedfiles = dir(directory_containing_bed, full.names = TRUE), variantfile = variantfile, genomefile = genomefile, scaling = scaling, bedfilepaths_to_names = bedfilepaths_to_names) %>%
    variant_density_compare_all_to_ref(variant_densities = ., reference_variant_density = ref_variant_density, comparator = comparator, output_checking_function = output_checking_function)

  return(densities_compared_to_ref)
}
#Next time::


#' Compare several variant densities to single reference
#'
#' Takes a bunch of raw variant densities + a single 'reference variant density' and applies some comparator function \strong{ f(raw_variant_density, reference_variant_density) = some_metric}
#' This comparator function is run separately for each raw_variant_density, resulting in a vector of outputs equal to length(variant_densities).
#'
#'
#' @param variant_densities vector of variant densities
#' @param reference_variant_density a single variant density you want to compare#'
#' @param comparator a function that takes two arguments. 1. the variant density of bedfile. 2. the variant density of bedfile2. The function should returns some metric of the difference. For example, magrittr::subtract would subtract the second argument from the first. (function)
#' @param output_checking_function a function that takes the output vector and returns TRUE if it is as expected and FALSE if it is not. Output of this function is fed into an assertion. (function w/ logi)
#'
#' @return for each variant density metric supplied, this function returns the output of a comparator function (output vector is equal in length to input variant_densities vector) (vector of any_value_type)
#' @export
#'
variant_density_compare_all_to_ref <- function(variant_densities, reference_variant_density, comparator, output_checking_function=is.numeric){
  utilitybelt::assert_that(is.double(variant_densities))
  utilitybelt::assert_that(assertthat::is.number(reference_variant_density))
  utilitybelt::assert_that(is.function(comparator))
  utilitybelt::assert_that(is.function(output_checking_function))
  differences.v = variant_densities %>%
    sapply(FUN = function(variant_density) { comparator(variant_density, reference_variant_density)}, USE.NAMES = TRUE)

  utilitybelt::assert_that(!is.list(differences.v))
  utilitybelt::assert_that(length(differences.v)==length(variant_densities))

  utilitybelt::assert_that(output_checking_function(differences.v))
  return(differences.v)
}
