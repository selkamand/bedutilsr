
#' Is Bedfile Sorted?
#'
#' @param bedfile_path path to bed file
#' @param genome path to bedtools genome file. To make one, see: https://www.biostars.org/p/70795/
#' @param warn Print warning messages if bed is unsorted (logical)
#'
#' @return TRUE if sorted as per genome file, FALSE if not (logical)
#' @importFrom rlang .data
#' @export
bed_is_sorted <- function(bedfile_path, genome, header = NULL, warn=TRUE) {

  options(dplyr.summarise.inform = FALSE)

  function_name = paste0("[",match.call()[[1]], "] ")

  utilitybelt::assert_that(file.exists(bedfile_path), msg = utilitybelt::fmterror(function_name, "could not find file: ", bedfile_path))
  utilitybelt::assert_that(file.exists(genome), msg = utilitybelt::fmterror(function_name, "could not find file: ", genome))
  utilitybelt::assert_that(is.logical(warn), msg = utilitybelt::fmterror(function_name, "warn argument [", warn, "] must be logical, not: ", class(warn)))

  if(is.null(header)){
    header <- file_has_header(bedfile_path)
  }

  #Read data
  bed.df <- vroom::vroom(bedfile_path, delim = "\t", col_names = c("chr", "start", "end"), col_select = 1:2, skip = as.integer(header)) %>% suppressMessages()
  chromosomes.v <- vroom::vroom(genome, delim = "\t", col_names = c("chr"), col_select = 1)[[1]] %>% suppressMessages()

  #Check chromosome order
  chrom_in_order_observed.v <- unique(bed.df$chr)
  order_of_chrom.v <- sapply(chrom_in_order_observed.v, function(x) {
    utilitybelt::assert_that(x %in% chromosomes.v, msg = utilitybelt::fmterror(function_name,"Bed file contains chromosome ", x, " which is not present in the genome file: ", genome))
    return(which(chromosomes.v==x))
    } )

  chromosomes_are_sorted <- !is.unsorted(order_of_chrom.v)
  if(chromosomes_are_sorted == FALSE){
    if (warn) message(utilitybelt::fmtwarning(function_name,"Bed file: ", bedfile_path, " is not sorted in accordance with genome file: ", genome, ". I would highly recommend running 'bedtools sort -i <bed> -g <genomefile>' before running analysis"))
    return(chromosomes_are_sorted)
  }
  else{
    chromosome_level_start_sort.df <- bed.df %>%
      dplyr::group_by(.data$chr) %>%
      dplyr::summarise(sorted=!is.unsorted(.data$start))

    if (any(chromosome_level_start_sort.df$sorted==FALSE)){
      if (warn)
        message(utilitybelt::fmtwarning(function_name,"Bed file: ", bedfile_path, " is sorted at the chromosome level but not positionally. I would highly recommend running 'bedtools sort -i <bed> -g <genomefile>' before running analysis"))
      return(FALSE)
    }
    else
      return(TRUE)
  }
}



# chromosomes_match_genomefile --------------------------------------------

#' Bedfile chromosome correctness
#'
#' @param bedfile_path path to bed file
#' @param genome path to bedtools genome file. To make one, see: https://www.biostars.org/p/70795/
#' @param warn Print warning messages if bed contains chromosomes not in the genomefile (logical)
#'
#' @return TRUE if all chromosomes in BEDfile are also present in the genome file
#' @export
bed_chromosomes_match_genome_file <- function(bedfile_path, genome, header = NULL, warn=TRUE) {

  options(dplyr.summarise.inform = FALSE)

  function_name = paste0("[",match.call()[[1]], "] ")

  utilitybelt::assert_that(file.exists(bedfile_path), msg = utilitybelt::fmterror(function_name, "could not find file: ", bedfile_path))
  utilitybelt::assert_that(file.exists(genome), msg = utilitybelt::fmterror(function_name, "could not find file: ", genome))
  utilitybelt::assert_that(is.logical(warn), msg = utilitybelt::fmterror(function_name, "warn argument [", warn, "] must be logical, not: ", class(warn)))

  if(is.null(header)){
    header=file_has_header(bedfile_path)
  }

  #Read data
  bed.df <- vroom::vroom(bedfile_path, delim = "\t", col_names = c("chr"), col_select = 1, skip = as.integer(header)) %>% suppressMessages()
  chromosomes.v <- vroom::vroom(genome, delim = "\t", col_names = c("chr"), col_select = 1)[[1]] %>% suppressMessages()

  #Check chromosome order
  chrom_in_order_observed.v <- unique(bed.df$chr)

  chrom_not_matching_genomefile.v <- chrom_in_order_observed.v[!(chrom_in_order_observed.v %in% chromosomes.v)]

  if(length(chrom_not_matching_genomefile.v) > 0){
    if(warn)
      message(utilitybelt::fmtwarning(function_name,"Bed file contains chromosomes [", paste0(chrom_not_matching_genomefile.v, collapse=","), "] which are not present in the genome file: ", genome))
    #browser()
    return(FALSE)
  }
  else
    return(TRUE)
}

#' Check if file is bed
#'
#' Check columns 2 and 3 are integers
#'
#' @param filepath path to putative bedfile
#' @param genomefile path to bedtools genome file. To make one, see: https://www.biostars.org/p/70795/
#'
#' @return (logical)
#' @export
#'
file_is_bed <- function(filepath){
  utilitybelt::assert_that(assertthat::is.string(filepath))
  utilitybelt::assert_that(file.exists(filepath))

  #browser()
  if(!file_delimiter_is(filepath = filepath, expected_delimiter = "\t")) {return(FALSE)}

  file_df <- data.table::fread(input = filepath, sep = "\t", nrows = 30)


 is.bed <- df_ncol_is_gt(file_df, n = 2) %>%
   if_true_then_run(df_col_is_integer(dataframe = file_df, n=2)) %>%
   if_true_then_run(df_col_is_integer(dataframe = file_df, n=3))

 utilitybelt::assert_that(assertthat::is.flag(is.bed))
 return(is.bed)
}

#' Is file a valid, sorted bed
#' Is file a valid, sorted bed
#'
#' @param filepath path to putative bedfile (string)
#' @param genomefile path to bedtools genome file. To make one, see: https://www.biostars.org/p/70795/ (string)
#'
#' @return (logical)
#' @export
#'
file_is_valid_bed <- function(filepath, genomefile, warn=TRUE){
  utilitybelt::assert_that(assertthat::is.string(filepath))
  utilitybelt::assert_that(file.exists(filepath))
  utilitybelt::assert_that(assertthat::is.string(genomefile))
  utilitybelt::assert_that(file.exists(genomefile))
  #browser()
  file_is_bed(filepath) %>%
    if_true_then_run(bed_chromosomes_match_genome_file(bedfile_path = filepath, genome = genomefile, warn = warn)) %>%
    if_true_then_run(bed_is_sorted(bedfile_path = filepath, genome = genomefile, warn = warn)) %>%
    if_true_then_run(!bed_intervals_overlap(bedfile_path = filepath, header = file_has_header(filepath), warn = warn)) %>%
    return()
}

#' Are files real, valid beds
#'
#' @param filepaths vector of paths to putative bedfiles (character)
#' @inheritParams  file_is_valid_bed
#'
#' @return
#' @export
#'
file_is_valid_bed_vectorized <- function(filepaths, genomefile, warn=TRUE){
  #browser()
  utilitybelt::assert_that(is.character(filepaths))
  utilitybelt::assert_that(utilitybelt::files_exist_all(filepaths))
  utilitybelt::assert_that(assertthat::is.string(genomefile))
  utilitybelt::assert_that(file.exists(genomefile))

  file_is_valid.vector <- purrr::map_lgl(filepaths, .f = function(filepath) { file_is_valid_bed(filepath = filepath, genomefile = genomefile, warn=warn) })
  utilitybelt::assert_that(length(file_is_valid.vector) == length(filepaths))
  return(file_is_valid.vector)
}

#' Guess if file has header
#'
#' uses data.table::fread to guess if file has header. TODO:::TESTTT
#'
#' @param filepath a filepath
#'
#' @return
#' @export
#'
file_has_header <- function(filepath){
  utilitybelt::assert_that(assertthat::is.string(filepath))
  utilitybelt::assert_that(file.exists(filepath))
  dataframe <- data.table::fread(input = filepath, nrows = 5)
  colnames_if_no_header <- paste0("V", 1:ncol(dataframe))
  #print(head(dataframe))
  return(!all(colnames(dataframe)==colnames_if_no_header))
}

df_ncol_map_lgl <- function(dataframe, ncol_function){
  utilitybelt::assert_that(is.data.frame(dataframe))
  utilitybelt::assert_that(is.function(ncol_function))

  res.lgl <- dataframe %>%
    ncol() %>%
    ncol_function()

  utilitybelt::assert_that(assertthat::is.flag(res.lgl))
}

df_ncol_is_gt <- function(dataframe, n){
  utilitybelt::assert_that(assertthat::is.number(n))
  df_ncol_map_lgl(dataframe = dataframe, function(ncols) { ncols > n} )
}

df_ncol_is_lt <- function(dataframe, n){
  utilitybelt::assert_that(assertthat::is.number(n))
  df_ncol_map_lgl(dataframe = dataframe, function(ncols) { ncols < n} )
}

df_col_is_integer = function(dataframe,n){
  utilitybelt::assert_that(assertthat::is.number(n))
  utilitybelt::assert_that(df_ncol_is_gt(dataframe = dataframe, n-1))

  is.numeric(dataframe[[n]]) %>%
    return()
}

if_true_then_run <- function(condition, expr){
  utilitybelt::assert_that(assertthat::is.flag(condition))

  if(condition)
    expr
  else
    return(FALSE)
}

file_delimiter_is <- function(filepath, expected_delimiter){
  utilitybelt::assert_that(assertthat::is.string(filepath))
  utilitybelt::assert_that(file.exists(filepath))
  utilitybelt::assert_that(assertthat::is.string(expected_delimiter))

  result.lgl <- reader::get.delim(filepath)==expected_delimiter
  utilitybelt::assert_that(assertthat::is.flag(result.lgl))

  return(result.lgl)
}
