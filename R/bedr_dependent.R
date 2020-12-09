
#NOTES:
#NOTE: ALL calls to bedr::somefunction() must be piped to withr::with_path(new = GLOBAL_BEDTOOLS_DIRECTORY, code = ., action = "prefix")

#' Check for Bedtools Binary
#' @param directory_containing_binary_search directory to search for binaries within (will be preppended onto PATH temporarily). If NULL, will search through PATH.
#' @param silent supress all output headed to stderr()
#' @return (invisible) TRUE if bedtools binary was found, throws error if it was not
#' @export
#'
#' @examples
#' assert_that_bedtools_binary_available()
assert_that_bedtools_binary_available <- function(directory_containing_binary_search = NULL, silent=FALSE) {
  tmpfile=tempfile()
  
  if (silent) 
    sink(tmpfile, type = c("output"))
  
  if (!is.null(directory_containing_binary_search)){
    
    value=withr::with_path(new = directory_containing_binary_search, 
                     code = utilitybelt::assert_that(
                       bedr::check.binary("bedtools"),
                       msg = utilitybelt::fmterror("Bedtools binary not found. Please install (https://bedtools.readthedocs.io/en/latest/) and add to path")
                     ))
      
      }
  else {
      value = utilitybelt::assert_that(
        bedr::check.binary("bedtools"), 
        msg = utilitybelt::fmterror("Bedtools binary not found. Please install (https://bedtools.readthedocs.io/en/latest/) and add to path")
        )
      }
  
  
  if (silent) {
    sink()
    unlink(tmpfile) 
  }
    
  return(invisible(value))
}



