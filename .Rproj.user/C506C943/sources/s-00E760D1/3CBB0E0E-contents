
#' bed_total_region_size
#'
#' @description
#' R equivalent of:
#' awk '{sum += $3-$2} END{print sum}'
#'
#' WARING: Any redundant segments of your bed-style dataframe will be counted multiple times.
#'
#' @param data a bed-style dataframe where 2nd column is start and third is stop
#'
#' @return summed length of all intervals
#' @export
bed_total_region_size <- function(data){
  utilitybelt::assert_that(is.data.frame(data))
  utilitybelt::assert_that(ncol(data) >= 3)
  utilitybelt::assert_that(is.numeric(data[[2]]), is.numeric(data[[3]]))

  #Add assertion that intervals don't overlap
  result=sum(data[[3]]-data[[2]])

  return(result)
}


