#' class_is
#'
#' Check if object has a particular class
#'
#' @param object object whose class you want to check (object)
#' @param tested_class class (string)
#'
#' @return TRUE if object class matches tested_class. FALSE if not.
#' @export
#'
class_is <- function(object, tested_class){
  utilitybelt::assert_that(is.character(tested_class))
  all(class(object) == tested_class)
}

#Partial application of function to create variations. Trying to find a way to do this compatible with roxygen
class_is_terminal_connection <- curry::partial(fun = class_is, args = list(tested_class = c("terminal", "connection")))