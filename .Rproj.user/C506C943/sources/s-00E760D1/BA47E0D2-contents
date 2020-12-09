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

#Check if stderr() / stdout() or some other connection to the terminal
#class_is_terminal_connection <- curry::partial(fun = class_is, args = list(tested_class = c("terminal", "connection")))



#' class_is_terminal_connection
#' 
#' Checks if an object is a connection to a terminal (e.g stderr() or stdout()).
#' 
#' @param object some object whose class you want to test
#'
#' @return TRUE/FALSE depending on whether the passed object is a connection to a terminal (logical)
#' @export
#'
#' @examples
#' class_is_terminal_connection(stdout()) #True
#' class_is_terminal_connection(stderr()) #True
#' class_is_terminal_connection("random text") #false
class_is_terminal_connection <- function(object){
  class_is(object, tested_class = c("terminal", "connection"))
}