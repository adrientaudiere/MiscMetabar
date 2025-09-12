#' Check ggplot2 compatibility for ComplexUpset
#'
#' @description
#' Checks if the current ggplot2 version is compatible with ComplexUpset.
#' ComplexUpset may not work properly with ggplot2 >= 4.0.0.
#'
#' @return Logical. TRUE if ggplot2 < 4.0.0, FALSE otherwise.
#' @export
#' @examples
#' is_ggplot2_compatible()
is_ggplot2_compatible <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    return(FALSE)
  }
  
  ggplot2_version <- packageVersion("ggplot2")
  return(ggplot2_version < "4.0.0")
}

#' Get compatibility message for ComplexUpset functions
#'
#' @description
#' Internal function to generate consistent compatibility messages.
#'
#' @param function_name Character. Name of the function being called.
#' @return Character. Compatibility message.
#' @keywords internal
.get_upset_compatibility_message <- function(function_name = "this function") {
  ggplot2_version <- if (requireNamespace("ggplot2", quietly = TRUE)) {
    as.character(packageVersion("ggplot2"))
  } else {
    "unknown"
  }
  
  paste0(
    function_name, " requires ggplot2 < 4.0.0 for ComplexUpset compatibility. ",
    "Current ggplot2 version: ", ggplot2_version, ". ",
    "Please downgrade ggplot2 or use alternative visualization methods."
  )
}