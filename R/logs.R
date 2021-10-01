# logs.R
# Functions for generating code execution notes

# ------------------------------------------------------------------------------
#' writelog
#'
#' Write a text message to a log file with the current time tag.
#'
#' @param logname Name of log file.
#' @param x Text content to be added.
#' @examples
#' writelog ("test_log", "Start generating data")
writelog <- function (logname,
                      x) {
  write (
    paste0 (
      format (Sys.time(), "%Y/%m/%d %H:%M:%S"),
      " ",
      x
    ),
    file   = logname,
    append = TRUE
  )

} # end of function -- writelog
# ------------------------------------------------------------------------------
