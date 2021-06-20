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


# ------------------------------------------------------------------------------
#' updateProgress (temporarily disabled)
#'
#' update the log file to show progress
#'
#' @param country Name of country.
#' @param cnt Numeric order of countries.
#' @param runs Total number of runs.
#' @param r A numeric order of the current run.
#' @param status A numeric variable for the stage of simulation progress: \code{1}
#' -Generate the input data, \code{2}-Start to run the fortran code, \code{3}-
#' Finish the fortran code, \code{4}-Return the processed results.
#' @examples
#' updateProgress ('BGD', 1, 200, 25, 1)

updateProgress <- function (country,
                            cnt,
                            runs,
                            r,
                            status) {
  ## wait until gavi_progress is yours (so cannot be overwritten by two parallel scripts at the same time)
  # while(
  #	 !(file.rename("gavi_progress",paste0("gavi_progress_locked_",cnt,"_",r)))
  # ){}
  ## update line
  # progress <- readLines("gavi_progress_locked",-1)
  # progress[((cnt-1)*runs+r)] <- paste0(
  #	 country,
  #	 " ",
  #	 paste0(c(rep(0,(nchar(runs)-nchar(r))),r),collapse=""),
  #	 " ",
  #	 status
  # )
  # writeLines(progress,"gavi_progress_locked")
  # file.rename(paste0("gavi_progress_locked_",cnt,"_",r),"gavi_progress")

} # end of function -- updateProgress
# ------------------------------------------------------------------------------
