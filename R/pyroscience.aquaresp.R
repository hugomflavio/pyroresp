#' Convert Respirometry Data from PyroScience and AquaResp Software to the FishResp Format
#'
#' The function is used to convert raw data from 'Pyro Oxygen Logger' or 'Pyroscience Workbench' (\href{https://www.pyro-science.com}{PyroScience}) and a summary file from 'AquaResp' (\href{www.aquaresp.com}{free software}) to 'FishResp' format. This function should be applied before usage of the functions \code{\link{import.test}} and \code{\link{import.meas}}. The output is a file containing raw respirometry data in the 'FishResp' format (see Details in \code{\link{import.test}} to read more information about the 'FishResp' format)
#'
#' @param input  a dataframe or name of a file which contains raw data obtained from the 'Pyro Oxygen Logger' software (\href{https://www.pyro-science.com}{PyroScience})
#' @param pyro.type The type of pyroscience software used to log the data (one of "logger" or "workbench")
#' @param aquaresp  a dataframe or name of a file which contains summary data obtained from the 'AquaResp' software (\href{www.aquaresp.com}{free software})
#' @param n.chamber  integer: the number of chambers used in an experiment (including empty ones)
#' @param date.format  string: date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' @param wait.phase  integer: duration of the wait phase (in seconds), see the 'AquaResp' summary file (row #5)
#' @param measure.phase  integer: duration of the measure phase (in seconds), see the 'AquaResp' summary file (row #6)
#'
#' @return The function returns a data frame containing raw data in the 'FishResp' format
#'
#' @examples
#'
#' pyroscience.path = system.file("extdata/pyroscience/pyroscience.txt",
#'                  package = "FishResp")
#' aquaresp.path = system.file("extdata/pyroscience/pyroscience-aquaresp.txt",
#'                  package = "FishResp")
#'
#' pyroscience.aquaresp(file = pyroscience.path,
#'                      pyro.type = "logger",
#'                      aquaresp = aquaresp.path,
#'                      date.format = "%m/%d/%Y",
#'                      n.chamber = 1,
#'                      wait.phase = 120,
#'                      measure.phase = 600)
#'
#' @export
pyroscience.aquaresp <- function(input,
                                 pyro.type = c("logger", "workbench"),
                                 aquaresp,
                                 n.chamber = 1:4,
                                 date.format,
                                 wait.phase = NA, measure.phase = NA){

  n.chamber <- match.arg(n.chamber)

  if (is.character(input)) {
    if (pyro.type == "logger")
      pyro.data <- load.pyroscience.logger.file(input, n.chamber, date.format)
    if (pyro.type == "workbench")
      pyro.data <- load.pyroscience.workbench.file(input, date.format)
  }
  else
    pyro.data <- check.pyrocience.data(input)

  if (is.character(aquaresp))
    aquaresp.data <- load.aquaresp.file(aquaresp, date.format)

  # Join datasets
  m <- 1 #counter

  for(i in aquaresp.data$Date.Time){
    a <- which(pyro.data$Date.Time == i)
    if(length(a) == 0){
      a1 <- substring(i, 1, 18)
      a2 <- as.numeric(substring(i, 19))
      if(a2 >= 9){
        a2 <- a2 - 1
        }
      else{
        a2 <- a2 + 1
        }
      b <- paste(a1, a2, sep = "")
      a <- which(pyro.data$Date.Time == b)
      }
    if(length(pyro.data$Date.Time) - a > measure.phase){
      pyro.data$Phase[a:(a + (measure.phase-1))] <- paste("M", m, sep = "")
      }
    else{
      }
    if(a - wait.phase > 0){
      pyro.data$Phase[(a - wait.phase):(a-1)] <- paste("W", m, sep = "")
      }
    else{
      }
      m = m + 1
    }

  return(pyro.data)
}

#' Load a file containing summary phase data obtained from the 'AquaResp' software
#' 
#' @param file The name of a file which contains summary data obtained from the 'AquaResp' software (\href{www.aquaresp.com}{free software})
#' @param date.format  A string indicating date format used in the input (e.g. "%Y-%m-%d")
#' 
#' @return A dataframe containing the summary data on the recorded phases
#' 
#' @export
#' 
load.aquaresp.file <- function(file, date.format) {
  V1 <- NULL
  # Loading AquaResp data
  aquaresp.data <- utils::read.table(file, sep = ";", skip = 16, header = FALSE, strip.white = TRUE)
  aquaresp.data <- subset(aquaresp.data, select=V1) # recommend replacing subset with functions that do not require calling column names as objects.
  colnames(aquaresp.data) <- "Date.Time"
  aquaresp.data$Date.Time <- as.POSIXct(as.character(aquaresp.data$Date.Time), format = paste(date.format, "%H:%M:%S"))
  return(aquaresp.data)
}

