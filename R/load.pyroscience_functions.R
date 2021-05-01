#' Load a respirometry data file originating from pyroscience's 'Pyro Oxygen Logger' software.
#' 
#' @param file The name of a file which contains raw data obtained from the 'Pyro Oxygen Logger' software (\href{https://www.pyro-science.com}{PyroScience}) 
#' @param n.chamber  integer: the number of chambers used in an experiment (including empty ones)
#' @param date.format  string: date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' 
#' @return A dataframe containing the recorded oxygen data
#' 
#' @export
#' 
load.pyroscience.logger.file <- function(file, n.chamber = 1:4, date.format) {
  # NOTE: The number of chambers could likely be obtained from the file contents.
  #       Would need an example to verify.
  n.chamber <- match.arg(n.chamber)

  if (!file.exists(file))
    stop("Could not find target file.", call. = FALSE)

  column.vector <- c(1, 2, 9, 5, 10, 6, 11, 7, 12, 8)[1:(2 + 2 * n.chamber)]
  column.names <- c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")[1:(2 + 2 * n.chamber)]

  pyro.data <- utils::read.table(file, sep = "\t", skip = 12 + (2 * n.chamber), header = FALSE, strip.white = TRUE)
  pyro.data <- pyro.data[, column.vector]
  names(pyro.data) <- column.names
  pyro.data$Date.Time <- paste(pyro.data$Date, pyro.data$Time)
  pyro.data$Phase <- "F"
  pyro.data[pyro.data == "---"] <- NA
  
  pyro.data <- pyro.data[, c(ncol(pyro.data) - 1, ncol(pyro.data), 3:(ncol(pyro.data)-2))]
  
  pyro.data$Date.Time <- as.POSIXct(pyro.data$Date.Time, format = paste(date.format, "%H:%M:%S"))
  return(pyro.data)
}

#' Load a respirometry data file originating from pyroscience's 'PyroScience Workbench' software.
#' 
#' @param file The name of a file which contains raw data obtained from the 'PyroScience Workbench' software (\href{https://www.pyro-science.com}{PyroScience}) 
#' @param date.format  A string indicating date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' @param o2_from Optional: a string describing the units in which o2 was measured. All options are available at \code{\link[respirometry]{conv_o2}}
#' @param o2_to Optional: a string describing the unit to which the o2 measurements should be converted. All options are available at \code{\link[respirometry]{conv_o2}}
#' 
#' @return A dataframe containing the recorded oxygen data
#' 
#' @export
#' 
load.pyroscience.workbench.file <- function(file, date.format, o2_from, o2_to) {
  if (!missing(o2_from) & !missing(o2_to))
    convert_o2 <- TRUE
  else
    convert_o2 <- FALSE

  if (!file.exists(file))
    stop("Could not find target file.", call. = FALSE)

  aux <- readLines(file, n = 100)
  preamble <- sum(grepl("^#", aux))
  last.chamber <- tail(aux[grepl("Ch.[1-4]\\] - Oxygen Sensor", aux)], 1)
  n.chambers <- as.numeric(sub("\\D*(\\d+).*", "\\1", last.chamber))

  rm(aux)

  if (convert_o2) {
    column.vector <- c(1, 2, 12, 17, 4, 30, 35, 22, 48, 53, 40, 66, 71, 58)[1:(2 + 3 * n.chambers)]
    column.names <- c("Date", "Time", "Temp.1", "Pressure.1", "Ox.1", "Temp.2", "Pressure.2", "Ox.2", "Temp.3", "Pressure.3", "Ox.3", "Temp.4", "Pressure.4", "Ox.4")[1:(2 + 3 * n.chambers)]
  } else {
    column.vector <- c(1, 2, 12, 4, 30, 22, 48, 40, 66, 58)[1:(2 + 2 * n.chambers)]
    column.names <- c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")[1:(2 + 2 * n.chambers)]
  }

  pyro.data <- as.data.frame(data.table::fread(file, sep = "\t", skip = preamble, strip.white = TRUE, tz = ""))

  pyro.data <- pyro.data[, column.vector]
  names(pyro.data) <- column.names

  pyro.data$Date.Time <- paste(pyro.data$Date, pyro.data$Time)
  pyro.data$Phase <- NA
  pyro.data[pyro.data == "---"] <- NA
  
  pyro.data <- pyro.data[, c(ncol(pyro.data) - 1, ncol(pyro.data), 3:(ncol(pyro.data)-2))]
  
  pyro.data$Date.Time <- as.POSIXct(pyro.data$Date.Time, format = paste(date.format, "%H:%M:%S"), tz = Sys.timezone())

  # failsafe in case the file comes with a trailing NA line
  if (all(is.na(pyro.data[nrow(pyro.data), ])))
    pyro.data <- pyro.data[-nrow(pyro.data), ]

  if (convert_o2) {
    for (i in 1:n.chambers) {
      pyro.data[, paste0("Ox.", i)] <- 
        respirometry::conv_o2(o2 = pyro.data[, paste0("Ox.", i)],
                              from = o2_from,
                              to = o2_to,
                              temp = pyro.data[, paste0("Temp.", i)],
                              sal = 0, atm_pres = pyro.data[, paste0("Pressure.", i)])
    }
    pyro.data <- pyro.data[,!grepl("Pressure", colnames(pyro.data))]
  }

  return(pyro.data)
}

#' Confirm that the pyroscience input contains the expected columns
#' 
#' @param x the input data.frame
#' 
#' @return a data frame with ordered columns
#' 
#' @keywords internal
#' 
check.pyrocience.data <- function(x) {
  column.names <- c("Date.Time", "Phase", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")
  check <- match(colnames(x), column.names[1:ncol(x)])

  if (any(is.na(check)))
    stop("The PyroScience dataframe columns do not match the expected.\n       The column names should be: '", paste(column.names[1:ncol(x)], collapse = "', '"), "'.", call. = FALSE)

  if (!any(grepl("POSIXct", class(x$Date.Time))))
    stop("Please format the column 'Date.Time' as POSIXct in the PyroScience input")

  return(x[, column.names[1:ncol(x)]])
}

