#' Convert Respirometry Data from PyroScience and AquaResp Software to the FishResp Format
#'
#' The function is used to convert raw data from 'Pyro Oxygen Logger' (\href{https://www.pyroscience.com/}{PyroScience}) and a summary file from 'AquaResp' (\href{https://www.aquaresp.com}{free software}) to 'FishResp' format. This function should be applied before usage of the functions \code{\link{import.test}} and \code{\link{import.meas}}. The output is a file containing raw respirometry data in the 'FishResp' format (see Details in \code{\link{import.test}} to read more information about the 'FishResp' format)
#'
#' @usage
#' pyroscience.aquaresp(pyroscience.file,
#'                      aquaresp.file,
#'                      fishresp.file,
#'                      n.chamber = c(1,2,3,4),
#'                      date.format = c("DMY", "MDY", "YMD"),
#'                      wait.phase = NA, measure.phase = NA)
#'
#' @param pyroscience.file  the name of a file which contains raw data obtained from the 'Pyro Oxygen Logger' software (\href{https://www.pyroscience.com/}{PyroScience})
#' @param aquaresp.file  the name of a file which contains summary data obtained from the 'AquaResp' software (\href{https://www.aquaresp.com}{free software})
#' @param fishresp.file  the name of an exported file containing raw data in the 'FishResp' format
#' @param n.chamber  integer: the number of chambers used in an experiment (including empty ones)
#' @param date.format  string: date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' @param wait.phase  integer: duration of the wait phase (in seconds), see the 'AquaResp' summary file (row #5)
#' @param measure.phase  integer: duration of the measure phase (in seconds), see the 'AquaResp' summary file (row #6)
#'
#' @return The function returns a data frame containing raw data in the 'FishResp' format
#'
#' @examples
#' \dontrun{
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
#' }
#' @export
pyroscience.aquaresp <- function(input,
                                 pyro.type = c("logger", "workbench"),
                                 aquaresp,
                                 n.chamber = 1:4,
                                 date.format,
                                 wait.phase = NA, measure.phase = NA){

  # Loading PyroScience data
  V1 <- V2 <- V3 <- V4 <- V5 <- V6 <- V7 <- V8 <- V9 <- V10 <- V11 <- V12 <- NULL

  if (n.chamber == 1){
    pyroscience.data <- read.table(pyroscience.file, sep = "\t", skip=14, header=F, strip.white=T)
    pyroscience.data <- subset(pyroscience.data, select=c(V1, V2, V9, V5))
    colnames(pyroscience.data)<-c("Date", "Time", "Temp.1", "Ox.1")
    for(i in 3:ncol(pyroscience.data)){pyroscience.data[,i] = as.numeric(gsub(',','.', pyroscience.data[,i]))}
    pyroscience.data$Date.Time <- paste(pyroscience.data$Date, pyroscience.data$Time, sep="/")
    pyroscience.data$Phase <- "F"
    pyroscience.data[pyroscience.data == "---"] <- NA
    pyroscience.data <- pyroscience.data[,c(5,6,3,4)]
    }

  else if (n.chamber == 2){
    pyroscience.data <- read.table(pyroscience.file, sep = "\t", skip=16, header=F, strip.white=T)
    pyroscience.data <- subset(pyroscience.data, select=c(V1, V2, V9, V5, V10, V6))
    colnames(pyroscience.data)<-c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2")
    for(i in 3:ncol(pyroscience.data)){pyroscience.data[,i] = as.numeric(gsub(',','.', pyroscience.data[,i]))}
    pyroscience.data$Date.Time <- paste(pyroscience.data$Date, pyroscience.data$Time, sep="/")
    pyroscience.data$Phase <- "F"
    pyroscience.data[pyroscience.data == "---"] <- NA
    pyroscience.data <- pyroscience.data[,c(7,8,3,4,5,6)]
    }

  else if (n.chamber == 3){
    pyroscience.data <- read.table(pyroscience.file, sep = "\t", skip=18, header=F, strip.white=T)
    pyroscience.data <- subset(pyroscience.data, select=c(V1, V2, V9, V5, V10, V6, V11, V7))
    colnames(pyroscience.data)<-c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3")
    for(i in 3:ncol(pyroscience.data)){pyroscience.data[,i] = as.numeric(gsub(',','.', pyroscience.data[,i]))}
    pyroscience.data$Date.Time <- paste(pyroscience.data$Date, pyroscience.data$Time, sep="/")
    pyroscience.data$Phase <- "F"
    pyroscience.data[pyroscience.data == "---"] <- NA
    pyroscience.data <- pyroscience.data[,c(9,10,3,4,5,6,7,8)]
    }

  else if (n.chamber == 4){
    pyroscience.data <- read.table(pyroscience.file, sep = "\t", skip=20, header=F, strip.white=T)
    pyroscience.data <- subset(pyroscience.data, select=c(V1, V2, V9, V5, V10, V6, V11, V7, V12, V8))
    colnames(pyroscience.data)<-c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")
    for(i in 3:ncol(pyroscience.data)){pyroscience.data[,i] = as.numeric(gsub(',','.', pyroscience.data[,i]))}
    pyroscience.data$Date.Time <- paste(pyroscience.data$Date, pyroscience.data$Time, sep="/")
    pyroscience.data$Phase <- "F"
    pyroscience.data[pyroscience.data == "---"] <- NA
    pyroscience.data <- pyroscience.data[,c(11,12,3,4,5,6,7,8,9,10)]
    }

  else{
    print("Please, choose the number of chambers between 1 and 4")
    }


  if(any(date.format == "DMY")){
    pyroscience.data$Date.Time<-strptime(as.character(pyroscience.data$Date.Time), "%d/%m/%Y/%H:%M:%S")
    ts<-times(strftime(pyroscience.data$Date.Time, "%H:%M:%S"))
    ds <- format(pyroscience.data$Date.Time, format="%d/%m/%Y")
    pyroscience.data$Date.Time <- paste(ds,ts, sep = "/")
    }

  else if(any(date.format == "MDY")){
    pyroscience.data$Date.Time<-strptime(as.character(pyroscience.data$Date.Time), "%m/%d/%Y/%H:%M:%S")
    ts<-times(strftime(pyroscience.data$Date.Time, "%H:%M:%S"))
    ds <- format(pyroscience.data$Date.Time, format="%m/%d/%Y")
    pyroscience.data$Date.Time <- paste(ds,ts, sep = "/")
    }

  else if(any(date.format == "YMD")){
    pyroscience.data$Date.Time<-strptime(as.character(pyroscience.data$Date.Time), "%Y/%m/%d/%H:%M:%S")
    ts<-times(strftime(pyroscience.data$Date.Time, "%H:%M:%S"))
    ds <- format(pyroscience.data$Date.Time, format="%Y/%m/%d")
    pyroscience.data$Date.Time <- paste(ds,ts, sep = "/")
    }


  # Loading AquaResp data
  aquaresp.data <- read.table(aquaresp.file, sep = ";", skip=16, header=F, strip.white=T)
  aquaresp.data <- subset(aquaresp.data, select=V1)
  colnames(aquaresp.data) <- "Date.Time"
  aquaresp.data$Date.Time<-strptime(as.character(aquaresp.data$Date.Time), "%Y-%m-%d %H:%M:%S")
  tms<-times(strftime(aquaresp.data$Date.Time, "%H:%M:%S"))

  if(any(date.format == "DMY")){
    dts <- format(aquaresp.data$Date.Time, format="%d/%m/%Y")
    x <- paste(dts,tms, sep = "/")
    }

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

