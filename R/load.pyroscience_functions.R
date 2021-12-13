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

  if (length(file) == 0 || !file.exists(file))
    stop("Could not find target file.", call. = FALSE)

  column.vector <- c(1, 2, 9, 5, 10, 6, 11, 7, 12, 8)[1:(2 + 2 * n.chamber)]
  column.names <- c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")[1:(2 + 2 * n.chamber)]

  pyro <- utils::read.table(file, sep = "\t", skip = 12 + (2 * n.chamber), header = FALSE, strip.white = TRUE)
  pyro <- pyro[, column.vector]
  names(pyro) <- column.names
  pyro$Date.Time <- paste(pyro$Date, pyro$Time)
  pyro$Phase <- "F"
  pyro[pyro == "---"] <- NA
  
  pyro <- pyro[, c(ncol(pyro) - 1, ncol(pyro), 3:(ncol(pyro)-2))]
  
  pyro$Date.Time <- as.POSIXct(pyro$Date.Time, format = paste(date.format, "%H:%M:%S"))
  return(pyro)
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
  preamble <- suppressWarnings(max(grep("^#", aux)))
  
  aux <- suppressWarnings(aux[grepl("Ch.[1-4]\\] - Oxygen Sensor", aux)])
  n.chambers <- length(sub("\\D*(\\d+).*", "\\1", aux))

  rm(aux)

  if (convert_o2) {
    column.vector <- c(1, 2, 12, 17, 4, 30, 35, 22, 48, 53, 40, 66, 71, 58)[1:(2 + 3 * n.chambers)]
    column.names <- c("Date", "Time", "Temp.1", "Pressure.1", "Ox.1", "Temp.2", "Pressure.2", "Ox.2", "Temp.3", "Pressure.3", "Ox.3", "Temp.4", "Pressure.4", "Ox.4")[1:(2 + 3 * n.chambers)]
  } else {
    column.vector <- c(1, 2, 12, 4, 30, 22, 48, 40, 66, 58)[1:(2 + 2 * n.chambers)]
    column.names <- c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")[1:(2 + 2 * n.chambers)]
  }

  pyro <- as.data.frame(data.table::fread(file, sep = "\t", skip = preamble, strip.white = TRUE, tz = ""))

  pyro <- pyro[, column.vector]
  names(pyro) <- column.names

  pyro$Date.Time <- paste(pyro$Date, pyro$Time)
  pyro$Phase <- NA
  pyro[pyro == "---"] <- NA
  
  pyro <- pyro[, c(ncol(pyro) - 1, ncol(pyro), 3:(ncol(pyro)-2))]
  
  pyro$Date.Time <- as.POSIXct(pyro$Date.Time, format = paste(date.format, "%H:%M:%S"), tz = Sys.timezone())

  # failsafe in case the file comes with a trailing NA line
  if (all(is.na(pyro[nrow(pyro), ])))
    pyro <- pyro[-nrow(pyro), ]

  if (convert_o2) {
    for (i in 1:n.chambers) {
      pyro[, paste0("Ox.", i)] <- 
        respirometry::conv_o2(o2 = pyro[, paste0("Ox.", i)],
                              from = o2_from,
                              to = o2_to,
                              temp = pyro[, paste0("Temp.", i)],
                              sal = 0, atm_pres = pyro[, paste0("Pressure.", i)])
    }
    pyro <- pyro[,!grepl("Pressure", colnames(pyro))]
  }

  if (any(is.na(pyro[,-1:-2])))
    warning('NA values found in the data!', immediate. = TRUE, call. = FALSE)

  return(pyro)
}

#' Dummy documentation
#' 
#' @export
#' 
patch.NAs <- function(input, method = c('linear', 'before', 'after')) {
  db <- input[,-1:-2]
  logical_db <- apply(db, 2, function(x) is.na(x))
  rle_list <- apply(logical_db, 2, rle)

  capture <- lapply (1:ncol(db), function(coln) {
    aux <- cumsum(rle_list[[coln]]$lengths)
    breaks <- data.frame(Value = rle_list[[coln]]$values,
                         Start = c(1, aux[-length(aux)] + 1),
                         Stop = aux)

    if (any(breaks$Value)) {
      nas <- breaks[breaks$Value, ]
      for (i in 1:nrow(nas)) {
        if (method == 'linear') {
          if (nas$Start[i] == 1) {
            warning("NAs found at the start of a column. Using method = 'after' for this instance.", immediate. = TRUE, call. = FALSE)
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Stop[i] + 1, coln]
          } 
          else if (nas$Stop[i] == nrow(db)) {
            warning("NAs found at the end of a column. Using method = 'before' for this instance.", immediate. = TRUE, call. = FALSE)
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Start[i] - 1, coln]
          } 
          else {
            fake_values <- seq(from = db[nas$Start[i] - 1, coln],
                               to = db[nas$Stop[i] + 1, coln],
                               length.out = nas$Stop[i] - nas$Start[i] + 3) # 1 for the before value, 1 for the after, and 1 for the value that is eliminated by the subtraction
            db[nas$Start[i]:nas$Stop[i], coln] <<- fake_values          
          }
        }
        if (method == 'before') {
          if (nas$Start[i] == 1) {
            warning("NAs found at the start of a column. Using method = 'after' for this instance.", immediate. = TRUE, call. = FALSE)
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Stop[i] + 1, coln]
          } 
          else {
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Start[i] - 1, coln]
          }
        }
        if (method == 'after') {
          if (nas$Stop[i] == nrow(db)) {
            warning("NAs found at the end of a column. Using method = 'before' for this instance.", immediate. = TRUE, call. = FALSE)
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Start[i] - 1, coln]
          }
          else {
            db[nas$Start[i]:nas$Stop[i], coln] <<- db[nas$Stop[i] + 1, coln]            
          }         
        }
      }
    }
  })

  return(as.data.frame(cbind(input[,1:2], db)))
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




#' Load a single channel file
#' 
#' @param file the input file
#' 
#' @return the imported channel file
#' 
#' @export
#' 
load.pyroscience.o2.file <- function(file, date.format, tz = Sys.timezone()) {
  if (length(file) == 0 || !file.exists(file))
    stop("Could not find target file.", call. = FALSE)

  aux <- utils::read.table(file, sep = "\t", skip = 24, header = FALSE, strip.white = TRUE)

  a <- gsub(" .*$", "", aux[1,])

  b <- gsub("^.* ", "", aux[1,])
  b <- gsub(']', '', b)

  output <- utils::read.table(file, sep = "\t", skip = 25, header = FALSE, strip.white = TRUE)

  colnames(output) <- paste0(a, '.', b)

  output$Date.Time <- as.POSIXct(paste(output$Date.Main, output$Time.Main), format = paste(date.format, "%H:%M:%S"), tz = tz)

  return(output)
}

#' Load a single channel file
#' 
#' @param file the input file
#' 
#' @return the imported channel file
#' 
#' @export
#' 
load.pyroscience.temp.file <- function(file, date.format, tz = Sys.timezone()) {
  if (length(file) == 0 || !file.exists(file))
    stop("Could not find target file.", call. = FALSE)

  output <- utils::read.table(file, sep = "\t", skip = 12, header = FALSE, strip.white = TRUE)

  colnames(output) <- c('Date', 'Time', 'ds', 'Temp', 'Status')

  output$Date.Time <- as.POSIXct(paste(output$Date, output$Time), format = paste(date.format, "%H:%M:%S"), tz = tz)

  output <- output[, c('Date.Time', 'ds', 'Temp', 'Status')]

  return(output)
}




#' dummy documentation
#' 
#' scan a target folder for a coolterm file and a pyroscience experiment file. import them.
#' 
#' @export
#' 
load_pyro_files <- function(folder) {
  if (!dir.exists(folder))
    stop('Could not find target folder')

  aux <- paste0(folder, "/", list.files(folder)[grepl(".txt", list.files(folder))])
  
  coolterm_file <- aux[grepl("CoolTerm", aux)]
  
  if (length(coolterm_file) == 0)
    stop('Could not find coolterm file')

  pyro_file <- aux[!grepl("CoolTerm", aux)]
  
  if (length(pyro_file) == 0)
    stop('could not find pyro file')

  if (length(pyro_file) > 1)
    stop('could not identify pyro file. too many txt files in folder?')

  output <- list()

  output$coolterm <- load.coolterm.file(coolterm_file)
  output$pyro <- load.pyroscience.workbench.file(pyro_file, date.format = "%d-%m-%Y")
  return(output)
}

#' dummy documentation
#' 
#' perform standard processing operations to the pyro/coolterm files
#' 
#' @export
#' 
process_pyro_files <- function(input, wait, chamber.info) {
  input$pyro <- patch.NAs(input$pyro, method = "linear")
  input$phased <- merge_pyroscience_coolterm(input$pyro, input$coolterm)
  input$phased <- input$phased[!is.na(input$phased$Phase),]
  input$meas <- clean.meas(input = input$phased, wait = wait)
  input$meas <- melt_resp(input$meas, chamber.info) 
  return(input)
}

#' dummy documentation
#' 
#' perform standard metabolic rate calculations in pyro datasets
#' 
#' @export
#' 
process_pyro_mr <- function(input, r2, smr.method = "calcSMR.low10pc") {
  input$all.slopes <- calc.slope(input$corrected)
  input$good.slopes <- extract.slope(input$all.slopes, r2 = r2)
  input$smr.slope <- extract.slope(input$good.slopes, method = smr.method, r2 = r2)
  input$mmr.slope <- extract.slope(input$good.slopes, method = "max", r2 = r2, n.slope = 1)

  input$mr <- calculate.MR(input$good.slopes)
  input$smr <- calculate.MR(input$smr.slope)
  input$mmr <- calculate.MR(input$mmr.slope)
  return(input)
}


