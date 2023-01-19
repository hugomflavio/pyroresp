#' Load a respirometry data file originating from pyroscience's 'Pyro Oxygen Logger' software.
#' 
#' @param file The name of a file which contains raw data obtained from the 'Pyro Oxygen Logger' software (\href{https://www.pyro-science.com}{PyroScience}) 
#' @param n.chamber	integer: the number of chambers used in an experiment (including empty ones)
#' @param date.format	string: date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' 
#' @return A dataframe containing the recorded oxygen data
#' 
#' @export
#' 
load.pyroscience.logger.file <- function(file, n.chamber = 1:4, date.format) {
	# NOTE: The number of chambers could likely be obtained from the file contents.
	#			 Would need an example to verify.
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
#' @param date.format	A string indicating date format used in raw data obtained from the 'Pyro Oxygen Logger' software (e.g. "\%Y-\%m-\%d")
#' 
#' @return A dataframe containing the recorded oxygen data
#' 
#' @export
#' 
load.pyroscience.workbench.file <- function(file, date.format) {
	if (!file.exists(file))
		stop("Could not find target file.", call. = FALSE)

	aux <- readLines(file, n = 100)
	preamble <- suppressWarnings(max(grep("^#", aux)))
	
	aux <- suppressWarnings(aux[grepl("Ch.[1-4]\\] - Oxygen Sensor", aux)])
	n.chambers <- length(sub("\\D*(\\d+).*", "\\1", aux))

	rm(aux)

	column.vector <- c(1, 2, 12, 17, 4, 30, 35, 22, 48, 53, 40, 66, 71, 58)[1:(2 + 3 * n.chambers)]
	column.names <- c("Date", "Time", "Temp.1", "Pressure.1", "Ox.1", "Temp.2", "Pressure.2", "Ox.2", "Temp.3", "Pressure.3", "Ox.3", "Temp.4", "Pressure.4", "Ox.4")[1:(2 + 3 * n.chambers)]

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

	if (any(is.na(pyro[,-1:-2])))
		warning('NA values found in the data!', immediate. = TRUE, call. = FALSE)

	return(pyro)
}

#' Dummy documentation
#' 
#' @export
#' 
patch.NAs <- function(input, method = c('linear', 'before', 'after'), verbose = TRUE) {
	columns.to.check <- colnames(input)[-c(1, grep("Phase", colnames(input)))]
	logical_input <- apply(input, 2, function(x) is.na(x))
	rle_list <- apply(logical_input, 2, rle)

	capture <- lapply(columns.to.check, function(i) {
		# cat(i, '\n')
		aux <- cumsum(rle_list[[i]]$lengths)
		breaks <- data.frame(Value = rle_list[[i]]$values,
												 Start = c(1, aux[-length(aux)] + 1),
												 Stop = aux)

		if (any(breaks$Value)) {
			nas <- breaks[breaks$Value, ]
			for (j in 1:nrow(nas)) {
				# cat(j, '\n')
				if (method == 'linear') {
					if (nas$Start[j] == 1) {
						if (verbose) warning("NAs found at the start of a column. Using method = 'after' for this instance.", immediate. = TRUE, call. = FALSE)
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Stop[j] + 1, i]
					} 
					else if (nas$Stop[j] == nrow(input)) {
						if (verbose) warning("NAs found at the end of a column. Using method = 'before' for this instance.", immediate. = TRUE, call. = FALSE)
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Start[j] - 1, i]
					} 
					else {
						fake_values <- seq(from = input[nas$Start[j] - 1, i],
															 to = input[nas$Stop[j] + 1, i],
															 length.out = nas$Stop[j] - nas$Start[j] + 3) # 1 for the before value, 1 for the after, and 1 for the value that is eliminated by the subtraction
						input[(nas$Start[j] - 1):(nas$Stop[j] + 1), i] <<- fake_values					
					}
				}
				if (method == 'before') {
					if (nas$Start[j] == 1) {
						if (verbose) warning("NAs found at the start of a column. Using method = 'after' for this instance.", immediate. = TRUE, call. = FALSE)
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Stop[j] + 1, i]
					} 
					else {
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Start[j] - 1, i]
					}
				}
				if (method == 'after') {
					if (nas$Stop[j] == nrow(input)) {
						if (verbose) warning("NAs found at the end of a column. Using method = 'before' for this instance.", immediate. = TRUE, call. = FALSE)
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Start[j] - 1, i]
					}
					else {
						input[nas$Start[j]:nas$Stop[j], i] <<- input[nas$Stop[j] + 1, i]						
					}				 
				}
			}
		}
	})

	return(input)
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
		stop("The PyroScience dataframe columns do not match the expected.\n			 The column names should be: '", paste(column.names[1:ncol(x)], collapse = "', '"), "'.", call. = FALSE)

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

	# grab first row to compile col names
	aux <- utils::read.table(file, sep = "\t", skip = 24, header = FALSE, strip.white = TRUE, nrows = 1)
	# grab rest of rows with the actual data
	output <- utils::read.table(file, sep = "\t", skip = 25, header = FALSE, strip.white = TRUE)

	a <- gsub(" .*$", "", aux[1,]) # grab first word only
	b <- gsub("^.* ", "", aux[1,]) # grab last word only
	b <- gsub(']', '', b) # remove lingering bracket

	colnames(output) <- paste0(a, '.', b) # combine first and last words

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
load.pyroscience.pH.file <- function(file, date.format, tz = Sys.timezone()) {
	if (length(file) == 0 || !file.exists(file))
		stop("Could not find target file.", call. = FALSE)

	# grab first row to compile col names
	aux <- utils::read.table(file, sep = "\t", skip = 20, header = FALSE, strip.white = TRUE, nrows = 1)
	# grab rest of rows with the actual data
	output <- utils::read.table(file, sep = "\t", skip = 21, header = FALSE, strip.white = TRUE)

	a <- gsub(" .*$", "", aux[1,]) # grab first word only
	b <- gsub("^.* ", "", aux[1,]) # grab last word only
	b <- gsub(']', '', b) # remove lingering bracket

	colnames(output) <- paste0(a, '.', b) # combine first and last words

	output <- output[, !duplicated(colnames(output))]

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
load_pyro_files <- function(folder, coolterm.legacy = FALSE) {
	if (!dir.exists(folder))
		stop('Could not find target folder')

	coolterm_file <- paste0(folder, "/", list.files(folder)[grepl("CoolTerm", list.files(folder))])
	
	if (length(coolterm_file) == 0)
		stop('Could not find coolterm file')

	output <- list()

	if (coolterm.legacy)
		output$coolterm <- load.coolterm.file(coolterm_file)
	else
		output$coolterm <- load_wide_coolterm_file(coolterm_file)

	output$pyro <- compile_run(folder)
	return(output)
}



#' dummy documentation
#' 
#' scan a target folder for pyroscience files matching a pattern and combine them.
#' 
#' @param folder the pyroscience run folder, containing a "ChannelData" folder inside
#' @param pattern The type of variable to look for. Currently accepted values: "Oxygen", "pH"
#' 
#' @export
#' 
compile_run <- function(folder) {
	files <- list.files(paste0(folder, '/ChannelData/'))

	file.link <- grepl("Oxygen|pH", files)

	if (all(!file.link)) {
		stop('No probe files found')
	}

	files <- files[file.link]

	aux <- lapply(files, function(i) {
		ch <- stringr::str_extract(i,'(?<=Ch.)[0-9]')

		if (grepl("Oxygen", i)) {
			x <- load.pyroscience.o2.file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
			x <- x[, c('Date.Time', 'Sample.CompT', 'Pressure.CompP', 'Oxygen.Main')]
			colnames(x)[2:4] <- paste0(c('Temp.', 'Pressure.', 'Ox.'),	ch)
		}

	 
		if (grepl("pH", i)) {
			x <- load.pyroscience.pH.file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
			x <- x[, c('Date.Time', 'pH.Main')]
			colnames(x)[2] <- paste0(c('pH.'),	ch)
		}

		return(x)
	})


	very.start <- min(as.POSIXct(sapply(aux, function(i) {
		as.character(min(i$Date.Time))
	})))

	very.end <- max(as.POSIXct(sapply(aux, function(i) {
		as.character(max(i$Date.Time))
	})))

	recipient <- data.frame(Date.Time = seq(from = very.start, to = very.end, by = 1),
													Phase = NA_character_)

	for (i in aux) {
		recipient <- merge(recipient, i[!duplicated(i$Date.Time), ], by = 'Date.Time', all = TRUE)
	}

	head(recipient)

	return(recipient)
}



#' dummy documentation
#' 
#' scan a target folder for a coolterm file and a pyroscience experiment file. import them.
#' 
#' @export
#' 
load_pyro_raw_files <- function(folder) {
	if (!dir.exists(folder))
		stop('Could not find target folder')

	aux <- paste0(folder, "/", list.files(folder)[grepl(".txt", list.files(folder))])
	
	coolterm_file <- aux[grepl("CoolTerm", aux)]
	
	if (length(coolterm_file) == 0)
		stop('Could not find coolterm file')


	output <- list()

	output$coolterm <- load.coolterm.file(coolterm_file)
	output$pyro <- compile_run(folder)
	return(output)
}



#' dummy documentation
#' 
#' perform standard processing operations to the pyro/coolterm files
#' 
#' @export
#' 
process_pyro_files <- function(input, wait, chamber.info, O2_unit, min_temp, max_temp) {

	all_units = c("percent_a.s.", "percent_o2", "hPa", "kPa",																											
				 "torr", "mmHg", "inHg", "mg_per_l", "ug_per_l", "umol_per_l",																							
				 "mmol_per_l", "ml_per_l", "mg_per_kg", "ug_per_kg", "umol_per_kg",																				 
				 "mmol_per_kg", "ml_per_kg", "volumes_percent")																														 

	if (!(O2_unit %in% all_units))																																										
			stop("the 'O2_argument' argument is not an acceptable unit. Please choose one of the following: ", paste(all_units, collapse = ", ")) 

  message("Merging pyroscience and coolterm file")
	input$phased <- merge_pyroscience_wide_coolterm(input$pyro, input$coolterm)
  message("Patching NA's in the data")
	input$phased <- patch.NAs(input$phased, method = 'linear', verbose = FALSE)
	message("Melting resp data into computer-friendly format")
  	input$meas_raw <- melt_resp(input = input$phased, info.data = chamber.info, O2.unit = O2_unit)
	message("Removing flush and wait values")
  	input$meas <- clean.meas(input = input$meas_raw, wait = 30)

	message("Converting O2 units and measuring delta")
  	input$meas$O2.umol.l <- respirometry::conv_o2(input$meas$O2.hPa, from = O2_unit, to = 'umol_per_l', temp = input$meas$Temp, sal = 0, atm_pres = input$meas$Pressure)
	input$meas$O2.a.s <- respirometry::conv_o2(input$meas$O2.hPa, from = O2_unit, to = 'percent_a.s.', temp = input$meas$Temp, sal = 0, atm_pres = input$meas$Pressure)

	input$meas <- calc_delta(input$meas, O2_col = "O2.umol.l")

	if (!missing(min_temp)) {
		time_break <- input$meas$Date.Time[which(input$meas$Phase == input$meas$Phase[head(which(input$meas$Temp > min_temp), 1)])[1]]
		input$meas <- input$meas[input$meas$Date.Time >= time_break, ]
	}

	if (!missing(max_temp)) {
		time_break <- input$meas$Date.Time[which(input$meas$Phase == input$meas$Phase[head(which(input$meas$Temp < max_temp), 1)])[1]]
		input$meas <- input$meas[input$meas$Date.Time >= time_break, ]		
	}

	return(input)
}

#' dummy documentation
#' 
#' perform standard metabolic rate calculations in pyro datasets
#' 
#' @export
#' 
process_pyro_mr <- function(input, r2, O2_raw, smr.method = "calcSMR.low10pc", max.length = 99999) {
	input$all.slopes <- calc.slope(input$corrected, O2_raw = O2_raw, max.length = max.length)
	input$good.slopes <- extract.slope(input$all.slopes, r2 = r2)
	input$smr.slope <- extract.slope(input$good.slopes, method = smr.method, r2 = r2)
	input$mmr.slope <- extract.slope(input$good.slopes, method = "max", r2 = r2, n.slope = 1)

	input$mr <- calculate.MR(input$good.slopes)
	input$smr <- calculate.MR(input$smr.slope)
	input$mmr <- calculate.MR(input$mmr.slope)

	# convert O2/Kg to O2/g
	input$smr$MR.mass.umol.g <- input$smr$MR.mass/1000
	input$mmr$MR.mass.umol.g <- input$mmr$MR.mass/1000
	input$mr$MR.mass.umol.g <- input$mr$MR.mass/1000

	return(input)
}




# load.pyroscience.workbench.file <- function(file, date.format, o2_from, o2_to) {
#	 if (!missing(o2_from) & !missing(o2_to))
#		 convert_o2 <- TRUE
#	 else
#		 convert_o2 <- FALSE

#	 if (!file.exists(file))
#		 stop("Could not find target file.", call. = FALSE)

#	 aux <- readLines(file, n = 100)
#	 preamble <- suppressWarnings(max(grep("^#", aux)))
	
#	 aux <- suppressWarnings(aux[grepl("Ch.[1-4]\\] - Oxygen Sensor", aux)])
#	 n.chambers <- length(sub("\\D*(\\d+).*", "\\1", aux))

#	 rm(aux)

#	 if (convert_o2) {
#		 column.vector <- c(1, 2, 12, 17, 4, 30, 35, 22, 48, 53, 40, 66, 71, 58)[1:(2 + 3 * n.chambers)]
#		 column.names <- c("Date", "Time", "Temp.1", "Pressure.1", "Ox.1", "Temp.2", "Pressure.2", "Ox.2", "Temp.3", "Pressure.3", "Ox.3", "Temp.4", "Pressure.4", "Ox.4")[1:(2 + 3 * n.chambers)]
#	 } else {
#		 column.vector <- c(1, 2, 12, 4, 30, 22, 48, 40, 66, 58)[1:(2 + 2 * n.chambers)]
#		 column.names <- c("Date", "Time", "Temp.1", "Ox.1", "Temp.2", "Ox.2", "Temp.3", "Ox.3", "Temp.4", "Ox.4")[1:(2 + 2 * n.chambers)]
#	 }

#	 pyro <- as.data.frame(data.table::fread(file, sep = "\t", skip = preamble, strip.white = TRUE, tz = ""))

#	 pyro <- pyro[, column.vector]
#	 names(pyro) <- column.names

#	 pyro$Date.Time <- paste(pyro$Date, pyro$Time)
#	 pyro$Phase <- NA
#	 pyro[pyro == "---"] <- NA
	
#	 pyro <- pyro[, c(ncol(pyro) - 1, ncol(pyro), 3:(ncol(pyro)-2))]
	
#	 pyro$Date.Time <- as.POSIXct(pyro$Date.Time, format = paste(date.format, "%H:%M:%S"), tz = Sys.timezone())

#	 # failsafe in case the file comes with a trailing NA line
#	 if (all(is.na(pyro[nrow(pyro), ])))
#		 pyro <- pyro[-nrow(pyro), ]

#	 if (convert_o2) {
#		 for (i in 1:n.chambers) {
#			 pyro[, paste0("Ox.", i)] <- 
#				 respirometry::conv_o2(o2 = pyro[, paste0("Ox.", i)],
#															 from = o2_from,
#															 to = o2_to,
#															 temp = pyro[, paste0("Temp.", i)],
#															 sal = 0, atm_pres = pyro[, paste0("Pressure.", i)])
#		 }
#		 pyro <- pyro[,!grepl("Pressure", colnames(pyro))]
#	 }

#	 if (any(is.na(pyro[,-1:-2])))
#		 warning('NA values found in the data!', immediate. = TRUE, call. = FALSE)

#	 return(pyro)
# }


#' dummy documentation
#' 
#' \
#' 
#' @export
#' 
remove_phase <- function(input, chamber, phase) {
	if (missing(chamber))
		input[!(input$Phase %in% phase), ]
	else
		input[!(input$Chamber.No %in% chamber & input$Phase %in% phase), ]
}
