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
load.pyroscience.o2.file <- function(file, date.format, skip = 0, tz = Sys.timezone()) {
	if (length(file) == 0 || !file.exists(file))
		stop("Could not find target file.", call. = FALSE)

	# grab first row to compile col names
	aux <- utils::read.table(file, sep = "\t", skip = 24, header = FALSE, strip.white = TRUE, nrows = 1)
	# grab rest of rows with the actual data
	output <- utils::read.table(file, sep = "\t", skip = 25 + skip, header = FALSE, strip.white = TRUE)

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
load.pyroscience.pH.file <- function(file, date.format, skip = 0, tz = Sys.timezone()) {
	if (length(file) == 0 || !file.exists(file))
		stop("Could not find target file.", call. = FALSE)

	# grab first row to compile col names
	aux <- utils::read.table(file, sep = "\t", skip = 20, header = FALSE, strip.white = TRUE, nrows = 1)
	# grab rest of rows with the actual data
	output <- utils::read.table(file, sep = "\t", skip = 21 + skip, header = FALSE, strip.white = TRUE)

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
load.pyroscience.temp.file <- function(file, date.format, skip = 0, tz = Sys.timezone()) {
	if (length(file) == 0 || !file.exists(file))
		stop("Could not find target file.", call. = FALSE)

	output <- utils::read.table(file, sep = "\t", skip = 12 + skip, header = FALSE, strip.white = TRUE)

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
load_pyro_files <- function(folder, legacy.coolterm = FALSE, legacy.device.name = TRUE, max.gap.fix = 1) {
	if (length(folder) == 0 || length(folder) > 1 || !dir.exists(folder))
		stop('Could not find target folder')

	coolterm_file <- list.files(folder)[grepl("CoolTerm", list.files(folder))]
	
	if (length(coolterm_file) == 0)
		stop('Could not find coolterm file')
	else
		coolterm_file <- paste0(folder, "/", coolterm_file)

	output <- list()

	coolterm <- lapply(coolterm_file, function(file) {
		if (legacy.coolterm)
			load_legacy_coolterm_file(file, max.gap.fix = max.gap.fix)
		else
			load_wide_coolterm_file(file, max.gap.fix = max.gap.fix)
	})

	if (legacy.device.name) {
		if(length(coolterm) > 1) {
			names(coolterm) <- stringr::str_extract(coolterm_file,'[A-Z](?=.txt)')
		} else {
			names(coolterm) <- "A"
		}
	} else {
		names(coolterm) <- stringr::str_extract(coolterm_file,'(?<=_)[^_]*(?=.txt)')
	}
	
	output$coolterm <- coolterm
	output$pyro <- compile_run(folder, legacy.device.name = legacy.device.name)
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
compile_run <- function(folder, legacy.device.name = TRUE) {
	files <- list.files(paste0(folder, '/ChannelData/'))

	file.link <- grepl("Oxygen|pH", files)

	if (all(!file.link)) {
		stop('No probe files found')
	}

	files <- files[file.link]

	sourcedata <- lapply(files, function(i) {
		if (legacy.device.name) {
			device <- stringr::str_extract(i,'[A-Z](?= Ch.)')
		} else {
			aux <- readLines(paste0(folder, '/ChannelData/', i), n = 10)
			aux <- aux[grepl("^#Device", aux)]
			device <- stringr::str_extract(aux,'(?<= )[^ ]*')
		}
		ch <- stringr::str_extract(i,'(?<=Ch.)[0-9]')

		if (grepl("Oxygen", i)) {
			x <- load.pyroscience.o2.file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
			x <- x[, c('Date.Time', 'Sample.CompT', 'Pressure.CompP', 'Oxygen.Main')]
			colnames(x)[2:4] <- paste0(c('Temp.', 'Pressure.', 'Ox.'),	device, ch)
			file.type <- "Oxygen"
		}

	 
		if (grepl("pH", i)) {
			x <- load.pyroscience.pH.file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
			x <- x[, c('Date.Time', 'pH.Main')]
			colnames(x)[2] <- paste0(c('pH.'),	device, ch)
			file.type <- "pH"
		}

		attributes(x)$sourcefile <- paste0(folder, '/ChannelData/', i)
		attributes(x)$device <- device
		attributes(x)$ch <- ch
		attributes(x)$filetype <- file.type
		return(x)
	})


	very.start <- min(as.POSIXct(sapply(sourcedata, function(i) {
		as.character(min(i$Date.Time))
	})))

	very.end <- max(as.POSIXct(sapply(sourcedata, function(i) {
		as.character(max(i$Date.Time))
	})))

	recipient <- data.frame(Date.Time = seq(from = very.start, to = very.end, by = 1),
													Phase = NA_character_)

	for (i in sourcedata) {
		recipient <- merge(recipient, i[!duplicated(i$Date.Time), ], by = 'Date.Time', all = TRUE)
	}

	attributes(recipient)$latestbatchstart <- 1

	output <- list(sourcedata = sourcedata, compileddata = recipient)
	return(output)
}



update_run <- function(input) {

	new_sourcedata <- lapply(input$sourcedata, function(i) {
		device <- attributes(i)$device
		ch <- attributes(i)$ch

		if (attributes(i)$filetype == "Oxygen") {
			x <- load.pyroscience.o2.file(attributes(i)$sourcefile, date.format = '%d-%m-%Y', skip = nrow(i))
			x <- x[, c('Date.Time', 'Sample.CompT', 'Pressure.CompP', 'Oxygen.Main')]
			colnames(x)[2:4] <- paste0(c('Temp.', 'Pressure.', 'Ox.'),	device, ch)
			# x <- rbind(i, x)
		}

		if (attributes(i)$filetype == "pH") {
			x <- load.pyroscience.pH.file(attributes(i)$sourcefile, date.format = '%d-%m-%Y', skip = nrow(i))
			x <- x[, c('Date.Time', 'pH.Main')]
			colnames(x)[2] <- paste0(c('pH.'),	device, ch)
			# x <- rbind(i, x)
		}

		attributes(x)$sourcefile <- attributes(i)$sourcefile
		attributes(x)$device <- device
		attributes(x)$ch <- ch
		attributes(x)$filetype <- attributes(i)$filetype
		return(x)
	})

	sourcedata <- lapply(1:length(new_sourcedata), function(i) {
		output <- rbind(input$sourcedata[[i]], new_sourcedata[[i]])
		attributes(output)$sourcefile <- attributes(input$sourcedata[[i]])$sourcefile
		attributes(output)$device <- attributes(input$sourcedata[[i]])$device
		attributes(output)$ch <- attributes(input$sourcedata[[i]])$ch
		attributes(output)$filetype <- attributes(input$sourcedata[[i]])$filetype
		return(output)
	})


	very.start <- min(as.POSIXct(sapply(new_sourcedata, function(i) {
		as.character(min(i$Date.Time))
	})))

	very.end <- max(as.POSIXct(sapply(new_sourcedata, function(i) {
		as.character(max(i$Date.Time))
	})))

	recipient <- data.frame(Date.Time = seq(from = very.start, to = very.end, by = 1),
													Phase = NA_character_)

	for (i in new_sourcedata) {
		recipient <- merge(recipient, i[!duplicated(i$Date.Time), ], by = 'Date.Time', all = TRUE)
	}

	recipient <- rbind(input$compileddata, recipient)
	attributes(recipient)$latestbatchstart <- nrow(input$compileddata)+1

	output <- list(sourcedata = sourcedata, compileddata = recipient)
	return(output)
}


#' dummy documentation
#' 
#' perform standard processing operations to the pyro/coolterm files
#' 
#' @export
#' 
process_pyro_files <- function(input, wait, chamber.info, O2_unit, min_temp, max_temp, start_time, stop_time, patch_NAs = TRUE) {

	all_units = c("percent_a.s.", "percent_o2", "hPa", "kPa",																											
				 "torr", "mmHg", "inHg", "mg_per_l", "ug_per_l", "umol_per_l",																							
				 "mmol_per_l", "ml_per_l", "mg_per_kg", "ug_per_kg", "umol_per_kg",																				 
				 "mmol_per_kg", "ml_per_kg", "volumes_percent")																														 

	if (!(O2_unit %in% all_units))																																										
			stop("the 'O2_argument' argument is not an acceptable unit. Please choose one of the following: ", paste(all_units, collapse = ", ")) 

	if (!missing("chamber.info")) {
		input$chamber.info <- chamber.info
	}

  	message("Merging pyroscience and coolterm file")
	input$phased <- merge_pyroscience_wide_coolterm(input)
  	if (patch_NAs) {
  		message("Patching NA's in the data")
		input$phased <- patch.NAs(input$phased, method = 'linear', verbose = FALSE)
	}
	message("Melting resp data into computer-friendly format")
  	input$meas_raw <- melt_resp(input = input$phased, info.data = chamber.info, O2.unit = O2_unit)
	message("Removing flush and wait values")
  	input$meas <- clean.meas(input = input$meas_raw, wait = wait)

	message("Converting O2 units and measuring delta")
	input$meas_raw$O2.umol.l <- respirometry::conv_o2(input$meas_raw[, paste0("O2.", O2_unit)], from = O2_unit, to = 'umol_per_l', temp = input$meas_raw$Temp, sal = 0, atm_pres = input$meas_raw$Pressure)
	input$meas$O2.umol.l <- respirometry::conv_o2(input$meas[, paste0("O2.", O2_unit)], from = O2_unit, to = 'umol_per_l', temp = input$meas$Temp, sal = 0, atm_pres = input$meas$Pressure)
	input$meas$O2.a.s <- respirometry::conv_o2(input$meas[, paste0("O2.", O2_unit)], from = O2_unit, to = 'percent_a.s.', temp = input$meas$Temp, sal = 0, atm_pres = input$meas$Pressure)
 
	input$meas <- calc_delta(input$meas, O2_col = "O2.umol.l")

	if (!missing(min_temp)) {
		cutoff <- head(which(input$meas$Temp > min_temp), 1)
		if (length(cutoff) == 0) {
			stop ("Temperature never rose above ", min_temp, ".")
		} else {
			first_true <- which(input$meas$Phase == input$meas$Phase[cutoff])[1]
			time_break <- input$meas$Date.Time[first_true]
			input$meas <- input$meas[input$meas$Date.Time >= time_break, ]
		}
	}

	if (!missing(max_temp)) {
		cutoff <- head(which(input$meas$Temp < max_temp), 1)
		if (length(cutoff) == 0) {
			stop ("Temperature never dropped below ", max_temp, ".")
		} else {
			first_true <- which(input$meas$Phase == input$meas$Phase[cutoff])[1]
			time_break <- input$meas$Date.Time[first_true]
			input$meas <- input$meas[input$meas$Date.Time >= time_break, ]
		}	
	}

	if (!missing(start_time)) {
		cutoff <- head(which(input$meas$Date.Time >= as.POSIXct(start_time)), 1)
		if (length(cutoff) == 0) {
			stop ("Data ends before ", start_time, ".")
		} else {
			first_phase <- input$meas$Phase[cutoff]
			first_true <- head(which(input$meas$Phase == first_phase), 1)
			time_break <- input$meas$Date.Time[first_true]
			input$meas <- input$meas[input$meas$Date.Time >= time_break, ]
		}
	}

	if (!missing(stop_time)) {
		cutoff <- tail(which(input$meas$Date.Time <= as.POSIXct(stop_time)), 1)
		if (length(cutoff) == 0) {
			stop ("Data starts after ", stop_time, ".")
		} else {
			last_phase <- input$meas$Phase[cutoff]
			last_true <- tail(which(input$meas$Phase == last_phase), 1)
			time_break <- input$meas$Date.Time[last_true]
			input$meas <- input$meas[input$meas$Date.Time <= time_break, ]
		}
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
	
	aux <- extract.slope(input$good.slopes, method = smr.method, r2 = r2)
	input$smr.slope <- merge(input$chamber.info, aux[, !(colnames(aux) %in% c("ID", "Mass", "Volume"))], by = "Probe", all = TRUE)
	
	aux <- extract.slope(input$good.slopes, method = "max", r2 = r2, n.slope = 1)
	input$mmr.slope <- merge(input$chamber.info, aux[, !(colnames(aux) %in% c("ID", "Mass", "Volume"))], by = "Probe", all = TRUE)
	
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
#' @export
#' 
remove_phase <- function(input, chamber, phase) {
	if (missing(chamber))
		input[!(input$Phase %in% phase), ]
	else
		input[!(input$Probe %in% chamber & input$Phase %in% phase), ]
}
