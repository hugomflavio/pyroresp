#' Fill in missing datapoints using the nearest available data.
#'
#' At high sampling rates, data points can occasionally be lost, causing odd
#' interruptions in the temperature and/or oxygen traces.
#' This function fills in those gaps.
#'
#' @param input A dataframe containing Timestamps on the first column and
#' 	matching data on the remaining. If input contains a "phase" column, it will
#' 	be transported to the output unchanged.
#' @param patch_method One of:
#' 	'linear' to capture the nearest before and after non-NA values and make a
#' 			 linear interpolation.
#' 	'before' to find the nearest non-NA value before the NA and use it to fill
#' 			 the gap.
#' 	'after'  to do the same as above but with the nearest value coming after
#' 			 the NA.
#' @param verbose Logical. Should exceptions be reported out loud.
#' 	Defaults to TRUE.
#'
#' @return The input table with the NAs filled in as requested.
#'
#' @export
#'
patch_NAs <- function(input, patch_method = c('linear', 'before', 'after'),
		verbose = TRUE) {

	patch_method <- match.arg(patch_method)

	columns_to_check <- colnames(input)[-c(1, grep("phase", colnames(input)))]
	logical_input <- apply(input, 2, is.na)
	rle_list <- apply(logical_input, 2, rle)

	# run the loop for each column separately
	capture <- lapply(columns_to_check, function(i) {
		# start by finding the gaps in the column and storing them in a table
		aux <- cumsum(rle_list[[i]]$lengths)
		breaks <- data.frame(Value = rle_list[[i]]$values,
							 Start = c(1, aux[-length(aux)] + 1),
							 Stop = aux)

		# if there are any breaks, start working on them
		if (any(breaks$Value)) {
			nas <- breaks[breaks$Value, ]
			# for every break found, apply the correction
			for (j in 1:nrow(nas)) {
				if (patch_method == 'linear') {
					# failsafe against NAs at the start when using linear
					if (nas$Start[j] == 1) {
						if (verbose) {
							warning("NAs found at the start of a column.",
								" Using method = 'after' for this instance.",
								immediate. = TRUE, call. = FALSE)
						}
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Stop[j] + 1, i]
					}
					# failsafe against NAs at the end when using linear
					if (nas$Stop[j] == nrow(input)) {
						if (verbose) {
							warning("NAs found at the end of a column.",
								" Using method = 'before' for this instance.",
								immediate. = TRUE, call. = FALSE)
						}
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Start[j] - 1, i]
					}
					if (nas$Start[j] != 1 && nas$Stop[j] != nrow(input)) {
						missing_interval <- (nas$Start[j] - 1):(nas$Stop[j] + 1)
						replacement <- seq(
							from = input[nas$Start[j] - 1, i],
							to = input[nas$Stop[j] + 1, i],
							length.out = nas$Stop[j] - nas$Start[j] + 3
						)
						# Explanation for the +3 in the length.out above:
						# +1 to repeat the last known value before the break
						# +1 to replace the first known value after the
						# break.
						# +1 to compensate for the one that is lost when
						# doing Stop - Start. E.g. row 3 - row 2 = 1, but both
						# row 3 and 2 are NAs
					}
				}
				if (patch_method == 'before') {
					if (nas$Start[j] == 1) {
						if (verbose) {
							warning("NAs found at the start of a column.",
								" Using method = 'after' for this instance.",
								immediate. = TRUE, call. = FALSE)
						}
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Stop[j] + 1, i]
					}
					else {
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Start[j] - 1, i]
					}
				}
				if (patch_method == 'after') {
					if (nas$Stop[j] == nrow(input)) {
						if (verbose) {
							warning("NAs found at the end of a column.",
								" Using method = 'before' for this instance.",
								immediate. = TRUE, call. = FALSE)
						}
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Start[j] - 1, i]
					}
					else {
						missing_interval <- nas$Start[j]:nas$Stop[j]
						replacement <- input[nas$Stop[j] + 1, i]
					}				
				}
				input[missing_interval, i] <<- replacement
			}
		}
	})
	return(input)
}


#' Load a single raw channel file
#'
#' @param file the path to a raw channel data file.
#' @param date_format the format used in the raw date values (locale dependent)
#' @param skip Number of data rows to skip. Defaults to 0. Note: This is not the
#' 	number of lines to skip from the top of the file. It is the number of lines
#' 	to skip once the actual raw data starts.
#' @param tz The time zone of the data. Defaults to the system time zone.
#'
#' @return A data frame containing the inported data
#'
#' @export
#'
read_pyro_raw_file <- function(file, date_format,
		skip = 0, tz = Sys.timezone()) {

	if (length(file) == 0 || !file.exists(file)) {
		stop("Could not find target file.", call. = FALSE)
	}

	# identify device and channel name
	aux <- readLines(file, n = 10)
	aux <- aux[grepl("^#Device", aux)]
	device <- stringr::str_extract(aux,'(?<= )[^ ]*')
	ch <- stringr::str_extract(file,'(?<=Ch.)[0-9]')

	if (grepl("Oxygen\\.txt$", file) || grepl("pH\\.txt$", file)) {

		base_skip <- if (grepl("Oxygen\\.txt$", file)) 24 else 20

		# grab first row to compile col names
		file_header <- utils::read.table(file, sep = "\t", skip = base_skip,
								 header = FALSE, strip.white = TRUE, nrows = 1)

		# grab rest of rows with the actual data
		output <- utils::read.table(file, sep = "\t",
									skip = base_skip + 1 + skip,
									header = FALSE, strip.white = TRUE)

		# grab first and last words of the headers.
		a <- gsub(" .*$", "", file_header[1,])
		b <- gsub("^.* ", "", file_header[1,])
		b <- gsub(']', '', b)

		# combine first and last words as new column names
		colnames(output) <- tolower(paste0(a, '_', b))

		# Ensure there are no duplicated colum names (data.table shenanigans)
		output <- output[, !duplicated(colnames(output))]
	}

	if (grepl("TempPT100Port\\.txt$", file)) {
		output <- utils::read.table(file, sep = "\t", skip = 12 + skip,
									header = FALSE, strip.white = TRUE)

		colnames(output) <- c('date_main', 'time_main', 'ds', 'temp', 'status')
	}

	# calculate POSIX timestamp
	date_time_aux <- paste(output$date_main, output$time_main)
	date_format <- paste(date_format, "%H:%M:%S")
	output$date_time <- as.POSIXct(date_time_aux, format = date_format, tz = tz)

	# reorganize columns
	if (grepl("Oxygen\\.txt$", file)) {
		output <- output[, c('date_time', 'sample_compt',
							 'pressure_compp', 'oxygen_main')]
		# units
		file_header_string <- paste(file_header, collapse = " ")
		o2_unit <- stringr::str_extract(file_header_string, 
										"(?<=Oxygen \\()[^\\)]*")
		units(output$oxygen_main) <- o2_unit

		pr_unit <- stringr::str_extract(file_header_string, 
										"(?<=Pressure \\()[^\\)]*")
		units(output$pressure_compp) <- pr_unit

		tp_unit <- stringr::str_extract(file_header_string, 
										"(?<=Sample Temp. \\()[^\\)]*")
		units(output$sample_compt) <- tp_unit

		colnames(output)[2:4] <- paste0(c('temp_', 'pressure_', 'ox_'),
										device, ch)
		file_type <- "oxygen"
	}

	if (grepl("pH\\.txt$", file)) {
		output <- output[, c('date_time', 'ph_main')]

		# units
		units(output$ph_main) <- "pH"

		colnames(output)[2] <- paste0(c('ph_'),	device, ch)
		file_type <- "ph"
	}

	if (grepl("TempPT100Port\\.txt$", file)) {
		output <- output[, c('date_time', 'temp')]

		# units
		file_header_string <- paste(file_header, collapse = " ")
		tp_unit <- stringr::str_extract(file_header_string, 
										"(?<=Sample Temp. \\()[^\\)]*")
		units(output$temp_main) <- tp_unit

		colnames(output)[2] <- paste0(c('temp_'), device, ch)

		file_type <- "temp"
	}

	attributes(output)$source_file <- file
	attributes(output)$device <- device
	attributes(output)$ch <- ch
	attributes(output)$file_type <- file_type

	return(output)
}


#' Discard readings
#' 
#' Discard the data from one or more phases of one or more probes.
#' 
#' @param input A computer-friendly data frame.
#' 	The output of \code{\link{melt_resp}} or any downstream function.
#' @param probe The probe(s) from which to discard data. Ommit to discard phases
#' 	from all probes.
#' @param phase The phase(s) to discard.
#' 
#' @return the input data frame without the discarded readings.
#'
#' @export
#'
discard_phase <- function(input, probe, phase) {
	target_phases <- input$cleaned$phase %in% phase 
	if (missing(probe)) {
		target_probes <- rep(TRUE, nrow(input$cleaned))
	} else {
		target_probes <- input$cleaned$probe %in% probe
	}

	to_keep <- !(target_phases & target_probes)

	input$cleaned <- input$cleaned[to_keep, ]
	return(input)
}
