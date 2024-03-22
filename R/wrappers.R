#' Wrapper to load the data for an experiment run (pyro and phases files).
#'
#' Scans the target folder for a phases file and the raw pyroscience
#' experiment files. Imports them.
#'
#' @param folder the path to the folder where the experiment data is stored
#' @inheritParams read_pyro_raw_file
#' @param phases_file a pattern to match for a phases file inside the folder
#' @inheritParams load_phases_file
#' @param probe_info a dataframe containing animal information. This
#' 	dataframe must contain the following columns:
#' 		- id     	    : The ID of the animal
#' 		- mass   	    : The mass of the animal, in grams
#' 		- volume 	    : The non-corrected volume of the chamber + tubing
#' 		- probe  	    : The device-channel combination for the probe
#' 		- first_cycle 	: The first cycle of valid data for that animal
#'
#' @return A list containing a phases dataframe and a pyro list with the
#' 	individual source data frames (in source_data), as well as a single,
#'  combined data frame organized by time (in compiled_data).
#'
#' @export
#'
load_experiment <- function(folder, date_format, tz = Sys.timezone(),
		phases_file = "CoolTerm", probe_info, fix_phases = TRUE) {

	if (length(folder) == 0 || !dir.exists(folder)) {
		stop('Could not find target folder')
	}

	if (length(folder) > 1) {
		stop('"folder" should be a string of length 1.')		
	}

	phases_file <- list.files(folder)[grepl(phases_file, list.files(folder))]

	if (length(phases_file) == 0) {
		stop('Could not find phases file')
	} else {
		phases_file <- paste0(folder, "/", phases_file)
	}

	if (!missing(probe_info)) {
		required_cols <- c("id", "mass", "volume", "probe", "first_cycle")
		cols_missing <- !(required_cols %in% colnames(probe_info))
		if (any(cols_missing)) {
			stop("The following required columns are missing ",
				 "from the probe_info input: ",
				 paste0(required_cols[cols_missing], collapse = ", "),
				 call. = FALSE)
		}
	}

	output <- list()

	phases <- lapply(phases_file, function(file) {
		load_phases_file(file, fix_phases = fix_phases)
	})

	names(phases) <- stringr::str_extract(phases_file, '(?<=_)[^_]*(?=.txt)')

	if (any(sapply(names(phases), nchar) > 4)) {
		warning("Long device names detected in the phases input. Are you sure",
				" you appended the device names correctly to the file name?",
				" These are the current device names: ", 
				paste(names(phases), collapse = ", "), ".")
	}

	output$phases <- phases
	output$pyro <- load_pyro_data(folder, date_format = date_format, tz = tz)

	if (!missing(probe_info)) {
		output$probe_info <- probe_info

		units(output$probe_info$mass) <- "g"
		units(output$probe_info$volume) <- "ml"
	}

	return(output)
}


#' Wrapper to scan a pyro folder and load raw data files.
#'
#' @param folder the pyroscience run folder,
#' 	containing a "ChannelData" folder inside
#' @inheritParams read_pyro_raw_file
#' @param type One of "Oxygen" to read only oxygen files, "pH" to read only pH
#' 	files, or "Oxygen|pH" to read both.
#'
#' @export
#'
load_pyro_data <- function(folder, date_format, tz, 
		type = c("Oxygen", "pH", "Oxygen|pH")) {
	type <- match.arg(type)

	files <- list.files(paste0(folder, '/ChannelData/'))

	file_link <- grepl(type, files)

	files <- files[file_link]

	source_data <- lapply(files, function(i) {
		read_pyro_raw_file(paste0(folder, '/ChannelData/', i),
						   date_format = date_format, tz = tz)
	})

	very_start <- min(as.POSIXct(sapply(source_data, function(i) {
		as.character(min(i$date_time))
	})))

	very_end <- max(as.POSIXct(sapply(source_data, function(i) {
		as.character(max(i$date_time))
	})))

	recipient <- data.frame(date_time = seq(from = very_start,
											to = very_end, by = 1))

	for (i in source_data) {
		new_piece <-  i[!duplicated(i$date_time), ]
		recipient <- merge(recipient, new_piece,
						   by = 'date_time', all = TRUE)
	}

	attributes(recipient)$latest_batch_start <- 1

	output <- list(source_data = source_data, compiled_data = recipient)
	return(output)
}


#' Wrapper to get experiment data ready for further analyses
#'
#' Perform standard processing operations to the pyro/phases files.
#' 
#' @param input The output of \code{\link{load_experiment}}
#' @param wait The number of seconds to discard as wait phase
#' @param convert_o2_unit_to The o2 unit desired for the final results
#' @param patch_NAs Logical. Should NA values found in the raw data be patched?
#' 	Defaults to TRUE.
#' @inheritParams patch_NAs
#' @param min_temp,max_temp 
#' 	For temperature ramp experiments. The minimum OR maximum temperatures
#' 	that must be reached before data is considered valid. Discards all phases
#' 	prior to this temperature being reached. Use only one of the two arguments
#' 	at a time.
#' @param start_time,stop_time
#' 	Trim the experiment to a specific time period. You may use one or both of
#' 	these arguments at the same time. Input must be a string in 
#' 	YYYY-MM-DD HH:MM:SS format.
#' @param from_cycle,to_cycle
#' 	Trim the experiment to a specific group of cycles. You may use one or both
#'  of these arguments at the same time. Input must be numeric.
#' @param verbose Logical. Should steps being taken be detailed with messages.
#' 	Defaults to TRUE.
#' 
#' @return An updated experiment list, containing a cleaned object with the
#' 	processed data.
#' 
#' @export 
#'
process_experiment <- function(input, wait, convert_o2_unit_to,
		patch_NAs = TRUE, patch_method = c("linear", "before", "after"),
		min_temp, max_temp, start_time, stop_time, from_cycle, to_cycle, 
		verbose = TRUE) {

	patch_method <- match.arg(patch_method)

	all_units <- c("hPa", "kPa", "torr", "mmHg", "inHg", "mg_per_l", 
				   "ug_per_l", "umol_per_l", "mmol_per_l", "ml_per_l",
				   "mg_per_kg", "ug_per_kg", "umol_per_kg", "mmol_per_kg", 
				   "ml_per_kg")																														

	if (!missing(convert_o2_unit_to) && !(convert_o2_unit_to %in% all_units)) {
		stop("the 'convert_o2_unit_to' argument is not an acceptable unit. ",
			 "Please choose one of the following: ", 
			 paste(all_units, collapse = ", "))
	}

	if (!missing(min_temp) & !missing(max_temp)) {
		stop("Please use only one of 'min_temp' or 'max_temp'",
				 " at a time. See function help for details.")
	}

  	if (verbose) message("M: Merging pyroscience and phases file.")
	input <- merge_pyro_phases(input)
  	if (patch_NAs) {
  		if (verbose) message("M: Patching NA's in the data.")
		input$phased <- patch_NAs(input$phased, patch_method = patch_method, 
								  verbose = FALSE)
	}
	
	if (verbose) message("M: Melting resp data into computer-friendly format")
  	input$melted <- melt_resp(input = input$phased, 
  							  probe_info = input$probe_info)

	if (verbose) message("M: Removing flush and wait values.")
  	input$cleaned <- clean_meas(input = input$melted, wait = wait)

	if (verbose) message("M: Calculating air saturation")

	o2_conv_cols <- c("o2", "temp", "sal", "pressure")
	not_NA <- complete.cases(input$cleaned[, o2_conv_cols])
	original_o2 <- sub("/", "_per_", units(input$cleaned$o2))

	input$cleaned$airsat <- NA
	input$cleaned$airsat[not_NA] <- 
		respirometry::conv_o2(
			o2 = as.numeric(input$cleaned$o2[not_NA]),
			from = original_o2,
			to = "percent_a.s.", 
			temp = as.numeric(input$cleaned$temp[not_NA]), 
			sal = as.numeric(input$cleaned$sal[not_NA]), 
			atm_pres = as.numeric(input$cleaned$pressure[not_NA])
		)
	units(input$cleaned$airsat) <- "percent"

	if (!missing(convert_o2_unit_to)) {
		
		if (verbose) {
			message("M: Converting oxygen unit from ", original_o2, 
					" to ", convert_o2_unit_to, ".")	
		}
		
		# if there is an o2 value but not all the others
		if (any(!is.na(input$cleaned$o2) & !not_NA)) {
			warning_cases <- sum(!is.na(input$cleaned$o2) & !not_NA)
			warning("Invalidating ", warning_cases, "oxygen value(s) as one ",
				"or more of the respective temperature, salinity, or pressure ",
				"values are missing, making it impossible to convert unit.")
			input$cleaned$o2[!not_NA] <- NA
		}

		input$cleaned$o2 <- as.numeric(input$cleaned$o2)
		input$cleaned$o2[not_NA] <- 
			respirometry::conv_o2(
				o2 = input$cleaned$o2[not_NA],
				from = original_o2,
				to = convert_o2_unit_to, 
				temp = as.numeric(input$cleaned$temp[not_NA]), 
				sal = as.numeric(input$cleaned$sal[not_NA]), 
				atm_pres = as.numeric(input$cleaned$pressure[not_NA])
			)
		units(input$cleaned$o2) <- gsub("_per_", "/", convert_o2_unit_to)

		# convert from melted too, which is used for plot_meas
		not_NA <- complete.cases(input$melted[, o2_conv_cols])
		input$melted$o2[!not_NA] <- NA
		input$melted$o2 <- as.numeric(input$melted$o2)
		input$melted$o2[not_NA] <-
			respirometry::conv_o2(
				o2 = input$melted$o2[not_NA],
				from = original_o2,
				to = convert_o2_unit_to, 
				temp = as.numeric(input$melted$temp[not_NA]), 
				sal = as.numeric(input$melted$sal[not_NA]), 
				atm_pres = as.numeric(input$melted$pressure[not_NA])
			)
		units(input$melted$o2) <- gsub("_per_", "/", convert_o2_unit_to)
	}

	if (verbose) message("M: Calculating deltas.")
	input$cleaned <- calc_delta(input$cleaned)

	cutoff <- NULL
	if (!missing(min_temp)) {
		units(min_temp) <- intToUtf8(c(176, 67))
		if (verbose) {
			message(paste0("M: Discarding phases under ", 
						   min_temp, intToUtf8(c(176, 67))))
		}
		cutoff <- head(which(input$cleaned$temp > min_temp), 1)
		if (length(cutoff) == 0) {
			stop ("Temperature never rose above ", min_temp, ".")
		}
	}
	if (!missing(max_temp)) {
		units(max_temp) <- intToUtf8(c(176, 67))
		if (verbose) {
			message(paste0("M: Discarding phases over ", 
						   max_temp, intToUtf8(c(176, 67))))
		}
		cutoff <- head(which(input$cleaned$temp < max_temp), 1)
		if (length(cutoff) == 0) {
			stop ("Temperature never dropped below ", max_temp, ".")
		}
	}
	if (!is.null(cutoff)) {
		the_matches <- which(input$cleaned$phase == input$cleaned$phase[cutoff])
		first_true <- head(the_matches, 1)
		time_break <- input$cleaned$date_time[first_true]
		input$cleaned <- input$cleaned[input$cleaned$date_time >= time_break, ]
	}

	if (!missing(start_time)) {
		if (verbose) {
			message(paste0("M: Discarding phases before ", start_time, "."))
		}
		the_matches <- which(input$cleaned$date_time >= as.POSIXct(start_time))
		cutoff <- head(the_matches, 1)
		if (length(cutoff) == 0) {
			stop ("Data ends before ", start_time, ".")
		} else {
			first_phase <- input$cleaned$phase[cutoff]
			the_matches <- which(input$cleaned$phase == first_phase)
			first_true <- head(the_matches, 1)
			break_ <- input$cleaned$date_time[first_true]
			input$cleaned <- input$cleaned[input$cleaned$date_time >= break_, ]
		}
	}

	if (!missing(stop_time)) {
		if (verbose) {
			message(paste0("M: Discarding phases after ", stop_time, "."))
		}
		the_matches <- which(input$cleaned$date_time <= as.POSIXct(stop_time))
		cutoff <- tail(the_matches, 1)
		if (length(cutoff) == 0) {
			stop ("Data starts after ", stop_time, ".")
		} else {
			last_phase <- input$cleaned$phase[cutoff]
			the_matches <- which(input$cleaned$phase == last_phase)
			last_true <- tail(the_matches, 1)
			break_ <- input$cleaned$date_time[last_true]
			input$cleaned <- input$cleaned[input$cleaned$date_time <= break_, ]
		}
	}

	if (!missing(from_cycle)) {
		if (verbose) {
			message(paste0("M: Discarding cycles prior to cycle ",
						   from_cycle, "."))
		}
		input$cleaned <- input$cleaned[input$cleaned$cycle > from_cycle, ]
	}

	if (!missing(to_cycle)) {
		if (verbose) {
			message(paste0("M: Discarding cycles after cycle ", to_cycle, "."))
		}
		input$cleaned <- input$cleaned[input$cleaned$cycle < from_cycle, ]
	}

	return(input)
}



#' Wrapper to perform metabolic rate calculations
#'
#' Calculates the slopes, filters by threshold R2, calculates metabolic rate for
#' each cycle, and from there calculates various SMR metrics and extracts MMR.
#'
#' @param input The output of \code{\link{subtract_bg}}
#' @inheritParams filter_r2
#' @inheritParams calc_smr
#' 
#' @return an updated input list containing the following new objects:
#' \itemize{
#'  \item \code{all_slopes}: The slopes calculated for each cycle.
#'  \item \code{good_slopes}: The slopes which pass the R2 threshold.
#'  \item \code{mr}: The metabolic rates calculated from the good slopes.
#'  \item \code{smr}: A data frame with the different SMR metrics calculated.
#' 		Relevant details for the different methods are saved in the attributes.
#'  \item \code{mmr}: A data frame containing the cycle with the highest
#' 		metabolic rate recorded.
#' }

#' @export
#'
process_mr <- function(input, r2 = 0.95, G = 1:4, 
					   q = c(0.2, 0.25), p = 0.1, n = 10) {

	input$all_slopes <- calc_slopes(input$cleaned)
	
	input$good_slopes <- filter_r2(input$all_slopes, r2 = r2)

	input$mr <- calc_mr(input$good_slopes)

	# convert seconds to hours (more common)
	the_seconds <- which(units(input$mr$mr_cor)$denominator == "s")
	if (length(the_seconds) > 0) {
		units(input$mr$mr_cor)$denominator[the_seconds] <- "h"
	}
	
	if (is.null(G) && is.null(q) && is.null(p) && is.null(n)) {
		input$smr <- NULL
	} else {
		input$smr <- calc_smr(input$mr, G = G, q = q, p = p, n = n)
		keep_these <- !(colnames(input$smr) %in% c("id", "mass", "volume"))
		smr_aux <- input$smr[, keep_these]
		input$smr <- merge(input$probe_info, smr_aux, 
							 				 by = "probe", all = TRUE)
	}
	
	input$mmr <- extract_mmr(input$mr)
	mmr_aux <- input$mmr[, !(colnames(input$mmr) %in% c("id", "mass", "volume"))]
	input$mmr <- merge(input$probe_info, mmr_aux,
					   				 by = "probe", all = TRUE)

	return(input)
}
