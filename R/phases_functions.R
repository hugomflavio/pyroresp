#' import flush and measurement phase times recorded by the program CoolTerm.
#' 
#' Updated function for the new pump controller with four individual controls
#' 
#' @param file The name of the file to be imported.
#' @param fix_phases Automatically correct any overlaps, gaps, or repeated
#' 	phase names found in the file.
#' 
#' @return A data frame containing the start and 
#' 	stop timestamps of the recorded phases.
#' 
#' @export
#' 
load_phases_file <- function(file, fix_phases) {
	if (!file.exists(file))
		stop("Could not find target file.", call. = FALSE)
	
	phases <- as.data.frame(data.table::fread(file, tz = ""))
	colnames(phases) <- tolower(colnames(phases))
	colnames(phases)[1] <- "timestamp"
	phases$timestamp <- as.POSIXct(phases$timestamp, tz = Sys.timezone())

	# extract start and stop times for each channel 
	output <- lapply(grep("phase", colnames(phases)), function(p) {
		breaks <- rle(phases[, p])
		starts <- c(1, cumsum(breaks$lengths) + 1)
		starts <- starts[-length(starts)]
		stops <- cumsum(breaks$lengths)

		data.frame(phase = breaks$values,
			   	   start = phases$timestamp[starts],
	   			   stop = phases$timestamp[stops])
	})
	names(output) <- paste0("ch", 1:length(output))

	output <- check_phases(output, fix_phases = fix_phases)

	attributes(output)$nlines <- nrow(phases) + 1
	attributes(output)$sourcefile <- file

	return(output)
}

#' Check if phases line up correctly and are not repeated
#' 
#' @param phases The phases list. The "phases" object in the list generated
#'  by \code{\link{load_experiment}}
#' @inheritParams load_phases_file
#' 
#' @return The phases list, with troublesome rows saved to the attributes of
#'  each probe's data frame.
#' 
#' @export
#' 
check_phases <- function(phases, fix_phases = FALSE) {
	# check start and stop times for gaps/overlaps
	overlap_instances <- 0
	max_overlap <- 0
	gap_instances <- 0
	max_gap <- 0
	repeat_instances <- 0

	output <- lapply(phases, function(x) {
		# The overlap vs gap code below could probably be optimized into
		# a single check.

		# OVERLAPS
		# find if phase starts before previous phase ends.
		overlap_check <- x$start[-1] <= x$stop[-nrow(x)]
		# the comparison above removes the first start and the last stop,
		# so that we start by comparing second start the first stop, and 
		# end by comparing the last start to the before-last stop.
		# If the first value is TRUE, then the second phase overlaps with
		# the first. 
		if (any(overlap_check)) {
			# add to counter
			overlap_instances <<- sum(overlap_instances, overlap_check)
			# Find how big the overlaps are.
			conflict_starts <- x$stop[which(overlap_check)]
			conflict_ends <- x$start[which(overlap_check) + 1]
			overlaps <- as.numeric(difftime(conflict_starts, conflict_ends)) + 1
			# add one second to the difftime as a time against itself results
			# in 0 seconds difference (i.e. 1 second overlap)

			# update max overlap if relevant
			max_overlap <<- max(max_overlap, overlaps)
			
			attributes(x)$overlaps <- which(overlap_check) + 1
			# +1 to blame the phase below. E.g. if the check is true, then
			# the second phase overlaps with the first (blame the second)
		}

		# GAPS
		# Find if there are gap seconds between 
		# a phase ending and another starting.
		gap_check <- x$start[-1] > (x$stop[-nrow(x)] + 1)
		# See explanations above for -1 and -nrow(x)
		# +1 to the stop so it should match the next start.
		if (any(gap_check)) {
			# update counter
			gap_instances <<- sum(gap_instances, gap_check)
			# Find how big the gaps are
			conflict_starts <- x$stop[which(gap_check)]
			conflict_ends <- x$start[which(gap_check) + 1]
			gaps <- as.numeric(difftime(conflict_ends, conflict_starts)) - 1
			# remove 1 second from the difftime as a time against the very next
			# second results in 1 second difference (i.e. 0 second gap)

			# update max gap if relevant
			max_gap <<- max(max_gap, gaps)

			attributes(x)$gaps <- which(gap_check) + 1
			# +1 to blame the phase below. E.g. if the check is true, then
			# the second phase overlaps with the first (blame the second)
		}

		# REPEATS
		repeat_check <- duplicated(x$phase)
		if (any(repeat_check)) {
			repeat_instances <- sum(repeat_instances, repeat_check)
			attributes(x)$repeated <- which(repeat_check)
		}

		return(x)
	})

	if (overlap_instances > 0) {
		if (fix_phases) {
			output <- fix_phase_overlaps(output)
		} else {
			warning("Found ", overlap_instances, 
				" overlapping phase(s)! Maximum overlap: ", max_overlap, 
				" second(s). This must be fixed before proceeding.",
				call. = FALSE, immediate. = TRUE)
			message("Overlaps may be automatically",
					" trimmed using fix_phase_overlaps().")
		}
	}

	if (gap_instances > 0) {
		if (fix_phases) {
			output <- fix_phase_gaps(output)
		} else {
			warning("Found ", gap_instances, 
				" gap(s) between phases! Maximum gap: ", max_gap, 
				" second(s). This must be fixed before proceeding.",
				call. = FALSE, immediate. = TRUE)
			message("Gaps may be automatically bridged using fix_phase_gaps().")
		}
	}

	if (repeat_instances > 0) {
		if (fix_phases) {
			output <- rename_phases(output)
		} else {
			warning("Found ", repeat_instances, 
				" repeated phase name(s)!. This must be fixed before proceeding.",
				call. = FALSE, immediate. = TRUE)
			message("Phase names may be automatically",
					" renamed using rename_phases().")
		}
	}

	return(output)
}

#' Fix phase overlaps
#' 
#' Automatically trims the end of phases so they do not
#' overlap with the beginning of the next ones.
#' 
#' @param input A list containing imported experiment data.
#' 	The output from \code{\link{load_experiment}}.
#' 
#' @return The input list with a corrected phases list.
#' 
#' @export
#' 
fix_phase_overlaps <- function(input) {
	# chose right input depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		aux <- list(dvc = input)
	} else {
		aux <- input$phases
	}

	fixed_overlaps <- 0

	output <- lapply(aux, function(dvc) {
		out2 <- lapply(dvc, function(prb) {
			# see check_phases for code explanations
			overlap_check <- prb$start[-1] <= prb$stop[-nrow(prb)]

			if (any(overlap_check)) {
				to_fix <- which(overlap_check)
				prb$stop[to_fix] <- prb$start[to_fix + 1] - 1
				fixed_overlaps <<- sum(fixed_overlaps, overlap_check)
			}

			return(prb)
		})
		return(out2)
	})

	if (fixed_overlaps > 0) {
		message("M: Fixed ", fixed_overlaps, 
				" overlaps by trimming the previous phase.")
	} else {
		message("M: No overlaps found.")
	}

	# chose right output depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		output <- output[[1]]
	} else {
		input$phases <- output
		output <- input
	}

	return(output)
}

#' Fix phase gaps
#' 
#' Automatically expands the end of phases so they end
#' just before the start of the next ones.
#' 
#' @param input A list containing imported experiment data.
#' 	The output from \code{\link{load_experiment}}
#' 
#' @return The input list with a corrected phases list.
#' 
#' @export
#' 
fix_phase_gaps <- function(input) {
	# chose right input depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		aux <- list(dvc = input)
	} else {
		aux <- input$phases
	}

	fixed_gaps <- 0

	output <- lapply(aux, function(dvc) {
		out2 <- lapply(dvc, function(prb) {
			# see check_phases for code explanations
			gap_check <- prb$start[-1] > (prb$stop[-nrow(prb)] + 1)

			if (any(gap_check)) {
				to_fix <- which(gap_check)
				prb$stop[to_fix] <- prb$start[to_fix + 1] - 1
				fixed_gaps <<- sum(fixed_gaps, gap_check)
			}

			return(prb)
		})
		return(out2)
	})

	if (fixed_gaps > 0) {
		message("M: Fixed ", fixed_gaps, 
				" gaps by extending the previous phase.")
	} else {
		message("M: No gaps found.")
	}

	# chose right output depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		output <- output[[1]]
	} else {
		input$phases <- output
		output <- input
	}

	return(output)
}

#' Rename phases from scratch
#' 
#' Replaces recorded phase names with new ones automatically starting from
#' M1 and F1. If the first phase is a flush, then that one gets renamed
#' to F0 instead.
#' 
#' @param input A list containing imported experiment data.
#' 	The output from \code{\link{load_experiment}}
#' 
#' @return The input list with a corrected phases list.
#' 
#' @export
#' 
rename_phases <- function(input) {
	# chose right input depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		aux <- list(dvc = input)
	} else {
		aux <- input$phases
	}

	output <- lapply(aux, function(dvc) {
		out2 <- lapply(dvc, function(prb) {
	
			aux <- substr(prb$phase, 1, 1)
			Fs <- aux == "F"
			if (Fs[1]) {
				prb$phase[Fs] <- paste0("F", 0:(sum(Fs) - 1))
			} else {
				prb$phase[Fs] <- paste0("F", 1:sum(Fs))
			}

			Ms <- aux == "M"
			prb$phase[Ms] <- paste0("M", 1:sum(Ms))			

			return(prb)
		})
		return(out2)
	})

	message("M: Phase names reset (starting at 1).")

	# chose right output depending on if function was
	# called by check_phases or by the user.
	if (is.null(input$phases)) {
		output <- output[[1]]
	} else {
		input$phases <- output
		output <- input
	}

	return(output)
}

#' Merge the pyro data with the phases.
#' 
#' Matches the two data sources by probe and timestamp and assigns the 
#' respective phases to each datapoint. The phases object must not have
#' gaps/overlaps/repeated phase names for downstream functions to work
#' properly.
#' 
#' @param input A list containing imported experiment data.
#' 	The output from \code{\link{load_experiment}}
#' 
#' @return An updated experiment list containing a phased data frame
#' 
#' @export
#' 
merge_pyro_phases <- function(input) {
	pyrodata <- input$pyro$compiled_data
	pyrodata$phase <- NULL

	phases <- input$phases

	new_col_order <- 1

	tmp <- lapply(1:length(phases), function(dvc) {
		lapply(1:length(phases[[dvc]]), function(prb) {
			# create column with placeholders
			new_col <- paste0("phase_", names(phases)[dvc], prb)
			pyrodata[, new_col] <- "F0"

			# assign phases
			for (i in 1:nrow(phases[[dvc]][[prb]])) {
				check1 <- pyrodata$date_time >= phases[[dvc]][[prb]]$start[i]
				check2 <- pyrodata$date_time <= phases[[dvc]][[prb]]$stop[i]
				this.phase <- check1 & check2

				pyrodata[this.phase, new_col] <- phases[[dvc]][[prb]]$phase[i]
			}

			NA_check <- is.na(pyrodata[, new_col])
			if (any(NA_check)) {
				warning(sum(NA_check), " measurement(s) in device ",dvc,
						", probe ", prb, " could not be assigned to a phase!",
						call. = FALSE, immediate. = TRUE)
			}

			pyrodata[, new_col] <- factor(pyrodata[, new_col], 
										  levels = unique(pyrodata[, new_col]))
		
			# export pyrodata object to outside of 
			# the lapply loop to save changes
			pyrodata <<- pyrodata

			# export location of columns pertaining to this
			# probe so they can all be grouped in final output
			dvc_prb <- paste0(names(phases)[dvc], prb)
			relevant_columns <- grep(dvc_prb, colnames(pyrodata))
			new_col_order <<- c(new_col_order, relevant_columns)
		})
		# export pyrodata object to outside of the lapply loop
		pyrodata <<- pyrodata
		new_col_order <<- new_col_order
	})
	rm(tmp)
	
	pyrodata <- pyrodata[, new_col_order]

	input$phased <- pyrodata

	return(input)
}


#' Replicate phases from one probe to others
#' 
#' Useful when using one flush pump for many chambers/probes
#' 
#' @param input A list containing imported experiment data.
#' 	The output from \code{\link{load_experiment}}
#' @param from_device Name of the device from which to copy
#' @param from_channel Name of the channel within "from_device" to duplicate
#' @param to_device Name of the device to which to copy (can be the same)
#' @param to_channel Name of the channel (or channels if vector) 
#'  to which to replicate.
#' 
#' @return The input list with a corrected phases list.
#' 
#' @export
#' 
replicate_phases <- function(input, 
		from_device, from_channel, to_device, to_channel){
	if(is.null(input$phases[[from_device]][[from_channel]]))
		stop("Could not find specified channel in device ", from_device)

	replacement <- input$phases[[from_device]][[from_channel]]

	for (todev in to_device) {
		for (toCH in to_channel) {
			input$phases[[todev]][[toCH]] <- replacement	 
		}
	}

	return(input)
}