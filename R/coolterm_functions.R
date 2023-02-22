#' import flush and measurement phase times recorded by the program CoolTerm.
#' 
#' @param file The name of the file to be imported.
#' 
#' @return A data.frame containing the start and stop timestamps of the recorded phases.
#' 
#' @export
#' 
load.coolterm.file <- function(file, max.gap.fix = 1) {
	if (!file.exists(file))
		stop("Could not find target file.", call. = FALSE)
	
	phases <- as.data.frame(data.table::fread(file, tz = ""))
	phases$V1 <- as.POSIXct(phases$V1, tz = Sys.timezone())

	breaks <- rle(phases$V2)
	starts <- c(1, cumsum(breaks$lengths) + 1)
	starts <- starts[-length(starts)]
	stops <- cumsum(breaks$lengths)

	x <- data.frame(Phase = breaks$values,
			   Start = phases$V1[starts],
			   Stop = phases$V1[stops])

	if (any(check <- x$Start[-1] <= x$Stop[-nrow(x)])) {
		maxoverlap <- max(as.numeric(difftime(x$Start[which(check) + 1], x$Stop[check])) + 1)
		warning("Found ", sum(check), " overlapping phase(s)! Maximum overlap: ", maxoverlap, " second(s). Saving troublesome rows in attributes", call. = FALSE, immediate. = TRUE)
		attributes(x)$overlaps <- which(check)
	}

	aux <- as.numeric(difftime(x$Start[-1], x$Stop[-nrow(x)])) - 1
	if (any(aux > 0)) {
		warning("Found ", sum(aux > 0), " time gap(s) between phases! Maximum gap: ", max(aux), " second(s). Saving troublesome rows in attributes.", call. = FALSE, immediate. = TRUE)
		auto.gaps <- which(aux > 0 & aux <= max.gap.fix)
		if (length(auto.gaps) > 0) {
			x$Stop[auto.gaps] <- x$Stop[auto.gaps] + aux[auto.gaps]
			message("M: Auto-fixed ", length(auto.gaps), " gaps by extending the previous phase.")
		}
		attributes(x)$gaps <- which(aux > 0)
	}

	aux <- table(rle(x$Phase)$values)

	if (any(aux > 1))
		warning('There are repeated phase names in the output!', immediate. = TRUE, call. = FALSE)

	return(x)
}


#' import flush and measurement phase times recorded by the program CoolTerm.
#' 
#' Updated function for the new pump controller with four individual controls
#' 
#' @param file The name of the file to be imported.
#' 
#' @return A data.frame containing the start and stop timestamps of the recorded phases.
#' 
#' @export
#' 
load_wide_coolterm_file <- function(file, max.gap.fix = 1) {
	if (!file.exists(file))
		stop("Could not find target file.", call. = FALSE)
	
	phases <- as.data.frame(data.table::fread(file, tz = ""))
	colnames(phases)[1] <- "Timestamp"
	phases$Timestamp <- as.POSIXct(phases$Timestamp, tz = Sys.timezone())
	# extract start and stop times for each channel 
	output <- lapply(grep("Phase", colnames(phases)), function(p) {
		breaks <- rle(phases[, p])
		starts <- c(1, cumsum(breaks$lengths) + 1)
		starts <- starts[-length(starts)]
		stops <- cumsum(breaks$lengths)

		data.frame(Phase = breaks$values,
			   	   Start = phases$Timestamp[starts],
	   			   Stop = phases$Timestamp[stops])
	})
	names(output) <- paste0("CH", 1:length(output))

	# check start and stop times for gaps/overlaps
	overlap_instances <- 0
	max_overlap <- 0

	gap_instances <- 0
	max_gap <- 0
	fixed_gaps <- 0
	output <- lapply(output, function(x) {
		if (any(check <- x$Start[-1] <= x$Stop[-nrow(x)])) {
			overlap_instances <<- sum(overlap_instances, check)
			conflict_starts <- x$Start[which(check) + 1]
			conflict_ends <- x$Stop[which(check)]
			overlaps <- as.numeric(difftime(conflict_starts, conflict_ends)) + 1
			max_overlap <<- max(max_overlap, overlaps)
		}

		aux <- as.numeric(difftime(x$Start[-1], x$Stop[-nrow(x)])) - 1
		if (any(aux > 0)) {
			gap_instances <<- sum(gap_instances, aux > 0)
			max_gap <<- max(max_gap, aux)
			
			auto.gaps <- which(aux > 0 & aux <= max.gap.fix)
			fixed_gaps <<- sum(fixed_gaps, length(auto.gaps))

			if (length(auto.gaps) > 0) {
				x$Stop[auto.gaps] <- x$Stop[auto.gaps] + aux[auto.gaps]
			}
		}

		aux <- table(rle(x$Phase)$values)

		if (any(aux > 1))
			warning('There are repeated phase names in the output!', immediate. = TRUE, call. = FALSE)

		return(x)
	})

	if (overlap_instances > 0) {
		warning("Found ", overlap_instances, " overlapping phase(s)! Maximum overlap: ", max_overlap, " second(s).", call. = FALSE, immediate. = TRUE)
	}

	if (gap_instances > 0) {
		warning("Found ", gap_instances, " time gap(s) between phases! Maximum gap: ", max_gap, " second(s).", call. = FALSE, immediate. = TRUE)
		message("M: Auto-fixed ", fixed_gaps, " gaps by extending the previous phase.")
	}

	return(output)
}



#' Merge the oxygen data obtained from the pyroscience software with the phases recorded by Coolterm.
#' 
#' @param pyroscience The dataframe containing the O2 measurements and imported using \code{\link{load.pyroscience.logger.file}} or \code{\link{load.pyroscience.workbench.file}}.
#' @param coolterm The dataframe containing the phase information, as imported by the function \code{\link{load.coolterm.file}}
#' 
#' @return a data.frame similar to the pyroscience input, but where the Phase column has been updated
#' 
#' @export
#' 
merge_pyroscience_coolterm <- function(pyroscience, coolterm) {
	pyroscience$Phase <- as.character(pyroscience$Phase)

	pyroscience$Phase[pyroscience$Date.Time < coolterm$Start[1]] <- "F0"

	for (i in 1:nrow(coolterm)) {
		this.phase <- pyroscience$Date.Time >= coolterm$Start[i] & pyroscience$Date.Time <= coolterm$Stop[i]
		pyroscience$Phase[this.phase] <- coolterm$Phase[i]
	}

	n.phases <- max(as.numeric(sub("[F|M]", "", coolterm$Phase)))
	pyroscience$Phase[pyroscience$Date.Time > coolterm$Stop[nrow(coolterm)]] <- paste0("F", n.phases + 1)

	if (any(check <- is.na(pyroscience$Phase)))
		warning(sum(check), " measurement(s) could not be assigned to a phase!", call. = FALSE, immediate. = TRUE)

	pyroscience$Phase <- factor(pyroscience$Phase, levels = unique(pyroscience$Phase))

	return(pyroscience)
}


#' Merge the oxygen data obtained from the pyroscience software with the phases recorded by Coolterm.
#' 
#' updated to work with the four channel pump controller
#' 
#' @param pyroscience The dataframe containing the O2 measurements and imported using \code{\link{load.pyroscience.logger.file}} or \code{\link{load.pyroscience.workbench.file}}.
#' @param coolterm The dataframe containing the phase information, as imported by the function \code{\link{load.coolterm.file}}
#' 
#' @return a data.frame similar to the pyroscience input, but where the Phase column has been updated
#' 
#' @export
#' 
merge_pyroscience_wide_coolterm <- function(pyroscience, coolterm) {
	pyroscience$Phase <- NULL

	new_col_order <- 1

	tmp <- lapply(1:length(coolterm), function(device) {
		lapply(1:length(coolterm[[device]]), function(probe) {
			# create column with placeholders
			pyroscience[, paste0("Phase.", names(coolterm)[device], probe)] <- "F0"

			# assign phases
			# cat(device) # debug messages
			# cat(probe) # debug messages
			for (i in 1:nrow(coolterm[[device]][[probe]])) {
				this.phase <- pyroscience$Date.Time >= coolterm[[device]][[probe]]$Start[i] & pyroscience$Date.Time <= coolterm[[device]][[probe]]$Stop[i]
				pyroscience[this.phase, paste0("Phase.", names(coolterm)[device], probe)] <- coolterm[[device]][[probe]]$Phase[i]
			}

			if (any(check <- is.na(pyroscience[, paste0("Phase.", names(coolterm)[device], probe)])))
				warning(sum(check), " measurement(s) in device ",device, ", probe ", probe, " could not be assigned to a phase!", call. = FALSE, immediate. = TRUE)

			pyroscience[, paste0("Phase.", names(coolterm)[device], probe)] <- factor(pyroscience[, paste0("Phase.", names(coolterm)[device], probe)], levels = unique(pyroscience[, paste0("Phase.", names(coolterm)[device], probe)]))
		
			# export pyroscience object to outside of the lapply loop
			pyroscience <<- pyroscience

			new_col_order <<- c(new_col_order, grep(paste0(names(coolterm)[device], probe), colnames(pyroscience)))
		})
		# export pyroscience object to outside of the lapply loop
		pyroscience <<- pyroscience
		new_col_order <<- new_col_order
	})
	rm(tmp)
	
	pyroscience <- pyroscience[, new_col_order]

	return(pyroscience)
}
