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
	phases$V2 <- factor(phases$V2, levels = unique(phases$V2))

	aux <- split(phases, phases$V2)

	x <- lapply(aux, function(x) {
		data.frame(Start = x$V1[1], Stop = x$V1[nrow(x)])
	})

	x <- as.data.frame(data.table::rbindlist(x))
	x$Phase <- names(aux)
	x <- x[, c("Phase", "Start", "Stop")]

	substrRight <- function(x, n){
	  substr(x, nchar(x)-n+1, nchar(x))
	}

	if (any(check <- x$Start[-1] <= x$Stop[-nrow(x)])) {
		maxoverlap <- max(as.numeric(difftime(x$Start[which(check) + 1], x$Stop[check])) + 1)
		warning("Found ", sum(check), " overlapping phase(s)! Maximum overlap: ", maxoverlap, " second(s). Saving troublesome rows in attributes", call. = FALSE, immediate. = TRUE)
		attributes(x)$overlaps <- which(check)
	}

	aux <- as.numeric(difftime(x$Start[-1], x$Stop[-nrow(x)])) - 1
	if (any(aux > 0)) {
		warning("Found ", sum(aux > 0), " gap(s) between phases! Maximum gap: ", max(aux), " second(s). Saving troublesome rows in attributes.", call. = FALSE, immediate. = TRUE)
		auto.gaps <- which(aux > 0 & aux <= max.gap.fix)
		if (length(auto.gaps) > 0) {
			x$Stop[auto.gaps] <- x$Stop[auto.gaps] + aux[auto.gaps]
			message("M: Auto-fixed ", length(auto.gaps), " gaps by extending the previous phase.")
		}
		attributes(x)$gaps <- which(aux > 0)
	}

	return(x)
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
