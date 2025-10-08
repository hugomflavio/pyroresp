#' suppress incomplete final line warning
#'
#' Data files do not end with an EOL. This is by design, so a warning should
#' not be thrown. This function catches that specific warning and suppresses
#' it, while letting other warnings be shown.
#'
#' @param expr the read.table function call
#'
#' @keywords internal
#'
#' @return The outcome of the expression
#'
suppress_EOL_warning <- function(expr) {
	.suppress_eol <- function(w) {
		if (!grepl("incomplete final line", w$message)) {
			warning(w)
		}
		tryInvokeRestart("muffleWarning")
	}
	x <- withCallingHandlers(expr, warning = function(w) .suppress_eol(w))
	return(x)
}

#' check that the requested argument value is present in the data
#'
#' @param arg the argument value
#' @param data A vector of the target data
#' @param name For messaging purposes; the type of data we're working with.
#' @param verbose Should a warning be sent if the argument does not fully
#' 	match the data?
#'
#' @return the matching elements of arg
#'
#' @keywords internal
#'
check_arg_in_data <- function(arg, data, name, verbose = TRUE) {
	arg_check <- arg %in% data

	if (all(!arg_check)) {
		stop("Could not find requested ", name, " in the input", call. = FALSE)
	}
	if (verbose && any(!arg_check)) {
		if (sum(!arg_check) < 5) {
			warning("Could not find ", name, " ",
					paste(arg[!arg_check], collapse = ", "),
					" in the input. Discarding those.",
					immediate. = TRUE, call. = FALSE)
		} else {
			warning("Could not find ", sum(!arg_check),
					" of the requested ", name, " ",
					"in the input. Discarding those.",
					immediate. = TRUE, call. = FALSE)
		}
	}

	arg <- arg[arg_check]
	return(arg)
}

#' convert weight to volume in ml
#'
#' @param w a units object in either grams or kilograms
#' @param density The density. defaults to 1 (1 g = 1 ml).
#'  Higher densities lead to lower volumes, e.g. if density = 2,
#' 	then 1 g = 0.5 ml.
#'
#' @return a units object in ml
#'
#' @export
#'
conv_w_to_ml <- function(w, density = 1) {
  if (!inherits(w, "units")) {
    stop("w must be a units class object")
  }
  if (!(as.character(units(w)) %in% c("g", "kg"))) {
    stop("w must be either in grams or kilograms.")
  }

  x <- as.numeric(w) / density

  if (as.character(units(w)) == "kg") {
    x <- x/1000
  }

  units(x) <- "ml"

  return(x)
}


#' Assign device names after completing experiment
#'
#' To use if you forgot to assign a device name at the start of the experiment.
#' Goes through the raw data files, renames the required device names, and
#' resaves the files.
#'
#' NOTE: This function will modify the files in your data folder!
#'
#' @inheritParams load_experiment
#' @param old_name The name of the device to be substituted.
#' @param new_name The new name to assign to the device.
#' @param old_letter The letter of the device to be substituted. To be used
#' 	instead of old_name, if preferable.
#' @param confirmed Logical: Have you reviewed the changes and confirm you
#'  want to apply them? Defaults to FALSE to perform a dry run.
#'
#' @return Nothing. Used for side effects.
#'
#' @export
#'
assign_device_names <- function(folder, old_name, new_name, old_letter,
																confirmed = FALSE, encoding = "ISO-8859-1") {
	found_none <- TRUE

	if (length(folder) != 1) {
		stop('"folder" should be a string of length 1.')
	}

	if (!dir.exists(folder)) {
		stop('Could not find target folder')
	}

	if (!missing(old_name) & !missing(old_letter)) {
		stop("Use only one of old_name or old_letter")
	}

	files <- list.files(paste0(folder, '/ChannelData/'))

	file_link <- grepl("Oxygen|pH", files)

	if (all(!file_link)) {
		stop('No probe files found in specified folder')
	}

	files <- files[file_link]

	capture <- lapply(files, function(i, old_name, new_name, old_letter) {
		the_file <- paste0(folder, '/ChannelData/', i)

		x <- readLines(the_file, warn = FALSE)
		x <- stringr::str_conv(x, encoding = encoding)

		r <- grep("^#Device", x)[1]

		# identify old device name and letter
		old_device_name <- stringr::str_extract(x[r],'(?<=Device: )[^\\[]*')
		old_device_name <- sub(" $", "", old_device_name)

		old_device_letter <- 	stringr::str_extract(x[r],'(?<=\\[)[^\\]]*')

		name_check <- !missing(old_name) && old_name == old_device_name
		letter_check <- !missing(old_letter) && old_letter == old_device_letter
		if (name_check | letter_check) {
			found_none <<- FALSE
			x[r] <- sub(old_device_name, new_name, x[r])
			if (confirmed) {
				writeLines(x, the_file)
				message("Renamed device '", old_device_name, "' [", old_device_letter,
					"] to '",	new_name, "' [", old_device_letter,
					"] in file '", i, "'.")
			} else {
				message("Would have renamed device '", old_device_name, "' [",
					old_device_letter , "] to '", new_name,
					"' [", old_device_letter,
					"] in file '", i, "'.")
			}
		}
	}, old_name = old_name, new_name = new_name, old_letter = old_letter)

	if (found_none) {
		message("No devices were found that matched the old name or letter.")
	}

	if (!confirmed & !found_none) {
		warning("No changes were made. Review the messages above.\n",
						"If you wish to apply those changes, run again with ",
						"confirmed = TRUE.", call. = FALSE)
	}
}

#' helper function to avoid attribute loss
#' caused by the split->lapply->rbind process
#'
#' @param from The original data frame
#' @param to The remade data frame
#'
#' @keywords internal
#'
#' @return "to", with any new attributes transferred from "from."
transfer_attributes <- function(from, to) {
	attr_from <- names(attributes(from))
	attr_to <- names(attributes(to))

	to_transfer <- attr_from[!(attr_from %in% attr_to)]

	if (length(to_transfer) > 0) {
		for (i in to_transfer) {
			attributes(to)[i] <- attributes(from)[i]
		}
	}
	return(to)
}

#' calculate standard error of the mean
#'
#' @param x a vector of numbers
#' @param na.rm Defaults to TRUE
#'
#' @return the SEM value
#'
#' @export
#'
sem <- function(x, na.rm = TRUE){
    a <- length(x)

    if(na.rm) {
        x <- x[!is.na(x)]
    }

    if(a != length(x)) {
        message("M: Omitted ", a - length(x), " missing values.")
		}

    output <- sd(x) / sqrt(length(x))

    return(output)
}


#' Calculate area under the curve for two vectors
#'
#' @param x a vector of x values
#' @param y a vector of y values
#' @param zero which value to consider as the y floor. defaults to 0.
#'
#' @return a dataframe containing the y and x values, as well as the calculated
#' 	auc for each pair of points and the cumulative auc as X progresses.
#'
#' @export
#'
auc <- function(x, y, zero = 0) {
	units(zero) <- units(y)
	if (length(x) != length(y)) {
		stop("x and y must have same length")
	}
	aux <- as_units(0, value = units(x))
	df <- data.frame(x = c(aux, x),
	                 y = y[c(1, 1:length(y))] - zero)
	df <- df[order(x), ]
	df$auc <- 0

	for (i in 2:nrow(df)) {
		y_avg <- (df$y[i] + df$y[i-1])/2 # e.g. km/h
		x_change <- df$x[i] - df$x[i-1] # e.g. h
		df$auc[i] <- y_avg * x_change # cuts the h, so only km is left
	}

	# transfer units
	dummy <- df$y[1] * df$x[1]
	units(df$auc) <- units(dummy)

	df$cumauc <- cumsum(df$auc)

	return(df)
}
