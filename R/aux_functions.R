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
		if (sum(arg_check) < 5) {
			warning("Could not find ", name, " ",
					paste(arg[arg_check], collapse = ", "),
					" in the input. Discarding those.",
					immediate. = TRUE, call. = FALSE)
		} else {
			warning("Could not find ", sum(arg_check), 
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
#' @param d The density. defaults to 1 (1g = 1ml). Defaults to 1.
#'  Higher densities lead to lower volumes, e.g. if d = 2, then 1g = 0.5ml.
#' 
#' @return a units object in ml
#' 
#' @export
#' 
conv_w_to_ml <- function(w, d = 1) {
  if (!inherits(w, "units")) {
    stop("w must be a units class object")
  }
  if (!(as.character(units(w)) %in% c("g", "kg"))) {
    stop("w must be either in grams or kilograms.")
  }

  x <- as.numeric(w) / d

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
#' @param assign_list A list of format device_letter = device_name for the
#'  devices to rename based on their letter.
#' 
#' @return Nothing. Used for side effects.
#' 
#' @export
#' 
assign_device_names <- function(folder, assign_list) {
	if (length(folder) == 0 || !dir.exists(folder)) {
		stop('Could not find target folder')
	}

	if (length(folder) > 1) {
		stop('"folder" should be a string of length 1.')		
	}

	files <- list.files(paste0(folder, '/ChannelData/'))

	file_link <- grepl("Oxygen|pH", files)

	if (all(!file_link)) {
		stop('No probe files found')
	}

	files <- files[file_link]

	capture <- lapply(files, function(i) {
		the_file <- paste0(folder, '/ChannelData/', i)

		x <- suppressWarnings(readLines(the_file))

		r <- grep("^#Device", x)[1]

		# identify device name and letter
		device_name <- stringr::str_extract(x[r],'(?<=Device: )[^\\[]*')
		device_name <- sub(" $", "", device_name)

		device_letter <- 	stringr::str_extract(x[r],'(?<=\\[)[^\\]]*')

		if (device_letter %in% names(assign_list)) {
			x[r] <- sub(device_name, assign_list[[device_letter]], x[r])
			writeLines(x, the_file)
			message("Renamed device ", device_name, " [", device_letter , "] to ",
							assign_list[[device_letter]], " in file '", i, "'.")
		} else {
			message("Could not find match for device ", device_name,
							" [", device_letter , "] in assign_list (file '", i, "'').")
		}
	})
}