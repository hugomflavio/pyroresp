#' Load a single raw channel file
#'
#' @param file the path to a raw channel data file.
#' @param date_format the format used in the raw date values (locale dependent)
#' @param skip Number of data rows to skip. Defaults to 0. Note: This is not the
#'   number of lines to skip from the top of the file. It is the number of lines
#'   to skip once the actual raw data starts.
#' @param tz The time zone of the data. Defaults to the system time zone.
#' @param encoding the encoding of the source file. Defaults to "iso-8859-1",
#'   which is the encoding used by the pyro workbench software.
#'
#' @return A data frame containing the inported data
#'
#' @export
#'
read_pyro_raw_file <- function(file, date_format,
    skip = 0, tz = Sys.timezone(), encoding = "ISO-8859-1") {

  if (length(file) == 0 || !file.exists(file)) {
    stop("Could not find target file.", call. = FALSE)
  }

  if(!(encoding %in% stringi::stri_enc_list())) {
    stop("encoding argument not recognized. See stringi::stri_enc_list()",
         " for a complete list of recognized values")
  }

  # identify device and channel name
  file_header <- readLines(file, n = 50, warn = FALSE)
  file_header <- stringr::str_conv(file_header, encoding = encoding)
  device_line <- file_header[grepl("^#Device", file_header)]
  device <- stringr::str_extract(device_line,'(?<=Device: )[^\\[]*')
  device <- sub(" $", "", device)
  ch <- stringr::str_extract(file,'(?<=Ch.)[0-9]')
  table_start <- grep("#--- Measurement Data", file_header)

  # grab first row to compile col names
  table_header <- suppress_EOL_warning(
                    read.table(file, sep = "\t", skip = table_start,
                               header = FALSE, strip.white = TRUE,
                               nrows = 1, fileEncoding = encoding)
                  )

  # grab rest of rows with the actual data
  output <- suppress_EOL_warning(
              read.table(file, sep = "\t",
                         skip = table_start + 1 + skip,
                         header = FALSE, strip.white = TRUE,
                         fileEncoding = encoding)
            )

  if (grepl("Oxygen\\.txt$", file) || grepl("pH\\.txt$", file)) {

    cal_line <- grep("lastCal", file_header)
    cal_headers <- unlist(strsplit(file_header[cal_line], "\t"))[-1]
    cal_settings <- unlist(strsplit(file_header[cal_line + 1], "\t"))[-1]
    names(cal_settings) <- cal_headers

    # grab first and last words of the table headers.
    a <- gsub(" .*$", "", table_header[1,])
    b <- gsub("^.* ", "", table_header[1,])
    b <- gsub(']', '', b)

    # combine first and last words as new column names
    colnames(output) <- tolower(paste0(a, '_', b))

    # Ensure there are no duplicated column names (data.table shenanigans)
    output <- output[, !duplicated(colnames(output))]

    # extract salinity
    settings_line <- grep("--- Settings & Calibration", file_header)
    sal_value <- stringr::str_extract(file_header[settings_line + 2],
                                      "(?<=\t)[^\t]*$")

    output$salinity <- as.numeric(sal_value)
  }

  if (grepl("TempPT100Port\\.txt$", file)) {

    cal_line <- grep("Calibration offset", file_header)
    cal_settings <- stringr::str_extract(file_header[cal_line],
                                         "(?<=Calibration offset = )[0-9|\\.]*")
    names(cal_settings) <- "Calibration offset"

    colnames(output) <- c('date_main', 'time_main', 'ds', 'temp', 'status')
  }

  # calculate POSIX timestamp
  date_time_aux <- paste(output$date_main, output$time_main)
  date_format <- paste(date_format, "%H:%M:%S")
  output$date_time <- as.POSIXct(date_time_aux, format = date_format, tz = tz)

  # reorganize columns and assign units
  if (grepl("Oxygen\\.txt$", file)) {
    output <- output[, c('date_time', 'sample_compt',
                         'pressure_compp', 'salinity', 'oxygen_main')]
    # assign units
    table_header_string <- paste(table_header, collapse = " ")
    tp_unit <- stringr::str_extract(table_header_string,
                                    "(?<=Sample Temp. \\()[^\\)]*")
    units(output$sample_compt) <- tp_unit

    pr_unit <- stringr::str_extract(table_header_string,
                                    "(?<=Pressure \\()[^\\)]*")
    units(output$pressure_compp) <- pr_unit

    # notice that the salinity unit is grabbed from a different place.
    sal_unit <- stringr::str_extract(file_header[settings_line + 1],
                                     "(?<=Salinity \\()[^\\)]*")
    units(output$salinity) <- sal_unit

    o2_unit <- stringr::str_extract(table_header_string,
                                    "(?<=Oxygen \\()[^\\)]*")

    if (grepl("%", o2_unit)) {
      o2_unit <- "%"
    }
    units(output$oxygen_main) <- o2_unit

    colnames(output)[2:5] <- paste0(c('temp_', 'pressure_', 'sal_', 'ox_'),
                                    device, ch)
    file_type <- "oxygen"
  }

  if (grepl("pH\\.txt$", file)) {
    output <- output[, c('date_time', 'ph_main')]
    # assign unit
    units(output$ph_main) <- "pH"

    colnames(output)[2] <- paste0(c('ph_'),  device, ch)
    file_type <- "ph"
  }

  if (grepl("TempPT100Port\\.txt$", file)) {
    output <- output[, c('date_time', 'temp')]
    # assign unit
    table_header_string <- paste(table_header, collapse = " ")
    tp_unit <- stringr::str_extract(table_header_string,
                                    "(?<=Sample Temp. \\()[^\\)]*")
    units(output$temp) <- tp_unit

    colnames(output)[2] <- paste0(c('temp_'), device)
    file_type <- "temp"
  }

  attributes(output)$source_file <- file
  attributes(output)$device <- device
  attributes(output)$ch <- ch
  attributes(output)$file_type <- file_type
  attributes(output)$cal_settings <- cal_settings

  return(output)
}

#' Fill in missing datapoints using the nearest available data.
#'
#' At high sampling rates, data points can occasionally be lost, causing odd
#' interruptions in the temperature and/or oxygen traces.
#' This function fills in those gaps.
#'
#' @param input A dataframe containing Timestamps on the first column and
#'   matching data on the remaining. If input contains a "phase" column, it will
#'   be transported to the output unchanged.
#' @param patch_method One of:
#'   'linear' to capture the nearest before and after non-NA values and make a
#'        linear interpolation.
#'   'before' to find the nearest non-NA value before the NA and use it to fill
#'        the gap.
#'   'after'  to do the same as above but with the nearest value coming after
#'        the NA.
#' @param verbose Logical. Should exceptions be reported out loud.
#'   Defaults to TRUE.
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
    breaks <- data.frame(value = rle_list[[i]]$values,
                         start = c(1, aux[-length(aux)] + 1),
                         stop = aux)

    # if there are any breaks, start working on them
    if (any(breaks$value)) {
      nas <- breaks[breaks$value, ]

      # check for wide gaps. Complain if needed
      nas$n <- nas$stop - nas$start + 1
      if (any(nas$n > 5)) {
        ngaps <- sum(nas$n > 5)
        warning("Found ",
                ifelse(ngaps == 1,
                     "a wide gap",
                     paste(ngaps, "wide gaps")),
                " of NAs in column ", i, " (n > 5).",
                " Won't auto-fix wide gaps. Check data manually.",
                immediate. = TRUE, call. = FALSE)
        nas <- nas[nas$n <= 5, ]
      }

      if (nrow(nas) > 0) {
        input <<- aux_fun_process_NAs(nas = nas, input = input,
                                      patch_method = patch_method,
                                      col = i, verbose = verbose)
      }
    }
  })
  return(input)
}

#' auxiliary function of patch_NAs
#'
#' Internal. Performs the actual patching for one column at a time.
#'
#' @param nas A table with the start and stop rows of the NA intervals
#' @param col the column being patched
#' @inheritParams patch_NAs
#'
#' @keywords internal
#'
aux_fun_process_NAs <- function(nas, col, input, patch_method, verbose) {
  # for every break found, apply the correction
  for (j in 1:nrow(nas)) {
    if (patch_method == 'linear') {
      # failsafe against NAs at the start when using linear
      if (nas$start[j] == 1) {
        if (verbose) {
          message("NAs found at the start of column ", col, ".",
                  " Using method = 'after' for this instance.")
        }
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$stop[j] + 1, col]
      }
      # failsafe against NAs at the end when using linear
      if (nas$stop[j] == nrow(input)) {
        if (verbose) {
          message("NAs found at the end of column ", col, ".",
                  " Using method = 'before' for this instance.")
        }
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$start[j] - 1, col]
      }
      if (nas$start[j] != 1 && nas$stop[j] != nrow(input)) {
        missing_interval <- (nas$start[j] - 1):(nas$stop[j] + 1)
        replacement <- seq(from = input[nas$start[j] - 1, col],
                           to = input[nas$stop[j] + 1, col],
                           length.out = nas$stop[j] - nas$start[j] + 3
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
      if (nas$start[j] == 1) {
        if (verbose) {
          message("NAs found at the start of of column ", col, ".",
            " Using method = 'after' for this instance.")
        }
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$stop[j] + 1, col]
      }
      else {
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$start[j] - 1, col]
      }
    }
    if (patch_method == 'after') {
      if (nas$stop[j] == nrow(input)) {
        if (verbose) {
          message("NAs found at the end of of column ", col, ".",
            " Using method = 'before' for this instance.")
        }
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$start[j] - 1, col]
      }
      else {
        missing_interval <- nas$start[j]:nas$stop[j]
        replacement <- input[nas$stop[j] + 1, col]
      }        
    }
    input[missing_interval, col] <- replacement
  }
  return(input)
}

#' Discard readings
#'
#' Discard the data from one or more phases of one or more probes.
#'
#' @param input A computer-friendly data frame.
#'   The output of \code{\link{melt_resp}} or any downstream function.
#' @param probes The probe(s) from which to discard data. Ommit to discard phases
#'   from all probes.
#' @param cycles The cycles(s) to discard.
#'
#' @return the input data frame without the discarded readings.
#'
#' @export
#'
discard_cycles <- function(input, probes, cycles) {
  cycles <- check_arg_in_data(cycles, input$trimmed$cycle, "cycles")
  target_cycles <- input$trimmed$cycle %in% cycles

  if (all(!target_cycles)) {
    stop("Couldn't find specified cycles in the input.")
  }

  if (missing(probes)) {
    target_probes <- rep(TRUE, nrow(input$trimmed))
  } else {
    probes <- check_arg_in_data(probes, input$trimmed$probe, "probes")
    target_probes <- input$trimmed$probe %in% probes
    if (all(!target_probes)) {
      stop("Couldn't find specified probes in the input.")
    }
  }

  to_keep <- !(target_cycles & target_probes)

  if (all(to_keep)) {
    stop("Couldn't find specified probe-cycle combination in the input.")
  }

  input$trimmed <- input$trimmed[to_keep, ]
  return(input)
}
