#' Wrapper to load the data for an experiment run (pyro and phases files).
#'
#' Scans the target folder for a phases file and the raw pyroscience
#' experiment files. Imports them.
#'
#' @param folder the path to the folder where the experiment data is stored
#' @inheritParams read_pyro_raw_file
#' @param phases a phases object. The output of \code{\link{process_phases}} 
#' @inheritParams load_phases
#' @param probe_info a dataframe containing animal information. This
#'   dataframe must contain the following columns:
#'     - id           : The ID of the animal
#'     - mass         : The mass of the animal, in grams
#'     - volume       : The non-corrected volume of the chamber + tubing
#'     - probe        : The device-channel combination for the probe
#'     - first_cycle  : The first cycle of valid data for that animal
#'
#' @return A list containing a phases dataframe and a pyro list with the
#'   individual source data frames (in source_data), as well as a single,
#'  combined data frame organized by time (in compiled_data).
#'
#' @export
#'
load_experiment <- function(folder, date_format, tz = Sys.timezone(),
    phases, probe_info, encoding = "ISO-8859-1") {

  if (length(folder) == 0 || !dir.exists(folder)) {
    stop('Could not find target folder')
  }

  if (length(folder) > 1) {
    stop('"folder" should be a string of length 1.')    
  }

  if (!missing(probe_info)) {
    required_cols <- c("id", "mass", "volume", "probe")
    cols_missing <- !(required_cols %in% colnames(probe_info))
    if (any(cols_missing)) {
      stop("The following required columns are missing ",
         "from the probe_info input: ",
         paste0(required_cols[cols_missing], collapse = ", "),
         call. = FALSE)
    }
  }

  if (any(sapply(names(phases), nchar) > 4)) {
    warning("Long device names detected in the phases input. Are you sure",
            " you appended the device names correctly to the file name?",
            " These are the current device names: ", 
            paste(names(phases), collapse = ", "), ".")
  }

  pyro <- load_pyro_data(folder, date_format = date_format, tz = tz,
               encoding = encoding)

  output <- list(phases = phases, pyro = pyro)
  
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
#'   containing a "ChannelData" folder inside
#' @inheritParams read_pyro_raw_file
#' @param type One of "Oxygen" to read only oxygen files, "pH" to read only pH
#'   files, or "Oxygen|pH" to read both.
#'
#' @export
#'
load_pyro_data <- function(folder, date_format, tz, 
    type = c("Oxygen", "pH", "Oxygen|pH"), encoding = "ISO-8859-1") {
  type <- match.arg(type)

  files <- list.files(paste0(folder, '/ChannelData/'))

  file_link <- grepl(type, files)

  files <- files[file_link]

  source_data <- lapply(files, function(i) {
    read_pyro_raw_file(paste0(folder, '/ChannelData/', i),
               date_format = date_format, tz = tz,
               encoding = encoding)
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
#' @inheritParams trim_resp
#' @param convert_o2_unit_to The o2 unit desired for the final results
#' @param patch_NAs Logical. Should NA values found in the raw data be patched?
#'   Defaults to TRUE.
#' @inheritParams patch_NAs
#' @inheritParams calc_delta
#' @param min_temp,max_temp 
#'   For temperature ramp experiments. The minimum OR maximum temperatures
#'   that must be reached before data is considered valid. Discards all phases
#'   prior to this temperature being reached. Use only one of the two arguments
#'   at a time.
#' @param start_time,stop_time
#'   Trim the experiment to a specific time period. You may use one or both of
#'   these arguments at the same time. Input must be a string in 
#'   YYYY-MM-DD HH:MM:SS format.
#' @param from_cycle,to_cycle
#'   Trim the experiment to a specific group of cycles. You may use one or both
#'  of these arguments at the same time. Input must be numeric.
#' @param verbose Logical. Should steps being taken be detailed with messages.
#'   Defaults to TRUE.
#' 
#' @return An updated experiment list, containing a cleaned object with the
#'   processed data.
#' 
#' @export 
#'
process_experiment <- function(input, wait, cycle_max = Inf, convert_o2_unit_to,
    patch_NAs = TRUE, patch_method = c("linear", "before", "after"),
    zero_buffer = 3, first_cycle = 1,
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

  if (verbose) {
    message("M: Merging pyroscience and phases file.")
  }
  input <- assign_phases(input)

  if (patch_NAs) {
      if (verbose) {
        message("M: Patching NA's in the data.")
      }
      input$phased <- patch_NAs(input$phased, patch_method = patch_method, 
                                verbose = FALSE)
  }
  
  if (verbose) message("M: Melting resp data into computer-friendly format")
    input <- melt_resp(input = input)

  if (verbose) message("M: Removing unwanted data.")
    input <- trim_resp(input = input, wait = wait,
                       cycle_max = cycle_max,
                       first_cycle = first_cycle)

  if (verbose) message("M: Calculating air saturation")

  o2_conv_cols <- c("o2", "temp", "sal", "pressure")
  not_NA <- complete.cases(input$trimmed[, o2_conv_cols])
  original_o2 <- sub("/", "_per_", units(input$trimmed$o2))

  input$trimmed$airsat <- NA
  input$trimmed$airsat[not_NA] <- 
    respirometry::conv_o2(
      o2 = as.numeric(input$trimmed$o2[not_NA]),
      from = original_o2,
      to = "percent_a.s.", 
      temp = as.numeric(input$trimmed$temp[not_NA]), 
      sal = as.numeric(input$trimmed$sal[not_NA]), 
      atm_pres = as.numeric(input$trimmed$pressure[not_NA])
    )
  units(input$trimmed$airsat) <- "percent"

  if (!missing(convert_o2_unit_to)) {
    
    if (verbose) {
      message("M: Converting oxygen unit from ", original_o2, 
          " to ", convert_o2_unit_to, ".")  
    }
    
    # if there is an o2 value but not all the others
    if (any(!is.na(input$trimmed$o2) & !not_NA)) {
      warning_cases <- sum(!is.na(input$trimmed$o2) & !not_NA)
      warning("Invalidating ", warning_cases, "oxygen value(s) as one ",
        "or more of the respective temperature, salinity, or pressure ",
        "values are missing, making it impossible to convert unit.")
      input$trimmed$o2[!not_NA] <- NA
    }

    input$trimmed$o2 <- as.numeric(input$trimmed$o2)
    input$trimmed$o2[not_NA] <- 
      respirometry::conv_o2(
        o2 = input$trimmed$o2[not_NA],
        from = original_o2,
        to = convert_o2_unit_to, 
        temp = as.numeric(input$trimmed$temp[not_NA]), 
        sal = as.numeric(input$trimmed$sal[not_NA]), 
        atm_pres = as.numeric(input$trimmed$pressure[not_NA])
      )
    units(input$trimmed$o2) <- gsub("_per_", "/", convert_o2_unit_to)

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
  input$trimmed <- calc_delta(input$trimmed, zero_buffer = zero_buffer)

  cutoff <- NULL
  if (!missing(min_temp)) {
    units(min_temp) <- intToUtf8(c(176, 67))
    if (verbose) {
      message(paste0("M: Discarding phases under ", 
               min_temp, intToUtf8(c(176, 67))))
    }
    cutoff <- head(which(input$trimmed$temp > min_temp), 1)
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
    cutoff <- head(which(input$trimmed$temp < max_temp), 1)
    if (length(cutoff) == 0) {
      stop ("Temperature never dropped below ", max_temp, ".")
    }
  }
  if (!is.null(cutoff)) {
    the_matches <- which(input$trimmed$phase == input$trimmed$phase[cutoff])
    first_true <- head(the_matches, 1)
    time_break <- input$trimmed$date_time[first_true]
    input$trimmed <- input$trimmed[input$trimmed$date_time >= time_break, ]
  }

  if (!missing(start_time)) {
    if (verbose) {
      message(paste0("M: Discarding phases before ", start_time, "."))
    }
    the_matches <- which(input$trimmed$date_time >= as.POSIXct(start_time))
    cutoff <- head(the_matches, 1)
    if (length(cutoff) == 0) {
      stop ("Data ends before ", start_time, ".")
    } else {
      first_phase <- input$trimmed$phase[cutoff]
      the_matches <- which(input$trimmed$phase == first_phase)
      first_true <- head(the_matches, 1)
      break_ <- input$trimmed$date_time[first_true]
      input$trimmed <- input$trimmed[input$trimmed$date_time >= break_, ]
    }
  }

  if (!missing(stop_time)) {
    if (verbose) {
      message(paste0("M: Discarding phases after ", stop_time, "."))
    }
    the_matches <- which(input$trimmed$date_time <= as.POSIXct(stop_time))
    cutoff <- tail(the_matches, 1)
    if (length(cutoff) == 0) {
      stop ("Data starts after ", stop_time, ".")
    } else {
      last_phase <- input$trimmed$phase[cutoff]
      the_matches <- which(input$trimmed$phase == last_phase)
      last_true <- tail(the_matches, 1)
      break_ <- input$trimmed$date_time[last_true]
      input$trimmed <- input$trimmed[input$trimmed$date_time <= break_, ]
    }
  }

  if (!missing(from_cycle)) {
    if (verbose) {
      message(paste0("M: Discarding cycles prior to cycle ",
               from_cycle, "."))
    }
    input$trimmed <- input$trimmed[input$trimmed$cycle > from_cycle, ]
  }

  if (!missing(to_cycle)) {
    if (verbose) {
      message(paste0("M: Discarding cycles after cycle ", to_cycle, "."))
    }
    input$trimmed <- input$trimmed[input$trimmed$cycle < from_cycle, ]
  }

  return(input)
}

#' Wrapper to get calculate slopes, correct them, and filter them
#'
#' @param input The output of \code{\link{process_experiment}}
#' @inheritParams subtract_bg
#' @inheritParams filter_r2
#' 
#' @return An updated experiment list, containing two new objects; one with
#'   all the slopes, and another with the slopes that pass the r2 threshold.
#' 
#' @export 
#'
process_slopes <- function(input, r2 = 0.95, pre, post, method) {
    input <- calc_slopes(input)
  
    input <- subtract_bg(input = input, pre = pre,
                         post = post, method = method)

    input$good_slopes <- filter_r2(input$slopes, r2 = r2)
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
#'     Relevant details for the different methods are saved in the attributes.
#'  \item \code{mmr}: A data frame containing the cycle with the highest
#'     metabolic rate recorded.
#' }
#'
#' @export
#'
process_mr <- function(input, G = 1:4, 
             q = c(0.2, 0.25), p = 0.1, n = 10) {

  input$mr <- calc_mr(input$good_slopes)

  # convert seconds to hours (more common)
  the_seconds <- which(units(input$mr$mr_abs)$denominator == "s")
  if (length(the_seconds) > 0) {
    units(input$mr$mr_abs)$denominator[the_seconds] <- "h"
  }
  the_seconds <- which(units(input$mr$mr_g)$denominator == "s")
  if (length(the_seconds) > 0) {
    units(input$mr$mr_g)$denominator[the_seconds] <- "h"
  }
  
  if (is.null(G) && is.null(q) && is.null(p) && is.null(n)) {
    input$smr <- NULL
  } else {
    input$smr <- calc_smr(input$mr, G = G, q = q, p = p, n = n)
    keep_these <- !(colnames(input$smr) %in% c("id", "mass", "volume"))
    smr_aux <- input$smr[, keep_these]
    input$smr <- merge(input$probe_info, smr_aux, 
                        by = "probe", all = TRUE)

    smr_cols <- colnames(input$smr)[grepl("_mr_g", colnames(input$smr))]
    for (i in smr_cols) {
      prefix <- sub("_mr_g", "", i)
      if (!is.null(input$bg$pre)) {
        new_col <- paste0(prefix, "pre_bg_pct")
        aux <- input$smr[, i]
        aux <- -aux * input$smr$mass / input$smr$volume
        units(aux) <- units(input$bg$pre$bg$slope)
        bg_link <- match(input$smr$probe, input$bg$pre$bg$probe)
        input$smr[, new_col] <- input$bg$pre$bg$slope[bg_link] / aux
        # units is now "1"; changing to percent automatically multiplies by 100
        units(input$smr[, new_col]) <- "percent"
      }
      if (!is.null(input$bg$post)) {
        new_col <- paste0(prefix, "post_bg_pct")
        aux <- input$smr[, i]
        aux <- -aux * input$smr$mass / input$smr$volume
        units(aux) <- units(input$bg$post$bg$slope)
        bg_link <- match(input$smr$probe, input$bg$post$bg$probe)
        input$smr[, new_col] <- input$bg$post$bg$slope[bg_link] / aux
        # units is now "1"; changing to percent automatically multiplies by 100
        units(input$smr[, new_col]) <- "percent"
      }
    }
  }
  
  input$mmr <- extract_mmr(input$mr)
  mmr_aux <- input$mmr[, !(colnames(input$mmr) %in% c("id", "mass", "volume"))]
  input$mmr <- merge(input$probe_info, mmr_aux,
                      by = "probe", all = TRUE)

  if (!is.null(input$smr) && !is.null(input$bg$pre)) {
    aux <- -input$mmr$mr_g * input$mmr$mass / input$smr$volume
    units(aux) <- units(input$bg$pre$bg$slope)
    bg_link <- match(input$mmr$probe, input$bg$pre$bg$probe)
    input$mmr$pre_bg_pct <- input$bg$pre$bg$slope[bg_link] / aux
    # units is now "1"; changing to percent automatically multiplies by 100
    units(input$mmr$pre_bg_pct) <- "percent"
  }
  if (!is.null(input$smr) && !is.null(input$bg$post)) {
    aux <- -input$mmr$mr_g * input$mmr$mass / input$smr$volume
    units(aux) <- units(input$bg$post$bg$slope)
    bg_link <- match(input$mmr$probe, input$bg$post$bg$probe)
    input$mmr$post_bg_pct <- input$bg$post$bg$slope[bg_link] / aux
    # units is now "1"; changing to percent automatically multiplies by 100
    units(input$mmr$post_bg_pct) <- "percent"
  }

  return(input)
}
 