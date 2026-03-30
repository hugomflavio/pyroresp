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
#'     - animal_id    : The ID of the animal
#'     - animal_mass  : The mass of the animal, in grams
#'     - animal_vol   : Optional: The volume of the animal.
#'                      Assumed to be equal to mass if missing.
#'     - chamber_vol  : The non-corrected volume of the chamber
#'     - probe        : The device-channel combination for the probe
#'     - first_cycle  : The first cycle of valid data for that animal
#'
#' @return A list containing a phases dataframe and a pyro list with the
#'   individual source data frames (in source_data), as well as a single,
#'   combined data frame organized by time (in compiled_data).
#'
#' @export
#'
load_experiment <- function(folder, date_format, tz = Sys.timezone(),
    phases, probe_info, encoding = "ISO-8859-1", mass_unit = "g",
    vol_unit = "ml") {

  if (length(folder) == 0 || !dir.exists(folder)) {
    stop('Could not find target folder')
  }

  if (length(folder) > 1) {
    stop('"folder" should be a string of length 1.')    
  }

  if (!missing(probe_info)) {
    probe_info <- process_probe_info(probe_info,
                                     vol_unit = vol_unit,
                                     mass_unit = mass_unit)
  } else {
    probe_info <- NULL
  }

  if (any(sapply(names(phases), nchar) > 4)) {
    warning("Long device names detected in the phases input.",
            " Confirm these are correct: ", 
            paste(names(phases), collapse = ", "), ".")
  }

  pyro <- load_pyro_folder(folder, date_format = date_format,
                           tz = tz, encoding = encoding)

  output <- list(phases = phases,
                 pyro = pyro,
                 probe_info = probe_info)

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
load_pyro_folder <- function(folder, date_format, tz, 
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
  names(source_data) <- files

  compiled_data <- compile_sources(source_data)

  output <- list(source_data = source_data,
                 compiled_data = compiled_data)
  return(output)
}

#' merge a list of data sources into a single table using date_time
#' as the link. Ensures there is one observation per second.
#' 
#' @param input a list of equally formatted tables
#' 
#' @return A table with the compiled data
#' 
#' @export
#' 
compile_sources <- function(input) {
  recipient <- input[[1]]
  head(recipient)
  recipient <- recipient[!duplicated(recipient$date_time), ]
  if (length(input) > 1) {
    for (i in input[-1]) {
      i <- i[!duplicated(i$date_time), ]
      recipient <- merge(recipient, i, by = 'date_time', all = TRUE)
    }
  }

  attributes(recipient)$latest_batch_start <- 1
  return(recipient)
}

#' Wrapper to get experiment data ready for further analyses
#'
#' Perform standard processing operations to the pyro/phases files.
#' 
#' @inheritParams trim_resp
#' @inheritParams assign_phases
#' @inheritParams patch_NAs
#' @inheritParams calc_delta
#' @param input The output of \code{\link{load_experiment}}
#' @param original_o2 The o2 unit the data was captured in.
#' @param convert_o2_to The o2 unit desired for the final results.
#' @param conv_with_mean_temp Calculate the mean temperature experienced
#'   during a given phase and use that for O2 unit conversion.
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
#'   of these arguments at the same time. Input must be numeric.
#' @param verbose Logical. Should steps being taken be detailed with messages.
#'   Defaults to TRUE.
#' 
#' @return An updated experiment list, containing a cleaned object with the
#'   processed data.
#' 
#' @export 
#'
process_experiment <- function(input, wait = 0, tail_trim = 0,
    meas_max = Inf, meas_min = 60, original_o2, convert_o2_to,
    conv_with_mean_temp = FALSE,
    na_action = c("ignore", "linear", "before", "after", "remove"),
    zero_buffer = 3, first_cycle = 1,
    min_temp, max_temp, start_time, stop_time, from_cycle, to_cycle,
    verbose = TRUE) {

  na_action <- match.arg(na_action)

  all_units <- c("hPa", "kPa", "torr", "mmHg", "inHg", "mg_per_l", 
           "ug_per_l", "umol_per_l", "mmol_per_l", "ml_per_l",
           "mg_per_kg", "ug_per_kg", "umol_per_kg", "mmol_per_kg", 
           "ml_per_kg")                                                            

  if (!missing(convert_o2_to) && !(convert_o2_to %in% all_units)) {
    stop("the 'convert_o2_to' argument is not an acceptable unit. ",
         "Please choose one of the following: ", 
         paste(all_units, collapse = ", "))
  }

  if (!missing(min_temp) & !missing(max_temp)) {
    stop("Please use only one of 'min_temp' or 'max_temp'",
         " at a time. See function help for details.")
  }

  if (!missing(start_time)) {
    if (verbose) {
      message(paste0("M: Discarding readings before ", start_time, "."))
    }
    keep <- input$pyro$compiled_data$date_time >= as.POSIXct(start_time)
    if (all(!keep)) {
      stop ("Data ends before ", start_time, ".")
    } else {
      input$pyro$compiled_data <- input$pyro$compiled_data[keep, ]
    }
  }

  if (!missing(stop_time)) {
    if (verbose) {
      message(paste0("M: Discarding readings after ", stop_time, "."))
    }
    keep <- input$pyro$compiled_data$date_time <= as.POSIXct(stop_time)
    if (all(!keep)) {
      stop ("Data starts after ", stop_time, ".")
    } else {
      input$pyro$compiled_data <- input$pyro$compiled_data[keep, ]
    }
  }

  if (verbose) {
    message("M: Merging pyroscience and phases file.")
  }
  input <- assign_phases(input, wait = wait, tail_trim = tail_trim)

  if (na_action %in% c("linear", "before", "after")) {
      if (verbose) {
        message("M: Addressing NA's in the data (na_action = \"",
                na_action, "\".")
      }
      input$phased <- patch_NAs(input$phased, na_action = na_action, 
                                verbose = FALSE)
  }
  
  if (verbose) {
    message("M: Melting resp data into computer-friendly format")
  }
  input <- melt_resp(input = input)

  if (na_action == "remove") {
    if (verbose) {
      message("M: Removing rows with no oxygen data (na_action = \"remove\").")
    }
    input$melted <- input$melted[!is.na(input$melted$o2), ]
  }
  if (verbose) message("M: Removing unwanted data (flush, wait, etc).")
    input <- trim_resp(input = input,
                       meas_max = meas_max,
                       meas_min = meas_min,
                       first_cycle = first_cycle)

  if (conv_with_mean_temp) {
    if (verbose) message("M: Calculating mean temp per phase.")
    aux <- aggregate(as.numeric(input$trimmed$temp),
                     list(probe = input$trimmed$probe,
                          phase = input$trimmed$phase),
                     mean, na.rm = TRUE)
    input_probe_phase <- paste(input$trimmed$probe, input$trimmed$phase)
    aux_probe_phase <- paste(aux$probe, aux$phase)
    link <- match(input_probe_phase, aux_probe_phase)
    input$trimmed$conv_temp <- aux$x[link]

    aux <- aggregate(as.numeric(input$melted$temp),
                     list(probe = input$melted$probe,
                          phase = input$melted$phase),
                     mean, na.rm = TRUE)
    input_probe_phase <- paste(input$melted$probe, input$melted$phase)
    aux_probe_phase <- paste(aux$probe, aux$phase)
    link <- match(input_probe_phase, aux_probe_phase)
    input$melted$conv_temp <- aux$x[link]

  } else {
    input$trimmed$conv_temp <- as.numeric(input$trimmed$temp)
    input$melted$conv_temp <- as.numeric(input$melted$temp)
  }

  if (verbose) message("M: Calculating air saturation.")

  o2_conv_cols <- c("o2", "temp", "sal", "pressure")
  not_NA <- complete.cases(input$trimmed[, o2_conv_cols])
  
  if (missing(original_o2)) {
    original_o2 <- sub("/", "_per_", units(input$trimmed$o2))
    original_o2 <- gsub("L", "l", original_o2)
  }

  input$trimmed$airsat <- NA
  input$trimmed$airsat[not_NA] <- 
    respirometry::conv_o2(
      o2 = as.numeric(input$trimmed$o2[not_NA]),
      from = original_o2,
      to = "percent_a.s.", 
      temp = as.numeric(input$trimmed$conv_temp[not_NA]), 
      sal = as.numeric(input$trimmed$sal[not_NA]), 
      atm_pres = as.numeric(input$trimmed$pressure[not_NA])
    )
  units(input$trimmed$airsat) <- "percent"

  if (!missing(convert_o2_to)) {
    
    if (verbose) {
      message("M: Converting oxygen unit from ", original_o2, 
          " to ", convert_o2_to, ".")  
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
        to = convert_o2_to, 
        temp = as.numeric(input$trimmed$conv_temp[not_NA]), 
        sal = as.numeric(input$trimmed$sal[not_NA]), 
        atm_pres = as.numeric(input$trimmed$pressure[not_NA])
      )
    units(input$trimmed$o2) <- gsub("_per_", "/", convert_o2_to)

    # convert from melted too, which is used for plot_meas
    not_NA <- complete.cases(input$melted[, o2_conv_cols])
    input$melted$o2[!not_NA] <- NA
    input$melted$o2 <- as.numeric(input$melted$o2)
    input$melted$o2[not_NA] <-
      respirometry::conv_o2(
        o2 = input$melted$o2[not_NA],
        from = original_o2,
        to = convert_o2_to, 
        temp = as.numeric(input$melted$conv_temp[not_NA]), 
        sal = as.numeric(input$melted$sal[not_NA]), 
        atm_pres = as.numeric(input$melted$pressure[not_NA])
      )
    units(input$melted$o2) <- gsub("_per_", "/", convert_o2_to)
  }

  # remove temporary temperature columns
  input$trimmed$conv_temp <- NULL
  input$melted$conv_temp <- NULL

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
#' @inheritParams subtract_bg
#' @inheritParams filter_r2
#' @inheritParams calc_slopes
#' @param input The output of \code{\link{process_experiment}}
#' 
#' @return An updated experiment list, containing two new objects; one with
#'   all the slopes, and another with the slopes that pass the r2 threshold.
#' 
#' @export 
#'
process_slopes <- function(input, r2 = 0.95, pre, post, method,
                           correct_for_water_vol = FALSE) {
    input <- calc_slopes(input, correct_for_water_vol = correct_for_water_vol)
  
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

  if (is.null(input$good_slopes)) {
    stop("input has no good slope data.",
         call. = FALSE)
  }
  input$mr <- calc_mr(input$good_slopes)

  # convert seconds to hours (more common)
  the_seconds <- which(units(input$mr$mr_abs)$denominator == "s")
  if (length(the_seconds) > 0) {
    units(input$mr$mr_abs)$denominator[the_seconds] <- "h"
  }
  if ("mr_g" %in% colnames(input$mr)) {
    the_seconds <- which(units(input$mr$mr_g)$denominator == "s")
    if (length(the_seconds) > 0) {
      units(input$mr$mr_g)$denominator[the_seconds] <- "h"
    }
  }

  if (is.null(G) && is.null(q) && is.null(p) && is.null(n)) {
    input$smr <- NULL
  } else {
    input$smr <- calc_smr(input$mr, G = G, q = q, p = p, n = n)
    keep_these <- !(colnames(input$smr) %in% c("id", "animal_mass", "water_vol"))
    smr_aux <- input$smr[, keep_these]
    input$smr <- merge(input$probe_info, smr_aux, 
                        by = "probe", all = TRUE)

    smr_cols <- colnames(input$smr)[grepl("_mr", colnames(input$smr))]
    for (i in smr_cols) {
      prefix <- sub("_mr", "", i)
      if (!is.null(input$bg$pre)) {
        new_col <- paste0(prefix, "pre_bg_pct")
        aux <- -input$smr[, i]
        if ("animal_mass" %in% colnames(input$smr)) {
          aux <- aux * input$smr$animal_mass
        }
        aux <- aux  / input$smr$water_vol
        units(aux) <- units(input$bg$pre$bg$slope)
        bg_link <- match(input$smr$probe, input$bg$pre$bg$probe)
        input$smr[, new_col] <- input$bg$pre$bg$slope[bg_link] / aux
        # units is now "1"; changing to percent automatically multiplies by 100
        units(input$smr[, new_col]) <- "percent"
      }
      if (!is.null(input$bg$post)) {
        new_col <- paste0(prefix, "post_bg_pct")
        aux <- -input$smr[, i]
        if ("animal_mass" %in% colnames(input$smr)) {
          aux <- aux * input$smr$animal_mass
        }
        aux <- aux  / input$smr$water_vol
        units(aux) <- units(input$bg$post$bg$slope)
        bg_link <- match(input$smr$probe, input$bg$post$bg$probe)
        input$smr[, new_col] <- input$bg$post$bg$slope[bg_link] / aux
        # units is now "1"; changing to percent automatically multiplies by 100
        units(input$smr[, new_col]) <- "percent"
      }
      # how to make a bg summary when using a reference chamber?
    }
  }
  
  input$mmr <- extract_mmr(input$mr)
  mmr_aux <- input$mmr[, !(colnames(input$mmr) %in% c("id", "mass", "volume"))]
  input$mmr <- merge(input$probe_info, mmr_aux,
                     by = "probe", all = TRUE)

  ## calculating bg as a percentage of mmr seems irrelevant.

  # if (!is.null(input$smr) && !is.null(input$bg$pre)) {
  #   aux <- -input$mmr$mr_g * input$mmr$mass / input$smr$volume
  #   units(aux) <- units(input$bg$pre$bg$slope)
  #   bg_link <- match(input$mmr$probe, input$bg$pre$bg$probe)
  #   input$mmr$pre_bg_pct <- input$bg$pre$bg$slope[bg_link] / aux
  #   # units is now "1"; changing to percent automatically multiplies by 100
  #   units(input$mmr$pre_bg_pct) <- "percent"
  # }
  # if (!is.null(input$smr) && !is.null(input$bg$post)) {
  #   aux <- -input$mmr$mr_g * input$mmr$mass / input$smr$volume
  #   units(aux) <- units(input$bg$post$bg$slope)
  #   bg_link <- match(input$mmr$probe, input$bg$post$bg$probe)
  #   input$mmr$post_bg_pct <- input$bg$post$bg$slope[bg_link] / aux
  #   # units is now "1"; changing to percent automatically multiplies by 100
  #   units(input$mmr$post_bg_pct) <- "percent"
  # }
  
  return(input)
}
 