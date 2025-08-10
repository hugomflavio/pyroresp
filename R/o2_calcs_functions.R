#' Prepare imported measurements for further analyses
#'
#' @param input the data frame containing imported oxygen measurements
#' @param meas_min integer: Minimum number of seconds a measurement phase must
#'   have to be considered valid. Defaults to 10.
#' @param meas_max integer: Max number of seconds to keep from each
#'  measurement phase. Defaults to Inf (i.e. the whole phase is used for
#'  metabolic rate estimations).
#' @param cut_last Should the last recorded phase be automatically cut?
#'  useful if the experiment was terminated mid-measurement.
#' @param first_cycle either a vector of length one with the first cycle for all
#'  probes, or a named vector containing the first cycle for each individual
#'  probe (with the respective probe name).
#'
#' @return A data frame containing measurements valid for further analyses
#'
#' @export
#'
trim_resp <- function(input, meas_max = Inf, meas_min = 60,
                      cut_last = FALSE, first_cycle = 1){

  if (is.null(input$melted)) {
    stop("Couldn't find object 'melted' inside input.",
         " Have you run melt_resp?", call. = FALSE)
  }

  units(meas_max) <- "s"
  units(meas_min) <- "s"

  # simplify object name
  melted <- input$melted

  # keep only data for probes for which we have info, if possible
  if (!is.null(input$probe_info)) {
    target_probes <- input$probe_info$probe
  } else {
    target_probes <- unique(melted$probe)
  }

  # expand firsy_cycle argument as needed
  if (length(first_cycle) == 1) {
    first_cycle <- rep(first_cycle, length(target_probes))
    names(first_cycle) <- target_probes
  } else {
    link <- !(target_probes %in% names(first_cycle))
    if (any(link)) {
      stop("Length of argument 'first_cycle' is not 1 but named values",
           " not supplied for all probes (missing",
           paste0(target_probes[link], collapse = ", "), ").",
           call. = FALSE)
    }
  }

  # split date and time
  melted$date <- as.Date(melted$date_time)
  melted$real_time <- chron::times(strftime(melted$date_time, "%H:%M:%S"))

  # keep only M data
  melted <- melted[grepl("^M", melted$phase), ]

  # this is used just to ensure that the phases maintain their order,
  # even it they don't start at one or are not sorted at the start.
  phase_order <- as.numeric(gsub("M", "", unique(melted$phase)))
  phase_order <- order(phase_order)

  melted$phase <- factor(melted$phase, levels = unique(melted$phase)[phase_order])

  # the rest has to be done on a probe by probe basis.
  by_probe <- split(melted, melted$probe)

  recipient <- lapply(target_probes, function(the_probe) {

    # simplify object name
    trimmed_db <- by_probe[[the_probe]]

    # trim away the cycles that happen before the first_cycle
    trimmed_db <- trimmed_db[trimmed_db$cycle >= first_cycle[the_probe], ]

    # Remove the final measurement phase if desired
    if (cut_last) {
      all_but_last <- trimmed_db$phase != tail(unique(trimmed_db$phase), 1)
      trimmed_db <- trimmed_db[all_but_last, ]
      trimmed_db$phase <- droplevels(trimmed_db$phase)
    }

    # remove data beyond maximum phase duration
    if (nrow(trimmed_db) > 0) {
      keep <- trimmed_db$phase_time <= meas_max
      trimmed_db <- trimmed_db[keep, ]
    }

    # remove cycles that are too short
    if (nrow(trimmed_db) > 0) {
      check <- aggregate(trimmed_db$phase_time, list(trimmed_db$phase), max)
      colnames(check) <- c("phase", "time")
      check$check <- check$time >= meas_min
      keep <- trimmed_db$phase %in% check$phase[check$check]
      trimmed_db <- trimmed_db[keep, ]
    } 

    if (nrow(trimmed_db) > 0) {
      row.names(trimmed_db) <- 1:nrow(trimmed_db)
    }

    return(trimmed_db)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  if (nrow(output) == 0) {
    stop("No valid measurement data found. Is wait too long?",
         " Is the phase information correct?", call. = FALSE)
  }

  input$trimmed <- output

  return(input)
}

#' Calculate oxygen delta for each cycle
#'
#' Subtracts the initial oxygen value to the remaining,
#' for each probe*cycle combination.
#'
#' @param input a data frame with trimmed measurements.
#'  The output of \code{\link{trim_resp}}.
#' @param zero_buffer when calculating the delta, a value must be assigned to 0.
#'  Traditionally, this is the very first value of the cycle. However, due to
#'  natural probe noise, this can cause an upward or downard shift to the whole
#'  delta line. This can be countered by calculating the mean of a few initial
#'  values, instead of relying only on the very first. zero_buffer sets how many
#'  values should be used to estimate the starting O2 concentration for the
#'  cycle.
#'
#' @return The input data frame with calculated deltas for o2 and airsat.
#'
#' @export
#'
calc_delta <- function(input, zero_buffer = 3) {

  # remove any existing delta column to avoid conflicts
  input$o2_delta <- NULL

  by_probe <- split(input, input$probe)

  recipient <- lapply(names(by_probe), function(the_probe) {

    trimmed_db <- by_probe[[the_probe]]

    by_phase <- split(trimmed_db, trimmed_db$phase)

    recipient <- lapply(by_phase, function(the_phase) {
      if (all(is.na(the_phase$o2))) {
        # if everything is NA, then we're only going to get NAs
        # in the end, so there's no point in trimming.
        aux <- the_phase
      } else {
        # but if there's data among the NAs, we don't want
        # NAs at the start to mess up our starting points.
        aux <- the_phase[!is.na(the_phase$o2), ]
      }
      first_o2 <- mean(aux$o2[1:zero_buffer])
      first_airsat <- mean(aux$airsat[1:zero_buffer])
      the_phase$o2_delta <- the_phase$o2 - first_o2
      the_phase$airsat_delta <- the_phase$airsat - first_airsat
      return(the_phase)
    })

    output <- data.table::rbindlist(recipient)

  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  output <- transfer_attributes(input, output)

  return(output)
}
