#' Prepare imported measurements for further analyses
#'
#' @param input the data frame containing imported oxygen measurements
#' @param wait integer: the number of first rows for each measurement phase (M)
#'  which should be reassigned to the wait phase (W). Note: If your
#'  phase-tracking device already assigns a wait phase, set this to 0.
#' @param cycle_max integer: Max allowed length (in number of rows) that each
#'  cycle is allowed to have. Defaults to Inf (i.e. the whole cycle is used for
#'  metabolic rate estimations).
#' @param auto_cut_last Should the last recorded phase be automatically cut?
#'  useful if the experiment was terminated mid-measurement.
#' @param first_cycle either a vector of length one with the first cycle for all
#'  probes, or a named vector containing the first cycle for each individual
#'  probe (with the respective probe name).
#'
#' @return A data frame containing measurements valid for further analyses
#'
#' @export
#'
trim_resp <- function(input, wait = 0, cycle_max = Inf,
                      auto_cut_last = FALSE, first_cycle = 1){

  if (is.null(input$melted)) {
    stop("Couldn't find object 'melted' inside input.",
         " Have you run melt_resp?", call. = FALSE)
  }
  if (wait > cycle_max) {
    stop("argument 'wait' must not be bigger than argument 'cycle_max'.")
  }
  melted <- input$melted

  if (!is.null(input$probe_info)) {
    target_probes <- input$probe_info$probe
  } else {
    target_probes <- unique(melted$probe)
  }

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

  melted$date <- as.Date(melted$date_time)
  melted$real_time <- chron::times(strftime(melted$date_time, "%H:%M:%S"))

  # Removing Non-Measurement Data
  melted <- melted[grepl("^M", melted$phase), ]

  # this is used just to ensure that the phases maintain their order,
  # even it they don't start at one or are not sorted at the start.
  phase_order <- as.numeric(gsub("[M]", "", unique(melted$phase)))
  phase_order <- order(phase_order)

  melted$phase <- factor(melted$phase, levels = unique(melted$phase)[phase_order])

  # the rest has to be done on a probe by probe basis.
  by_probe <- split(melted, melted$probe)

  recipient <- lapply(target_probes, function(the_probe) {

    trimmed_db <- by_probe[[the_probe]]

    # trim away the cycles that happen before the first_cycle
    trimmed_db <- trimmed_db[trimmed_db$cycle >= first_cycle[the_probe], ]

    # Remove the final measurement phase if necessary (tail error)
    if (nrow(trimmed_db) > 0) {
      rows_per_phase <- table(trimmed_db$phase)
      if (tail(rows_per_phase, 1) < wait | auto_cut_last) {
        all_but_last <- trimmed_db$phase != tail(levels(trimmed_db$phase), 1)
        trimmed_db <- trimmed_db[all_but_last, ]
        trimmed_db$phase <- droplevels(trimmed_db$phase)
      }
    }

    # remove the wait phases
    if (nrow(trimmed_db) > 0 && wait != 0) {
      # the code below grabs the 1:nrow vector, breaks it out by phase,
      # and then selects only the rows that fall within the desired
      # measurement periods.
      index <- unlist(tapply(X = 1:nrow(trimmed_db),
                             INDEX = trimmed_db$phase,
                             FUN = function(x) {
                                x[(wait + 1):min(length(x), cycle_max)]
                             }),
                      use.names = FALSE)
      trimmed_db <- trimmed_db[index, ]
    }

    # reset rownames
    if (nrow(trimmed_db) > 0) {
      row.names(trimmed_db) <- 1:nrow(trimmed_db)

      # and now, by phase, calculate the passing time
      aux <- split(trimmed_db, trimmed_db$phase)

      aux <- aux[sapply(aux, nrow) > 0]

      aux <- lapply(aux, function(x) {
        x$phase_time <- as.numeric(difftime(
                          time1 = x$date_time,
                          time2 = x$date_time[1],
                          units = 's'
                        ))
        units(x$phase_time) <- "s"

        return(x)
      })
      trimmed_db <- as.data.frame(data.table::rbindlist(aux))
    }

    return(trimmed_db)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  if (nrow(output) == 0) {
    stop("No measurement data found. Is wait too long?",
         " Is the phase data correct?", call. = FALSE)
  }

  attributes(output)$wait <- paste(wait, "data points")

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
      the_phase$o2_delta <- the_phase$o2 - mean(the_phase$o2[1:zero_buffer],
                                                na.rm = TRUE)
      the_phase$airsat_delta <- the_phase$airsat - mean(the_phase$airsat[1:zero_buffer],
                                                        na.rm = TRUE)
      return(the_phase)
    })

    output <- data.table::rbindlist(recipient)

  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  output <- transfer_attributes(input, output)

  return(output)
}
