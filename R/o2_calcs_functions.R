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
#' 
#' @return A data frame containing measurements valid for further analyses
#' 
#' @export
#' 
clean_meas <- function(input, wait = 0, cycle_max = Inf, auto_cut_last = FALSE){
  if (wait > cycle_max) {
    stop("wait must not be bigger than cycle_max.")
  }
  
  input$date <- as.Date(input$date_time)
  input$real_time <- chron::times(strftime(input$date_time, "%H:%M:%S"))

  # Removing Non-Measurement Data  
  input <- input[grepl("^M", input$phase), ]
  
  # this is used just to ensure that the phases maintain their order,
  # even it they don't start at one or are not sorted at the start.
  phase_order <- as.numeric(gsub("[M]", "", unique(input$phase)))
  phase_order <- order(phase_order)
  
  input$phase <- factor(input$phase, levels = unique(input$phase)[phase_order])


  #the rest has to be done on a probe by_probe basis.

  by_probe <- split(input, input$probe)

  recipient <- lapply(names(by_probe), function(the_probe) {

    trimmed_db <- by_probe[[the_probe]]

    # Removing the final measurement phase if necessary or forced (tail error)
    rows_per_phase <- table(trimmed_db$phase)

    if (length(rows_per_phase) > 1) {
      mean_rows_per_phase <- mean(rows_per_phase[-length(rows_per_phase)])
    } else {
      mean_rows_per_phase <- rows_per_phase
    }

    if (tail(rows_per_phase, 1) < wait | auto_cut_last) {
      all_but_last <- trimmed_db$phase != tail(levels(trimmed_db$phase), 1)
      trimmed_db <- trimmed_db[all_but_last, ]
      trimmed_db$phase <- droplevels(trimmed_db$phase)
    }


    # cut off first n rows from 'M' phases
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

      # and now, by phase, calculate the passing time and 
      # store the start and end points of that measurement
      aux <- split(trimmed_db, trimmed_db$phase)

      aux <- aux[sapply(aux, nrow) > 0]

      aux <- lapply(aux, function(x) {
        x$Start.Meas <- x$Real.Time[1]
        x$End.Meas <- x$Real.Time[nrow(x)]
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
    stop("No measurement data found. Is wait too long? ",
         "Is the phase data correct?")
  }
  
  units(wait) <- "seconds"
  attributes(output)$wait_time <- wait

  return(output)
}

#' Calculate oxygen delta for each cycle
#' 
#' Subtracts the initial oxygen value to the remaining, 
#' for each probe*cycle combination.
#' 
#' @param input a data frame with cleaned measurements.
#'  The output of \code{\link{clean_meas}}.
#' 
#' @return The input data frame with calculated deltas for o2 and airsat.
#' 
#' @export
#' 
calc_delta <- function(input) {

  # remove any existing delta column to avoid conflicts
  input$o2_delta <- NULL
  
  by_probe <- split(input, input$probe)
  
  recipient <- lapply(names(by_probe), function(the_probe) {

    trimmed_db <- by_probe[[the_probe]]

    by_phase <- split(trimmed_db, trimmed_db$phase)
    
    recipient <- lapply(by_phase, function(the_phase) {
      the_phase$o2_delta <- the_phase$o2 - the_phase$o2[1]
      the_phase$airsat_delta <- the_phase$airsat - the_phase$airsat[1]
      return(the_phase)
    })
    
    output <- data.table::rbindlist(recipient)

  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  output <- transfer_attributes(input, output)

  return(output)
}
