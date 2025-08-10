#' import phase times object from RData file created by Cricket
#' 
#' @param file The name of the file to be imported.
#' 
#' @return A data frame containing the phase information.
#' 
#' @export
#' 
load_phases <- function(file) {
  e <- new.env()
  load(file, envir = e)
  flush_data <- e$flush_data
  return(flush_data)
}

#' process the phase data
#' 
#' Calculates start and stop times for each phase
#' 
#' @param phases a phases data.frame. The output of \code{\link{load_phases}}
#' @param tz the timezone of the timestamp data. Defaults to the system timezone.
#' 
#' @return an updated phases object
#' 
#' @export
#' 
process_phases <- function(phases, tz = Sys.timezone()) {
  colnames(phases) <- tolower(colnames(phases))
  phases$timestamp <- as.POSIXct(phases$timestamp, tz = tz)
  phases <- split(phases, phases$chamber)
  phases <- lapply(phases, function(x) {
    x$start <- x$timestamp
    x$stop <- c(x$timestamp[-1], Inf)
    return(x)
  })
  return(phases)
}

#' Merge the pyro data with the phases.
#' 
#' Matches the two data sources by probe and timestamp and assigns the 
#' respective phases to each datapoint. The phases object must not have
#' gaps/overlaps/repeated phase names for downstream functions to work
#' properly.
#' 
#' @param input A list containing imported experiment data.
#'   The output from \code{\link{load_experiment}}
#' @param wait integer: the number of seconds at the start of each measurement 
#'   phase (M) that should be reassigned to the wait phase (W). Note: If your
#'   phase-tracking device already assigns a wait phase, set this to 0.
#' 
#' @return An updated experiment list containing a phased data frame
#' 
#' @export
#' 
assign_phases <- function(input, wait) {
  pyr <- input$pyro$compiled_data
  phases <- input$phases
  check_probes_match(colnames(pyr), phases)

  new_col_order <- 1
  tmp <- lapply(names(phases), function(prb) {
    # create column with placeholders if probe is in the data
    keyword <- paste0("_", prb, "$")
    if (any(grepl(keyword, colnames(pyr)))) {

      # initialize new columns
      new_phase_col <- paste0("phase_", prb)
      new_time_col <- paste0("phase_time_", prb)
      pyr[, new_phase_col] <- NA_character_
      pyr[, new_time_col] <- NA_real_

      # assign device-given phases
      for (i in 1:nrow(phases[[prb]])) {
        check1 <- pyr$date_time >= phases[[prb]]$start[i]
        check2 <- pyr$date_time < phases[[prb]]$stop[i]
        this_phase <- check1 & check2

        pyr[this_phase, new_phase_col] <- paste0(phases[[prb]]$phase[i],
                                                 phases[[prb]]$cycle[i])
      }

      # assign readings prior the first phase and after the last phase
      fzeros <- pyr$date_time < phases[[prb]]$start[1]
      pyr[fzeros, new_phase_col] <- "F0"
      fInfs <- pyr$date_time > phases[[prb]]$stop[nrow(phases[[prb]])]
      pyr[fInfs, new_phase_col] <- "FInf"

      # look for any orphan readings
      NA_check <- is.na(pyr[, new_phase_col])
      if (any(NA_check)) {
        warning(sum(NA_check), " reading(s) in probe ", prb,
                " could not be assigned to a phase!",
                call. = FALSE, immediate. = TRUE)
      }

      # prepare to calculate phase times and assign wait phase
      prb_data <- pyr[, c(1, which(grepl(keyword, colnames(pyr))))]
      # before splitting, save the row order. This is important
      # because F0 can happen across the whole dataset.
      prb_data$row_order <- 1:nrow(prb_data)
      by_phase <- split(prb_data, prb_data[, new_phase_col])

      recipient <- lapply(by_phase, function(x) {
        # calculate time within phase
        x[, new_time_col] <- difftime(x$date_time, x$date_time[1], units = "s")
        if (wait > 0 && grepl("M", x[1, new_phase_col])) {
          # assign wait phase
          these <- x[, new_time_col] < wait
          x[, new_phase_col][these] <- sub("M", "W", x[, new_phase_col][these])
          # adjust measurement phase times.
          aux <- x[!these, ]
          aux$new_times <- difftime(aux$date_time, aux$date_time[1], units = "s")
          x[, new_time_col][!these] <- aux$new_times
        }
        return(x)
      })
      recipient <- do.call(rbind, recipient)
      # remember to restore the original row order!
      recipient <- recipient[order(recipient$row_order), ]

      # transfer updated columns
      pyr[, new_phase_col] <- recipient[, new_phase_col]
      pyr[, new_time_col] <- recipient[, new_time_col]
      # convert difftime to units so everything's in the same format
      pyr[, new_time_col] <- as.numeric(pyr[, new_time_col])
      units(pyr[, new_time_col]) <- "s"

      # export pyrodata object to outside of 
      # the lapply loop to save changes
      pyr <<- pyr

      # export location of columns pertaining to this
      # probe so they can all be grouped in final output
      dvc_prb <- paste0(prb)
      relevant_columns <- grep(dvc_prb, colnames(pyr))
      new_col_order <<- c(new_col_order, relevant_columns)
    } else {
      warning("Phases present for probe ", prb,
              " but no matching pyro data found. ",
              "Disregarding this probe.",
              call. = FALSE, immediate. = TRUE)
    }
  })
  rm(tmp)

  # adjust column positions so probe info is all neatly organized
  pyr <- pyr[, new_col_order]

  input$phased <- pyr

  attributes(input$phased)$wait <- paste(wait, "seconds")

  return(input)
}

#' Replicate phases from one probe to others
#' 
#' Useful when using one flush pump for many chambers/probes
#' 
#' @param input A list containing imported experiment data.
#'   The output from \code{\link{load_experiment}}
#' @param from_device Name of the device from which to copy
#' @param from_channel Name of the channel within "from_device" to duplicate
#' @param to_device Name of the device to which to copy (can be the same)
#' @param to_channel Name of the channel (or channels if vector) 
#'  to which to replicate.
#' 
#' @return The input list with a corrected phases list.
#' 
#' @export
#' 
replicate_phases <- function(input, 
    from_device, from_channel, to_device, to_channel){
  if(is.null(input$phases[[from_device]][[from_channel]])) {
    stop("Could not find specified channel in device ", from_device)
  }

  replacement <- input$phases[[from_device]][[from_channel]]

  for (todev in to_device) {
    for (toCH in to_channel) {
      input$phases[[todev]][[toCH]] <- replacement   
    }
  }

  return(input)
}

#' confirm that the probes listed for each device in the pyro input
#' are present in the phases input
#' 
#' @param pyro_names A vector of column names from the compiled_data object
#' @param phases The phases list.
#' 
#' @return nothing. Used for side effects
#' 
#' @keywords internal
#' 
check_probes_match <- function(pyro_names, phases) {
  pyro_names <- pyro_names[!grepl("date_time", pyro_names)]
  pyro_names <- stringr::str_extract(pyro_names, "(?<=_)[^$]*$")
  probe_names <- unique(pyro_names)
  if (any(!(probe_names %in% names(phases)))) {
    these <- !(probe_names %in% names(phases))
    stop("The could not find all the pyro probe names",
         " in the phases input. Are you sure this is ",
         " the right phases input?",
         " Probes missing phase info: ",
       paste(probe_names[these], collapse = ", "))
  }
}
