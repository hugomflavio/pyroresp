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
  load(target, envir = e)
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
    x$stop <- c(x$timestamp[-1] - 1, Inf)
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
#' 
#' @return An updated experiment list containing a phased data frame
#' 
#' @export
#' 
assign_phases <- function(input) {
  pyr <- input$pyro$compiled_data
  phases <- input$phases
  check_probes_match(colnames(pyr), phases)

  new_col_order <- 1
  tmp <- lapply(names(phases), function(prb) {
    # create column with placeholders if probe is in the data
    keyword <- paste0("_", prb, "$")
    if (any(grepl(keyword, colnames(pyr)))) {

      new_col <- paste0("phase_", prb)
      pyr[, new_col] <- "F0"

      # assign phases
      for (i in 1:nrow(phases[[prb]])) {
        check1 <- pyr$date_time >= phases[[prb]]$start[i]
        check2 <- pyr$date_time <= phases[[prb]]$stop[i]
        this_phase <- check1 & check2

        pyr[this_phase, new_col] <- paste0(phases[[prb]]$phase[i],
                                           phases[[prb]]$cycle[i])
      }

      NA_check <- is.na(pyr[, new_col])
      if (any(NA_check)) {
        warning(sum(NA_check), " measurement(s) in device ",dvc,
                ", probe ", prb,
                " could not be assigned to a phase!",
                call. = FALSE, immediate. = TRUE)
      }

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
  
  pyr <- pyr[, new_col_order]

  input$phased <- pyr

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
