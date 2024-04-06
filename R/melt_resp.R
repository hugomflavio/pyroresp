#' Convert a human-friendly respirometry table to a computer-friendly format
#' 
#' Also appends additional probe info if provided.
#' 
#' @param input A dataframe with oxygen, temperature, pressure, and phase data.
#'  The output of \code{\link{merge_pyro_phases}}
#' @inheritParams load_experiment
#' 
#' @export
#' 
melt_resp <- function(input, probe_info) {

  ox_aux <- reshape2::melt(input, 
    id.vars = "date_time",
    measure.vars = colnames(input)[grepl("ox_", colnames(input))],
    value.name = "o2",
  )
  ox_aux$probe <- stringr::str_extract(ox_aux$variable, "(?<=ox_).*")
  ox_aux$variable <- NULL
  colnames(ox_aux)[1] <- "date_time"
  ox_aux <- ox_aux[, c("date_time", "probe", "o2")]

  pressure_aux <- reshape2::melt(input, 
    id.vars = "date_time",
    measure.vars = colnames(input)[grepl("pressure_", colnames(input))],
    value.name = "pressure"
  )
  pressure_aux$variable <- NULL
  pressure_aux$date_time <- NULL

  temp_aux <- reshape2::melt(input, 
    id.vars = "date_time",
    measure.vars = colnames(input)[grepl("temp_", colnames(input))],
    value.name = "temp"
  )
  temp_aux$variable <- NULL
  temp_aux$date_time <- NULL

  sal_aux <- reshape2::melt(input, 
    id.vars = "date_time",
    measure.vars = colnames(input)[grepl("sal_", colnames(input))],
    value.name = "sal"
  )
  sal_aux$variable <- NULL
  sal_aux$date_time <- NULL

  phase_aux <- reshape2::melt(input, 
    id.vars = "date_time",
    measure.vars = colnames(input)[grepl("phase_", colnames(input))],
    value.name = "phase"
  )
  phase_aux$variable <- NULL
  phase_aux$date_time <- NULL

  pre_output <- cbind(ox_aux, pressure_aux, temp_aux, sal_aux, phase_aux)

  if (any(grepl("ph_", colnames(input)))) {
    pre_output$ph <- as.vector(t(input[, grepl("ph_", colnames(input))]))
  }

  pre_output$cycle <- as.numeric(sub("F|M", "", as.character(pre_output$phase)))

  # if fish information is provided
  if (!is.null(probe_info)) {
    # discard any data that does not match the probes we want.
    pre_output <- pre_output[pre_output$probe %in% probe_info$probe, ]

    if (nrow(pre_output) == 0) {
      stop("The probes listed in probe_info don't match any probe listed on ",
           "the data. Aborting.")
    }

    # Then include the extra information
    link <- match(pre_output$probe, probe_info$probe)
    pre_output <- cbind(pre_output, probe_info[link, c("id", "mass", "volume")])
  
    # and trim away the cycles that happen before the first_cycle
    by_probe <- split(pre_output, pre_output$probe)
    
    recipient <- lapply(names(by_probe), function(the_probe, probe_info) {

      trimmed_db <- by_probe[[the_probe]]

      first_cycle <- probe_info$first_cycle[probe_info$probe == the_probe]
      if (!is.na(first_cycle)) {
        trimmed_db <- trimmed_db[trimmed_db$cycle >= first_cycle, ]
      } else {
        trimmed_db <- NULL
      }
      
      return(trimmed_db) 

    }, probe_info = probe_info)

    output <- as.data.frame(data.table::rbindlist(recipient))
  } else {
    output <- as.data.frame(pre_output)
  }

  output$phase <- as.character(output$phase)

  return(output)
}



