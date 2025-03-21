#' Calculate slopes from clean measurements
#'
#' @param input A dataframe of clean measurements.
#'  The output of \code{\link{trim_resp}}.
#'
#' @export
#'
calc_slopes <- function(input) {
  # The operation is done by cycle and by probe,
  # so the dataset is broken twice below
  by_probe <- split(input$trimmed, input$trimmed$probe)

  recipient <- lapply(by_probe, function(the_probe) {
    by_cycle <- split(the_probe, the_probe$cycle)
    by_cycle <- by_cycle[sapply(by_cycle, nrow) > 0]

    recipient <- lapply(by_cycle, function(the_cycle) {
        output <- data.frame(probe = the_cycle$probe[1],
                             id = the_cycle$id[1],
                             mass = the_cycle$mass[1],
                             volume = the_cycle$volume[1],
                             date_time = the_cycle$date_time[nrow(the_cycle)],
                             phase = the_cycle$phase[1],
                             cycle = the_cycle$cycle[1],
                             temp = mean(the_cycle$temp, na.rm = TRUE))

      if (all(is.na(the_cycle$o2_delta))) {
        output$slope <- NA
        output$se <- NA
        output$r2 <- NA
      } else {
        # lm and units don't play along well if you intend
        # to ask for a summary. Must drop units and reattach
        # them later until this is fixed.
        m <- lm(as.numeric(o2_delta) ~ as.numeric(phase_time),
                data = the_cycle)

        output$slope <- coef(m)[2]
        output$se <- summary(m)$coef[4] # why do we care about this? IDK...
        output$r2 <- summary(m)$adj.r.squared
      }
      # make a dummy variable to let {units} figure itself out
      dummy <- the_cycle$o2_delta[1] / the_cycle$phase_time[1]
      units(output$slope) <- units(dummy)

      return(output)
    })
    # start rebinding back to a data table
    output <- data.table::rbindlist(recipient)
    return(output)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  if (any(is.na(output$slope))) {
    warning("Could not calculate slopes for ", sum(is.na(output$slope)),
            " probe-cycle combinations.",
            immediate. = TRUE, call. = FALSE)
  }

  input$slopes <- output
  return(input)
}

#' Calculate slopes from clean measurements
#'
#' @param input A resp object with a trimmed table.
#' @param probe The probe for which to calculate a slope
#' @param cycle The cycle for which to calculate a slope
#' @param max_duration Number of seconds to skip at the start
#'   of the phase.
#' @param max_duration Only use points up to this time
#'   into the phase (in seconds). Defaults to Inf.
#'
#' @export
#'
calc_single_slope <- function(input, probe, cycle, skip = 0, max_duration = Inf) {

  if (length(probe) > 1) {
    stop("Choose only one probe", call. = FALSE)
  }
  if (!(probe %in% unique(input$trimmed$probe))) {
    stop("Chosen probe does not exist. Available probes: ",
         paste(names(by_probe), collapse = ", "),
         call. = FALSE)
  }
  my_probe <- input$trimmed[input$trimmed$probe == probe, ]

  if (length(cycle) > 1) {
    stop("Choose only one cycle", call. = FALSE)
  }
  if (!(cycle %in% unique(my_probe$cycle))) {
    stop("Chosen probe does not exist. Available probes: ",
         paste(names(by_probe), collapse = ", "),
         call. = FALSE)
  }

  the_cycle <- my_probe[my_probe$cycle == cycle, ]

  if ((max_duration - skip) < 10) {
    warning("Calculating very short slopes is unadvisable.",
            immediate. = TRUE, call. = FALSE)
  }

  start_here <- as.numeric(the_cycle$phase_time) >= skip
  end_here <- as.numeric(the_cycle$phase_time) <= max_duration
  the_cycle <- the_cycle[start_here & end_here, ]

  output <- data.frame(probe = the_cycle$probe[1],
                       id = the_cycle$id[1],
                       mass = the_cycle$mass[1],
                       volume = the_cycle$volume[1],
                       date_time = the_cycle$date_time[nrow(the_cycle)],
                       phase = the_cycle$phase[1],
                       cycle = the_cycle$cycle[1],
                       temp = mean(the_cycle$temp, na.rm = TRUE))

  if (all(is.na(the_cycle$o2_delta))) {
    output$slope <- NA
    output$se <- NA
    output$r2 <- NA
  } else {
    # lm and units don't play along well if you intend
    # to ask for a summary. Must drop units and reattach
    # them later until this is fixed.
    m <- lm(as.numeric(o2_delta) ~ as.numeric(phase_time),
            data = the_cycle)

    output$slope <- coef(m)[2]
    output$se <- summary(m)$coef[4] # why do we care about this? IDK...
    output$r2 <- summary(m)$adj.r.squared
  }
  # make a dummy variable to let {units} figure itself out
  dummy <- the_cycle$o2_delta[1] / the_cycle$phase_time[1]
  units(output$slope) <- units(dummy)

  return(output)
}

#' Filter  slopes by threshold R2
#'
#' @param slopes  a data frame obtained by using
#'  the function \code{\link{calc_slopes}}
#' @param r2  numeric: minimal coefficient of determination (\eqn{r^{2}})
#'  for valid slopes. Default \eqn{r^{2}} = 0.95.
#'
#' @return A data frame with the information about extracted slopes.
#'
#' @export
#'
filter_r2 <- function(slopes, r2 = 0.95){
  to_exclude <- sum(slopes$r2 < r2, na.rm = TRUE)
  total <- sum(!is.na(slopes$r2))

  to_keep <- slopes$r2 >= r2
  to_keep[is.na(to_keep)] <- FALSE

  slopes <- slopes[to_keep, ]

  message("Excluded ", to_exclude, " slope(s) (",
          round(to_exclude / total * 100, 1),
          "%).")

  attributes(slopes)$r2_threshold <- r2
  attributes(slopes)$n_discarded <- to_exclude

  return(slopes)
}

