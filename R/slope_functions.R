#' Calculate slopes from clean measurements
#' 
#' @param input A dataframe of clean measurements.
#'  The output of \code{\link{clean_meas}}.
#' 
#' @export
#' 
calc_slopes <- function(input) {
  # The operation is done by phase and by_probe,
  # so the dataset is broken twice below
  by_probe <- split(input, input$probe)

  recipient <- lapply(by_probe, function(the_probe) {
    by_phase <- split(the_probe, the_probe$phase)
    by_phase <- by_phase[sapply(by_phase, nrow) > 0]

    recipient <- lapply(by_phase, function(phase) {
        output <- data.frame(probe = phase$probe[1],
                             id = phase$id[1],
                             mass = phase$mass[1],
                             volume = phase$volume[1],
                             date_time = phase$date_time[nrow(phase)],
                             phase = phase$phase[1],
                             cycle = phase$cycle[1],
                             temp = mean(phase$temp))

      if (all(is.na(phase$o2_cordelta))) {
        output$slope_with_bg <- NA
        output$slope_cor <- NA
        output$se <- NA
        output$r2 <- NA
      } else {

        # lm and units don't play along well if you intend
        # to ask for a summary. Must drop units and reattach
        # them later until this is fixed.
        m_with_bg <- lm(as.numeric(o2_delta) ~ as.numeric(phase_time), 
                        data = phase)
        m_cor <- lm(as.numeric(o2_cordelta) ~ as.numeric(phase_time),
                        data = phase)

        output$slope_with_bg <- coef(m_with_bg)[2]
        output$slope_cor <- coef(m_cor)[2]
        output$se <- summary(m_cor)$coef[4]
        output$r2 <- summary(m_cor)$r.squared
      }
      #make a dummy variable to let {units} figure the units out by itself
      dummy <- phase$o2_delta[1] / phase$phase_time[1]
      units(output$slope_with_bg) <- units(dummy)
      
      dummy <- phase$o2_cordelta[1] / phase$phase_time[1]
      units(output$slope_cor) <- units(dummy)
     
      return(output)
    })
    # start rebinding back to a data table
    output <- data.table::rbindlist(recipient)
    return(output)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  if (any(is.na(output$slope_cor))) {
    warning("Could not calculate slopes for ", sum(is.na(output$slope_cor)), 
            " probe-cycle combinations.",
            immediate. = TRUE, call. = FALSE)
  }

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

