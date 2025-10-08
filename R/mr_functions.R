#' Calculate metabolic rates
#'
#' Calculates the absolute and mass-specific metabolic rates,
#' as well as the relative importance of background (in percentage).
#'
#' @param slope_data  a data frame obtained by using either the
#'  function \code{\link{calc_slopes}} or \code{\link{filter_r2}}
#' @inheritParams conv_w_to_ml
#'
#' @return A data frame with mr_abs, bg, and mr_cor.
#'
#' @export
#'
calc_mr <- function(slope_data, density = 1){
  if ("mass" %in% colnames(slope_data)) {
    vol <- slope_data$volume - conv_w_to_ml(slope_data$mass, density)
    slope_data$mr_abs <- -(slope_data$slope_cor * vol)

    # temporarily drop and reassign units.
    # this is needed to avoid {units} merging oxygen and animal
    # weight when O2 is measured in weight (e.g. mg).
    # see: https://github.com/r-quantities/units/issues/411
    slope_data$mr_g <- drop_units(slope_data$mr_abs) / drop_units(slope_data$mass)
    units(slope_data$mr_g) <- paste0(units(slope_data$mr_abs),
                                     "/", units(slope_data$mass))
  } else {
    slope_data$mr_abs <- -slope_data$slope_cor
  }

  return(slope_data)
}

#' Break down MR from a single cycle
#'
#' Allows determining how oxygen consumption evolved throughout a cycle.
#' Particularly useful for cycles performed after exercise, where metabolic
#' rate might significantly decrease throughout the measurement cycle.
#' Calculates a linear model for each combination of points.
#'
#' @param input The output of \code{\link{process_mr}}
#' @param probe Which probe to select
#' @param cycle which cycle to select
#' @param smoothing How many points should be gathered for
#'  each calculation of the rolling MR
#' @param r2 Minimal coefficient of determination (\eqn{r^{2}})
#'  for valid slopes. Defaults to \eqn{r^{2}} = 0.95.
#' @inheritParams conv_w_to_ml
#'
#' @export
#'
roll_mr <- function(input, probe, cycle, smoothing, density = 1, r2 = 0.95) {
  if (length(probe) > 1) {
    stop("Select only one probe.")
  }
  if (length(cycle) > 1) {
    stop("select only one cycle.")
  }
  if (smoothing < 5) {
    stop("smoothing value is too low.")
  }
  # gather o2 values
  this_probe <- input$trimmed$probe == probe
  if (all(!this_probe)) {
    stop("Could not find specified probe.")
  }

  this_cycle <- input$trimmed$cycle == cycle
  if (all(!this_cycle)) {
    stop("Could not find specified cycle.")
  }
  this_data <- input$trimmed[this_probe & this_cycle, ]
  this_id <- input$probe_info$id[input$probe_info$probe == probe]

  # calculate the lms for each group of seconds
  recipient <- lapply(smoothing:nrow(this_data), function(i) {
    m <- lm(as.numeric(o2_delta) ~ as.numeric(phase_time),
            data = this_data[(i - smoothing + 1):i, ])
    output <- data.frame(id = this_id,
                         probe = probe,
                         cycle = cycle,
                         phase_time = i,
                         smoothing = smoothing,
                         slope = coef(m)[2],
                         r2 = summary(m)$r.squared)
    return(output)
  })
  output <- data.table::rbindlist(recipient)

  # transfer units
  units(output$phase_time) <- units(input$trimmed$phase_time)
  dummy <- this_data$o2_delta[1] / this_data$phase_time[1]
  units(output$slope) <- units(dummy)

  # transfer bg slope for the cycle
  this_probe <- input$slopes$probe == probe
  this_cycle <- input$slopes$cycle == cycle
  output$slope_bg <- input$slopes$slope_bg[this_probe & this_cycle]
  output$slope_cor <- output$slope - output$slope_bg

  # convert to metabolic rate
  the_vol <- this_data$volume[1] - conv_w_to_ml(this_data$mass[1], density)
  output$mr_abs <- -(output$slope_cor * the_vol)
  output$mr_g <- output$mr_abs / this_data$mass[1]

  # look out for bad slopes
  r2_link <- output$r2 >= r2
  if (any(!r2_link)) {
    message(sum(!r2_link),
            " slopes are below the set r2 threshold and",
            " will be ignored for determination of max mmr")
  }
  # find the max value
  trim_mr <- output[r2_link, ]
  index <- which.max(trim_mr$mr_g)
  max_mr <- trim_mr[index, , drop = FALSE]

  if (is.null(input$rolling_mr)) {
    input$rolling_mr <- list(values = output, max = max_mr)
  } else {
    input$rolling_mr$values <- rbind(input$rolling_mr$values, output)
    input$rolling_mr$max <- rbind(input$rolling_mr$max, max_mr)
  }

  return(input)
}
