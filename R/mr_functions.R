#' Calculate metabolic rates
#'
#' Calculates absolute and mass-specific metabolic rates.
#'
#' @param input a data frame with columns slope_cor (mandatory),
#'  water_vol (mandatory), and animal_mass (optional).
#'
#' @return The input with columns mr_abs and mr_g if animal_mass provided.
#'
#' @export
#'
calc_mr <- function(input){
  input$mr_abs <- -(input$slope_cor * input$water_vol)
  if ("animal_mass" %in% colnames(input)) {
    # see: https://github.com/r-quantities/units/issues/411
    prev_option <- units_options("simplify")
    units_options(simplify = FALSE)
    input$mr_g <- input$mr_abs / input$animal_mass
    units_options(simplify = prev_option)
  }

  return(input)
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
  output <- calc_mr(output)

  # look out for bad slopes
  r2_link <- output$r2 >= r2
  if (any(!r2_link)) {
    message(sum(!r2_link),
            " slopes are below the set r2 threshold and",
            " will be ignored for determination of max mmr")
  }
  # find the max value
  trim_mr <- output[r2_link, ]

  if ("mr_g" %in% colnames(trim_mr)) {
    index <- which.max(trim_mr$mr_g)
  } else {
    index <- which.max(trim_mr$mr_abs)    
  }
  max_mr <- trim_mr[index, , drop = FALSE]

  if (is.null(input$rolling_mr)) {
    input$rolling_mr <- list(values = output, max = max_mr)
  } else {
    input$rolling_mr$values <- rbind(input$rolling_mr$values, output)
    input$rolling_mr$max <- rbind(input$rolling_mr$max, max_mr)
  }

  return(input)
}
