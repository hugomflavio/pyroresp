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
calc_mr <- function(slope_data, d = 1){
  slope_data$mr_abs <- -(slope_data$slope_cor *
                         (slope_data$volume - conv_w_to_ml(slope_data$mass, d))
                        )

  slope_data$bg <- (slope_data$slope_with_bg - slope_data$slope_cor) / 
                    slope_data$slope_with_bg
  units(slope_data$bg) <- "percent" # automatically multiplies by 100

  slope_data$mr_cor <- slope_data$mr_abs / slope_data$mass

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
#' @inheritParams conv_w_to_ml
#' 
#' @export
#' 
roll_mr <- function(input, probe, cycle, smoothing, d = 1) {
  if (length(probe) > 1) {
    stop("Select only one probe.")
  }
  if (length(cycle) > 1) {
    stop("select only one cycle.")
  }

  this_probe <- input$cleaned$probe == probe
  this_cycle <- input$cleaned$cycle == cycle

  if (all(!this_probe)) {
    stop("Could not find specified probe.")
  }
  if (all(!this_cycle)) {
    stop("Could not find specified cycle.")
  }
  if (all(!this_probe & !this_cycle)) {
    stop("Could not find specified probe and cycle combination.")
  }

  x <- input$cleaned[this_probe & this_cycle, ]
  the_id <- input$probe_info$id[input$probe_info$probe == probe]

  recipient <- lapply(smoothing:nrow(x), function(i) {
    m <- lm(as.numeric(o2_cordelta) ~ as.numeric(phase_time),
            data = x[(i - smoothing + 1):i, ])
    output <- data.frame(id = the_id,
                         probe = probe,
                         cycle = cycle,
                         phase_time = i,
                         smoothing = smoothing,
                         slope_cor = coef(m)[2],
                         r2 = summary(m)$r.squared)
    return(output)
  })
  output <- data.table::rbindlist(recipient)

  units(output$phase_time) <- units(input$cleaned$phase_time)

  dummy <- x$o2_cordelta[1] / x$phase_time[1]
  units(output$slope_cor) <- units(dummy)

  output$mr_abs <- -(output$slope_cor * 
                      (x$volume[1] - conv_w_to_ml(x$mass[1], d))
                    )
  output$mr_cor <- output$mr_abs / x$mass[1]

  # find the max value
  index <- which.max(output$mr_cor)
  max_mr <- output[index, , drop = FALSE]
  
  if (is.null(input$rolling_mr)) {
    input$rolling_mr <- list(values = output, max = max_mr)
  } else {
    input$rolling_mr$values <- rbind(input$rolling_mr$values, output)
    input$rolling_mr$max <- rbind(input$rolling_mr$max, max_mr)
  } 

  return(input)
}

#' Remove part of the mr trace
#' 
#' Useful when there are artifacts along the cycle that could disrupt fractioned
#'  mr calculation. Refinds max mmr after removing the desired segment.
#' 
#' @param input An experiment list containing 
#'  the output of \code{\link{roll_mr}}.
#' @param probe which probe to act upon
#' @param cycle Which phase to act upon
#' @param from second from which to start exclusion
#' @param to second where to end exclusion
#' 
#' @export
#' 
exclude_rolling_mr_segment <- function(input, probe, cycle, from, to) {
  if (length(probe) > 1) {
    stop("Select only one probe.")
  }
  if (length(cycle) > 1) {
    stop("select only one cycle.")
  }

  if (is.null(input$rolling_mr$values)) {
    stop("Could not find rolling mr values in input.")
  }

  this_probe <- input$rolling_mr$values$probe == probe
  this_cycle <- input$rolling_mr$values$cycle == cycle

  if (all(!this_probe)) {
    stop("Could not find specified probe.")
  }
  if (all(!this_cycle)) {
    stop("Could not find specified cycle.")
  }
  if (all(!this_probe & !this_cycle)) {
    stop("Could not find specified probe and cycle combination.")
  }

  phase_time_aux <- as.numeric(input$rolling_mr$values$phase_time)
  rows_out <- this_probe & this_cycle &
              phase_time_aux >= from & 
              phase_time_aux <= to

  cols_to_clean <- c("slope_cor", "r2", "mr_abs", "mr_cor")

  # set to NA instead of removing the rows so ggplot still works nicely.
  input$rolling_mr$values[rows_out, cols_to_clean] <- NA

  # find new max mmr for the probe/cycle
  # remake this_probe and this_cycle with the new row number
  this_probe <- input$rolling_mr$values$probe == probe
  this_cycle <- input$rolling_mr$values$cycle == cycle
  index <- which.max(input$rolling_mr$values$mr_cor[this_probe & this_cycle])

  # deposit it on the max side
  # the index is only right for this probe/cycle combo
  aux_table <- input$rolling_mr$values[this_probe & this_cycle, ] 
  this_row <- aux_table[index, ]

  this_probe <- input$rolling_mr$max$probe == probe
  this_cycle <- input$rolling_mr$max$cycle == cycle
  input$rolling_mr$max[this_probe & this_cycle, ] <- this_row

  return(input)
}
