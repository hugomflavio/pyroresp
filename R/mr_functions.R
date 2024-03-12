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
#' @param phase which phase to select
#' @param smoothing How many points should be gathered for 
#'  each calculation of the rolling MR
#' @inheritParams conv_w_to_ml
#' 
#' @export
#' 
roll_mr <- function(input, probe, phase, smoothing, d = 1) {
  if (length(probe) > 1) {
    stop("Select only one probe.")
  }
  if (length(phase) > 1) {
    stop("select only one phase.")
  }

  this_phase <- input$cleaned$phase == phase
  this_probe <- input$cleaned$probe == probe

  if (all(!this_phase)) {
    stop("Could not find specified phase.")
  }
  if (all(!this_probe)) {
    stop("Could not find specified probe.")
  }
  if (all(!this_phase & !this_probe)) {
    stop("Could not find specified phase and probe combination.")
  }

  x <- input$cleaned[this_phase & this_probe, ]
  the_id <- input$probe_info$id[input$probe_info$probe == probe]

  recipient <- lapply(smoothing:nrow(x), function(i) {
    m <- lm(as.numeric(o2_cordelta) ~ as.numeric(phase_time),
            data = x[(i - smoothing + 1):i, ])
    output <- data.frame(id = the_id,
                         probe = probe,
                         phase = phase,
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
