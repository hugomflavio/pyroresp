#' Calculation of Metabolic Rate
#'
#' The function is used to calculate and plot background respiration, absolute and mass-specific metabolic rates.
#'
#' @param slope.data  a data frame obtained by using the function \code{\link{extract.slope}}
#' @param density  numeric: the density of an animal body (\eqn{kg/m^{3}})
#'
#' @return The function returns a data frame with calculated background respiration, absolute and mass-specific metabolic rates.
#'
#' @export

calculate.MR  <- function(slope.data, density = 1000){

  BW = slope.data$Mass/1000 # convert Mass from 'g' to 'kg'
  V = slope.data$Volume/1000 - (BW/(density/1000)) #convert Volume from 'mL' to 'L' and Density from 'kg/M^3' to 'kg/L'

  slope.data$MR.abs.with.BR = -(slope.data$Slope.with.BR*V*3600) #3600 sec = 1 hour

  slope.data$BR = (slope.data$Slope.with.BR - slope.data$Slope)/slope.data$Slope.with.BR*100 # in %
  slope.data$MR.abs = -(slope.data$Slope*V*3600)
  slope.data$MR.mass = -(slope.data$Slope*V*3600/BW)

  # a <- slope.data$DO.unit[1]
  # MR.data <- slope.data[,-12]
  # MR.data$DO.unit <- a[1]
  MR.data <- slope.data

  attributes(MR.data)$selection.method <- attributes(slope.data)$selection.method
  return(MR.data)
}




#' Break down MR during a single cycle
#' 
#' Allows determining how oxygen consumption evolved throughout a cycle. Particularly useful for
#' cycles performed after exercise, where metabolic rate might significantly decrease throughout
#' the measurement cycle. Calculates a linear model for each combination of points.
#' 
#' @param input The output of process_pyro_mr
#' @param probe Which probe to select
#' @param phase which phase to select
#' @param smoothing How many points should be gathered for each calculation of the rolling MR
#' 
#' @export
#' 
fraction_mr <- function(input, probe, phase, smoothing) {
  x <- input$corrected[input$corrected$Phase == phase & input$corrected$Probe == probe, ]
  head(x)

  recipient <- lapply(smoothing:nrow(x), function(i) {
    aux_m <- lm(O2.delta.corrected ~ Phase.Time,
          data = x[(i - smoothing + 1):i, ])
    output <- data.frame(
      Probe = probe,
      Phase = phase,
      Sec = i,
      Smoothing = smoothing,
      Slope = coef(aux_m)[2],
      R2 = summary(aux_m)$r.squared)
    return(output)
  })

  output <- data.table::rbindlist(recipient)

  BW <- x$Mass[1]/1000
  V <- x$Volume[1]/1000 - BW

  output$MR.abs <- -(output$Slope * V * 3600)
  output$MR.mass <- -(output$Slope * V * 3600 / BW)

  return(output)
}

#' Wrapper of fraction_mr for MMR cycles
#'
#' @inheritParams fraction_mr
#' 
#' @export
#' 
compile_rolling_mmr <- function(input, smoothing) {
  detailed_mmr_list <- lapply(1:nrow(input$mmr), function(i) {
    if (!is.na(input$mmr$Phase[i])) {
      output <- fraction_mr(input, input$mmr$Probe[i], input$mmr$Phase[i], smoothing = smoothing)
      output$MR.mass.umol.g <- output$MR.mass/1000
      output$ID <- input$chamber.info$ID[i]
    } else {
      output <- NULL
    }
    return(output)
  })

  detailed_mmr <- as.data.frame(do.call(rbind, detailed_mmr_list))

  max_mmr_list <- lapply(1:nrow(input$mmr), function(i) {
    index <- which.max(detailed_mmr_list[[i]]$MR.mass.umol.g)
    detailed_mmr_list[[i]][index, , drop = FALSE]
  })

  max_mmr <- as.data.frame(do.call(rbind, max_mmr_list))

  max_mmr <- merge(input$mmr[, c("Probe", "ID")], max_mmr[, !(colnames(max_mmr) %in% "ID")], by = "Probe", all = TRUE)

  return(list(detailed_mmr = detailed_mmr, max_mmr = max_mmr))
}

#' Remove part of the mr trace
#' 
#' Useful when there are artifacts along the cycle that could disrupt fractioned mr calculation.
#' 
#' @param rolling_mmr the output of compile_rolling_mmr
#' @param probe which probe to act upon
#' @param start second from which to start exclusion
#' @param end second where to end exclusion
#' 
#' @export
#' 
exclude_rolling_mmr_segment <- function(rolling_mmr, probe, start, end) {
  to_exclude <- rolling_mmr$detailed_mmr$Probe == probe & 
                rolling_mmr$detailed_mmr$Sec >= start & 
                rolling_mmr$detailed_mmr$Sec <= end

  rolling_mmr$detailed_mmr[to_exclude, 5:9] <- NA

  to_check <- rolling_mmr$detailed_mmr[rolling_mmr$detailed_mmr$Probe == probe, ]
  index <- which.max(to_check$MR.mass.umol.g)
  
  rolling_mmr$max_mmr[rolling_mmr$max_mmr$Probe == probe, ] <- to_check[index, colnames(rolling_mmr$max_mmr)]

  return(rolling_mmr)
}


#' dummy doc
#' 
#' @export
#' 
calculate.bg <- function(input, O2_col, method = c('mean', 'first', 'last'), force.linear = TRUE, smoothing = 30){

  method <- match.arg(method)

  n.phases <- unique(input$Phase)

  if (method == 'first')
    input <- input[input$Phase == n.phases[1], ]
    
  if (method == 'last')
    input <- input[input$Phase == n.phases[length(n.phases)], ]

  chamber.lists <- split(input, input$Probe)

  bg.lists <- lapply(names(chamber.lists), function(chamber) {
    # cat(chamber, '\n')
    sub.data <- input[input$Probe == chamber, ]

    if (force.linear) {

      bg.lm <- eval(parse(text = paste("lm(", O2_col, "~ Phase.Time, data = sub.data)")))
      bg.lm$coefficients[1] <- 0 

      output <- data.frame(Phase.Time = 1:max(sub.data$Phase.Time))
      output$O2.background <- as.vector(stats::predict(bg.lm, output, type = "response", se.fit = FALSE))
    }
    else {
      output <- stats::aggregate(sub.data[, O2_col], by = list(sub.data$Phase.Time), mean, na.rm = TRUE)
      colnames(output) <- c('Phase.Time', 'O2.background')
      
      if (smoothing > 1) {
        x <- stats::filter(output$O2.background, rep(1/smoothing, smoothing), sides = 2)

        if (smoothing%%2 == 0) {
          for (i in 1:(smoothing/2-1)) {
            r <- round(i/2)
            x[i] <- mean(output$O2.background[(i-r):(i+r)])
          }

          last.value <- length(x) - round(smoothing/2)
          smooth.values <- (last.value - round(smoothing/2)):last.value
          smooth.slope <- mean(x[smooth.values-1]-x[smooth.values])

          for (i in (length(x)-(smoothing/2)):length(x)) {
            x[i] <- x[i-1] - smooth.slope
          }
        }
        output$O2.background <- x
      }
    }  
    return(output)
  })
  names(bg.lists) <- names(chamber.lists)

  output <- as.data.frame(data.table::rbindlist(bg.lists, idcol = 'Probe'))

  return(output)
}


#' Override background correction factor for any given
#' chambers using the values from another chamber.
#' 
#' @export
#' 
override_bg <- function(bg, replace, with) {
  if (!is.data.frame(bg))
    stop("bg must be a data.frame")

  if (all(!grepl("Probe", colnames(bg))))
    stop("bg must contain a 'Probe' column.")
  
  if (all(!grepl("O2.background", colnames(bg))))
    stop("bg must contain a 'O2.background' column.")
  
  if (any(!(replace %in% bg$Probe)))
    stop("Could not find some of the specified chambers to replace in bg")

  if (length(with) != 1)
    stop("Please chose only one chamber to use as replacement in 'with'.")
  
  if (!(with %in% bg$Probe))
    stop("Could not find the replacement chamber in bg")

  for (i in replace) {
    if (sum(bg$Probe == i) > sum(bg$Probe == with))
      stop("the cycle for the replacement chamber is shorted than the cycle for the chamber to be replaced.")

    bg$O2.background[bg$Probe == i] <- bg$O2.background[bg$Probe == with][1:sum(bg$Probe == i)]
    # the additional subset at the end ensures 
  }

  return(bg)
}



#' Dummy documentation
#' 
#' @export
#' 
calc_delta <- function(input, O2_col) {

  input[, paste0('O2.delta', sub('O2', '', O2_col))] <- NULL
  
  by.chamber <- split(input, input$Probe)
  
  recipient <- lapply(names(by.chamber), function(the.chamber, info.data) {

    trimmed.db <- by.chamber[[the.chamber]]

    by.phase <- split(trimmed.db, trimmed.db$Phase)
    
    recipient <- lapply(by.phase, function(the.phase) {
      the.phase$O2.delta <- the.phase[, O2_col] - the.phase[1, O2_col]
      return(the.phase)
    })
    
    output <- data.table::rbindlist(recipient)

  })

  output <- as.data.frame(data.table::rbindlist(recipient))
  colnames(output)[ncol(output)] <- paste0('O2.delta', sub('O2', '', O2_col))

  return(output)
}


#' Prepare imported measurements for further analyses
#' 
#' @param input the data frame containing imported oxygen measurements
#' @param wait integer: the number of first rows for each measurement phase (M) which should be reassigned to the wait phase (W). The parameter should be used when the wait phase (W) is absent (e.g. in 'Q-box Aqua' logger software) or not long enough to eliminate non-linear change in DO concentration over time from the measurement phase (M) after shutting off water supply from the ambient water source.
#' 
#' @return A data.frame containing measurements valid for further analyses
#' 
#' @export
#' 
clean.meas <- function(input, wait = 0, auto.cut.last = FALSE){

  input$Date <- as.Date(input$Date.Time)
  input$Real.Time <- chron::times(strftime(input$Date.Time, "%H:%M:%S"))

  # Removing Non-Measurement Data  
  input <- input[grepl("^M", input$Phase), ]
  
  # this is used just to ensure that the phases maintain their order,
  # even it they don't start at one or are not sorted at the start.
  phase.order <- as.numeric(gsub("[M]", "", unique(input$Phase)))
  phase.order <- order(phase.order)
  
  input$Phase <- factor(input$Phase, levels = unique(input$Phase)[phase.order])


  #the rest has to be done on a chamber by chamber basis.

  by.chamber <- split(input, input$Probe)

  recipient <- lapply(names(by.chamber), function(the.chamber) {

    trimmed.db <- by.chamber[[the.chamber]]

    # Removing the final measurement Phase if necessary or forced (tail error)
    rows.per.phase <- table(trimmed.db$Phase)

    if (length(rows.per.phase) > 1)
      mean.rows.per.phase <- mean(rows.per.phase[-length(rows.per.phase)])
    else
      mean.rows.per.phase <- rows.per.phase

    if (tail(rows.per.phase, 1) < wait | auto.cut.last) {
      trimmed.db <- trimmed.db[trimmed.db$Phase != tail(levels(trimmed.db$Phase), 1), ]
      trimmed.db$Phase <- droplevels(trimmed.db$Phase)
    }


    # cut off first n rows from 'M' phase
    if(wait != 0){
      # the code below grabs the 1:nrow vector, breaks it out by phase, and then uses tail() with
      # a negative n to grab all numbers but the first 30 that show up for each phase.
      index <- unlist(tapply(1:nrow(trimmed.db), trimmed.db$Phase, tail, -(wait)), use.names = FALSE)
      trimmed.db <- trimmed.db[index, ]
    }

    # reset rownames
    row.names(trimmed.db) <- 1:nrow(trimmed.db)

    # and now, by phase, calculate the passing time and store the start and end points of that measurement
    aux <- split(trimmed.db, trimmed.db$Phase)

    aux <- aux[sapply(aux, nrow) > 0]

    aux <- lapply(aux, function(x) {
      x$Start.Meas <- x$Real.Time[1]
      x$End.Meas <- x$Real.Time[nrow(x)]
      x$Phase.Time <- as.numeric(difftime(x$Date.Time, x$Date.Time[1], units = 's'))

      return(x)
    })
    trimmed.db <- as.data.frame(data.table::rbindlist(aux))
 
    return(trimmed.db)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  return(output)
}




#' Correction of Metabolic Rate Measurements
#'
#' The function is used to correct metabolic rate measurements for background respiration. To this end, oxygen consumption is estimated as the slope of the linear regression of measured \eqn{O_{2}} concentration over time, and is extracted for background respiration test and for each measurement phase. The correction is based on subtraction of oxygen consumption obtained during background respiration test from oxygen consumption obtained during metabolic rate measurements.
#'
#' @param pre.bg  a data frame obtained by using the function \code{\link{calculate.bg}} for a blank test before actual metabolic rate measurements
#' @param post.bg  a data frame obtained by using the function \code{\link{calculate.bg}} for a blank test after actual metabolic rate measurements
#' @param meas.data  a data frame obtained by using the function \code{\link{process_pyro_files}} for actual metabolic rate measurements
#' @param method  string: the name of the method used for background respiration correction:
#' @param O2_col the name of the column containing the raw O2 values
#' @param O2_delta_col The name of the column containing the O2 delta values
#'
#' #' \itemize{
#' \item  "pre.test" - subtracts oxygen consumption of pre.data from oxygen consumptions of meas.data
#' \item  "post.test" - subtracts oxygen consumption of post.data from oxygen consumptions of meas.data
#' \item  "average" - subtracts an averaged oxygen consumption of pre.data and\cr post.data from oxygen consumptions of meas.data
#' \item  "linear" - subtracts a vector of progressively changing microbial consumptions from oxygen consumptions of meas.data. The values of oxygen consumption are linearly predicted from two reference points: oxygen consumption of pre.data and oxygen consumption of post.data.
#' \item  "exponential" - subtracts a vector of progressively changing microbial consumptions from oxygen consumptions of meas.data. The values of oxygen consumption are exponentially predicted from two reference points: oxygen consumption of pre.data and oxygen consumption of post.data.
#' \item  "parallel" - subtracts oxygen consumption in an empty chamber from oxygen consumptions of meas.data for each chamber
#' \item  "none" - does not perform oxygen consumption subtraction. Not recommended for anything other than checking test data.
#' }
#' @param empty.chamber  string: the name of an empty chamber used only for the method 'parallel'
#'
#' @return  The function returns a data frame containing data of metabolic rate measurements corrected for background respiration.
#'
#'
#' AMR.clean <- correct.meas(post.data = post,
#'                           meas.data = AMR.raw,
#'                           method = "post.test")
#'
#' @references {Svendsen, M. B. S., Bushnell, P. G., & Steffensen, J. F. (2016). Design and setup of intermittent-flow respirometry system for aquatic organisms. Journal of Fish Biology, 88(1), 26-50.}
#'
#' @export
#'
correct.meas <- function (pre.bg, post.bg, meas.data, O2_col = 'O2.raw', O2_delta_col = "O2.delta.raw",
                        method = c("pre.test", "post.test", "average",
                                   "linear", "exponential", "parallel", "none"),
                        empty.chamber){

  method <- match.arg(method)

  M.total <- max(meas.data$Cycle)

  if (method == "pre.test") {
    my_bg <- pre.bg

    meas.data$temporary_index <- paste(meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Probe, my_bg$Phase.Time)
  }

  if (method == "post.test") {
    my_bg <- post.bg

    meas.data$temporary_index <- paste(meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Probe, my_bg$Phase.Time)
  }

  if (method == "average") {
    stop("average has not been updated yet")  
}

  if (method == "linear") {
    my_bg <- calculate_linear_bg_progression(pre.bg = pre.bg, post.bg = post.bg, meas.data = meas.data)

    meas.data$temporary_index <- paste(meas.data$Cycle, meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Cycle, my_bg$Probe, my_bg$Phase.Time)
  }

  if (method == "exponential") {
    stop("exponential has not been updated yet")
  }

  if (method == "parallel") {
    if (!any(empty.chamber %in% meas.data$Probe)) {
      stop("Could not find bg probe in data.")
    }
    my_bg <- meas.data[meas.data$Probe == empty.chamber, ]
    my_bg$O2.background <- my_bg[, O2_delta_col]
    my_bg$temporary_index <- paste(my_bg$Cycle, my_bg$Phase.Time)

    meas.data$temporary_index <- paste(meas.data$Cycle, meas.data$Phase.Time)
  }

  if(method == "none") {
    meas.data$O2.background <- 0
  } else {
    link <- match(meas.data$temporary_index, my_bg$temporary_index)

    meas.data$O2.background <- my_bg$O2.background[link]

    meas.data$temporary_index <- NULL
    
    if (any(is.na(meas.data$O2.background)))
      warning('Some measurement phases are longer than the background. Discarding overextended points', immediate. = TRUE, call. = FALSE)

    meas.data <- meas.data[!is.na(meas.data$O2.background), ]
  }

  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  meas.data$O2.corrected <- meas.data[, O2_col] - meas.data$O2.background

  aux <- split(meas.data, paste0(meas.data$Probe, meas.data$Phase))
  aux <- lapply(aux, function(x) {
    x$O2.delta.corrected <- x$O2.corrected - x$O2.corrected[1]
    return(x)
  })

  output <- as.data.frame(data.table::rbindlist(aux))
  
  attributes(output)$correction_method <- method
  
  return(output)
}






#' Internal function
#' 
#' @keywords internal
#' 
calculate_linear_bg_progression <- function(pre.bg, post.bg, meas.data) {
  cycles <- max(meas.data$Cycle)

  my_bg <- lapply(unique(meas.data$Probe), function(chamber) {
    # cat(chamber, "\n")
    pre_line <- pre.bg$O2.background[pre.bg$Probe == chamber]

    post_line <- post.bg$O2.background[post.bg$Probe == chamber]

    if (length(pre_line) != length(post_line)) {
      warning("The pre-bg and post-bg in chamber ", chamber, " have different lengths! Truncating longer vector.", immediate. = TRUE, call. = FALSE)
      
      if (length(pre_line) < length(post_line))
        post_line <- post_line[1:length(pre_line)]
      else
        pre_line <- pre_line[1:length(post_line)]
    }

    difference <- post_line - pre_line

    if (all(difference == 0)) {
      stop("The pre-bg and post-bg are exactly the same! Cannot calculate linear progression.", call. = FALSE)
    }

    # linear increments.
    increments <- difference / (cycles - 1)

    my_list <- lapply(1:length(pre_line), function(i) {
      seq(from = pre_line[i], to = post_line[i], by = increments[i])
    })

    my_matrix <- as.data.frame(do.call(rbind, my_list))

    return(my_matrix)
  })
  names(my_bg) <- unique(meas.data$Probe)



  my_bg_simplified <- lapply(names(my_bg), function(chamber) {
    x <- suppressMessages(reshape2::melt(my_bg[[chamber]]))
    colnames(x) <- c("Cycle", "O2.background")
    x$Phase.Time <- 1:nrow(my_bg[[chamber]])
    x$Cycle <- as.numeric(sub("V", "", x$Cycle))
    x$Probe <- chamber
    return(x)
  })

  my_bg_df <- do.call(rbind, my_bg_simplified)
  return(my_bg_df)
}



# calculate_plateau_bg_progression <- function(pre.bg, post.bg, meas.data) {
#   cycles <- max(meas.data$Cycle)

#   my_bg <- lapply(unique(meas.data$Probe), function(chamber) {
#     # cat(chamber, "\n")
#     pre_line <- pre.bg$O2.background[pre.bg$Probe == chamber]

#     post_line <- post.bg$O2.background[post.bg$Probe == chamber]

#     if (length(pre_line) != length(post_line)) {
#       warning("The pre-bg and post-bg in chamber ", chamber, " have different lengths! Truncating longer vector.", immediate. = TRUE, call. = FALSE)
      
#       if (length(pre_line) < length(post_line))
#         post_line <- post_line[1:length(pre_line)]
#       else
#         pre_line <- pre_line[1:length(post_line)]
#     }

#     difference <- post_line - pre_line

#     if (all(difference == 0)) {
#       stop("The pre-bg and post-bg are exactly the same! Cannot calculate linear progression.", call. = FALSE)
#     }

#     increments <- difference / (cycles - 1)

#     my_list <- lapply(1:length(pre_line), function(i) {
#       seq(from = pre_line[i], to = post_line[i], by = increments[i])
#     })

#     my_matrix <- as.data.frame(do.call(rbind, my_list))

#     return(my_matrix)
#   })
#   names(my_bg) <- unique(meas.data$Probe)



#   my_bg_simplified <- lapply(names(my_bg), function(chamber) {
#     x <- suppressMessages(reshape2::melt(my_bg[[chamber]]))
#     colnames(x) <- c("Cycle", "O2.background")
#     x$Phase.Time <- 1:nrow(my_bg[[chamber]])
#     x$Cycle <- as.numeric(sub("V", "", x$Cycle))
#     x$Probe <- chamber
#     return(x)
#   })

#   my_bg_df <- do.call(rbind, my_bg_simplified)
#   return(my_bg_df)
# }


# i = 100

# L = post_line[i] # maximum value
# k = 0.1 # growth rate/steepness
# x = 0:100 # the x values
# x0 = 0 # sigmoid midpoint

# y = L/(1+exp(-k*(x-x0)))
# plot(y)



# loadf("convertLink")