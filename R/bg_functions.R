#' Calculate background values
#' 
#' Can calculate linear background with forced 0 intercept, or rolling average
#' background if smoothing is specified.
#' 
#' @param input The output of \code{\link{process_experiment}}.
#' @param method How should multiple background cycles for a single probe
#'  be handled? One of:
#'    'mean'  : All the measured cycles are used for the calculations
#'    'first' : Only the first cycle is used for the calculations
#'    'last'  : Only the last cycle is used for the calculations
#' @param smoothing Optional. If included, the function calculates a rolling
#'  average background, instead of the classical linear model with 0-intercept.
#'  Value in seconds.
#' 
#' @return The input list, with a new "bg" object.
#' 
#' @export
#' 
calc_bg <- function(input, method = c('mean', 'first', 'last'),
    smoothing){

  method <- match.arg(method)

  # extract some variables we'll be using for convenience
  cleaned <- input$cleaned
  phases <- unique(cleaned$phase)

  # discard unnecessary phases if relevant
  if (method == 'first')
    cleaned <- cleaned[cleaned$phase == phases[1], ]
    
  if (method == 'last')
    cleaned <- cleaned[cleaned$phase == phases[length(phases)], ]

  # split by probe to start working
  probe_lists <- split(cleaned, cleaned$probe)

  if (missing(smoothing)) {
    bg_lists <- lapply(probe_lists, function(probe) {
      
      bg_lm <- lm(as.numeric(o2_delta) ~ as.numeric(phase_time), data = probe)

      # for the intercept to be 0
      bg_lm$coefficients[1] <- 0 

      # calculate the values for the relevant phase times
      output <- data.frame(phase_time = 1:max(probe$phase_time))

      output$o2_bg <- as.vector(
        stats::predict(bg_lm, output, type = "response", se.fit = FALSE)
      )

      # done
      return(output)
    })
  } else {
    bg_lists <- lapply(probe_lists, function(probe, smoothing) {
      # start by calculating the means by time point
      output <- stats::aggregate(x = probe$o2_delta, 
                                 by = list(probe$phase_time), 
                                 FUN = mean, na.rm = TRUE)
      colnames(output) <- c('phase_time', 'o2_bg')
      
      # then, if smoothing is being applied, do averages over time
      if (smoothing > 1) {
        # Bulk of the averages is automatically calculated
        x <- stats::filter(x = output$o2_bg, 
                           filter = rep(1/smoothing, smoothing), 
                           method = "convolution",
                           sides = 2)
        # but this will give you NAs on the head and tail of the vector.
        # While understandable, for our purpose we want to fill those gaps in.
        first_value <- head(which(!is.na(x)), 1)
        for (i in 1:(first_value - 1)) {
          # |  usable  range  |
          # |1 - - - i - - - -|2i
          usable_range <- 1:(2 * i - 1)
          x[i] <- mean(output$o2_bg[usable_range])
        }

        last_value <- tail(which(!is.na(x)), 1)
        for (i in (last_value + 1):length(x)) {
          # Determining usable ranges at the tail using length(x)
          #
          #  |     usable range     |
          #  | y  - - - i - - - l(x)|
          #  |    A    |-|     B    |
          #  
          #  A == B
          #  B = l(x) - i
          #  y = i - B

          #  Test, assume length(x) = 20 and i = 16
          #
          #  B = 20 - 16 = 4
          #  y = 16 - 4 = 12
          #  |    usable range    |
          #  |12 - - - 16 - - - 20| Correct.
          
          #  Test, assume length(x) = 20 and i = 17
          #
          #  B = 20 - 17 = 3
          #  y = 17 - 3 = 14
          #  |  usable range  |
          #  |14 - - 17 - - 20| Correct.

          B <- length(x) - i
          y <- i - B
          usable_range <- y:length(x)
          x[i] <- mean(output$o2_bg[usable_range])
        }

        output$o2_bg <- as.numeric(x)
      }
      return(output)
    }, smoothing = smoothing)
  }
  
  output <- as.data.frame(data.table::rbindlist(bg_lists, idcol = 'probe'))

  units(output$phase_time) <- units(cleaned$phase_time)
  units(output$o2_bg) <- units(cleaned$o2)

  input$bg <- output

  return(input)
}


#' Replace background readings for any given
#' probes using the values from another probe.
#' 
#' @param input The output of \code{\link{calc_bg}}.
#' @param replace A string of probe names to replace.
#' @param with The probe to use as a reference for replacement.
#' 
#' @return the updated bg dataframe
#' 
#' @export
#' 
replace_bg <- function(input, replace, with) {
  bg <- input$bg

  if (is.null(bg))
    stop("Could not find a bg object in the input")

  if (any(!(replace %in% bg$probe)))
    stop("Could not find some of the specified probes to replace in bg")

  if (length(with) != 1)
    stop("Please chose only one probe to use as replacement in 'with'.")
  
  if (!(with %in% bg$probe))
    stop("Could not find the replacement probe in bg")

  for (i in replace) {
    if (sum(bg$probe == i) > sum(bg$probe == with))
      stop("the cycle for the replacement probe is shorted than the cycle for the_probe to be replaced.")

    bg$o2_bg[bg$probe == i] <- bg$o2_bg[bg$probe == with][1:sum(bg$probe == i)]
    # that's some unholy sequential bracketing but I don't
    # recall why I did it so I am letting it stay :)
  }

  input$bg <- bg
  return(input)
}


#' Expand background readings linearly
#' 
#' Use this if your measurement phase is longer than your background and you
#' want ot expand your background linearly to accomodate the full duration
#' of the phase.
#' 
#' @param input an experiment list containing background readings.
#'  The output of \code{\link{calc_bg}}.
#' @param to The number of seconds to which to extend the background readings.
#' 
#' @return the updated input.
#' 
#' @export
#' 
extrapolate_bg <- function(input, to) {
  bg <- input$bg
  aux <- split(bg, bg$probe)

  recipient <- lapply(names(aux), function(P) {
    m <- lm(as.numeric(o2_bg) ~ as.numeric(phase_time),
            data = aux[[P]])
    output <- data.frame(probe = P,
                         phase_time = 1:to)
    units(output$phase_time) <- units(aux[[P]]$phase_time)
    output$o2_bg <- predict(m, output)
    output$o2_bg[1:nrow(aux[[P]])] <- aux[[P]]$o2_bg
    return(output)
  })

  output <- do.call(rbind, recipient)
  units(output$o2_bg) <- units(bg$o2_bg)

  input$bg <- output
  return(input)
}

#' Subtract background from slopes
#'
#'
#' @param input A list containing cleaned oxygen data.
#'  Obtained through \code{\link{process_experiment}}.
#' @param pre (optional) A data frame containing pre-test background readings.
#'  Obtained through \code{\link{calc_bg}}.
#' @param post (optional) A data frame containing post-test background readings.
#'  Obtained through \code{\link{calc_bg}}.
#' @param method  the name of the method used for background correction:
#' #' \itemize{
#' \item  "pre" - subtracts the pre-experiment background readings
#'          from the experiment oxygen consumption readings.
#' \item  "post" - subtracts the post-experiment background readings
#'          from the experiment oxygen consumption readings.
#' \item  "average" - averages the pre- and post-experiment background readings
#'          and subtracts them from the experiment oxygen consumption readings.
#' \item  "linear" - calculates a linear progression between pre- and
#'          post-experiment background readings and subtracts the respective
#'          background from the experiment oxygen consumption readings, using
#'          the phase number as an indicator of the position in the linear
#'          progression.
#' \item  "parallel" - subtracts the oxygen consumption readings of one probe
#'          (listed using ref_probe) from the remaining probes, matching both
#'          by cycle.
#' \item  "none" - does not perform oxygen consumption subtraction. 
#'          Not recommended for anything other than checking test data.
#' }
#' @param ref_probe  string: the name of an empty probe used 
#'  only for the method 'parallel'
#'
#' @return  A data frame containing data of metabolic rate measurements 
#'  corrected for background respiration.
#'
#' @references {Svendsen, M. B. S., Bushnell, P. G., & Steffensen, J. F. (2016). 
#' Design and setup of intermittent-flow respirometry system for aquatic 
#' organisms. Journal of Fish Biology, 88(1), 26-50.}
#'
#' @export
#'
subtract_bg <- function (input, pre, post,
                        method = c("pre", "post", "average",
                                   "linear", "parallel", "none"),
                        ref_probe){

  method <- match.arg(method)

  # bg object checks
  check <- method %in% c("pre", "average", "linear", "exponential")
  if (check && missing(pre)) {
    stop("method = '", method, "' but argument pre is missing.")
  }

  check <- method %in% c("post", "average", "linear", "exponential")
  if (check && missing(post)) {
    stop("method = '", method, "' but argument pre is missing.")
  }

  if (method %in% c("average", "linear", "exponential")) {
    split_pre <- split(pre$bg, pre$bg$probe)
    split_post <- split(post$bg, post$bg$probe)
    capture <- lapply(names(split_pre), function(probe) {
      if (nrow(split_pre[[probe]]) != nrow(split_post[[probe]])) {
        warning("The pre-bg and post-bg in probe ", probe, 
          " have different lengths! Truncating longer vector.", 
          immediate. = TRUE, call. = FALSE)
        
        if (nrow(split_pre[[probe]]) < nrow(split_post[[probe]])) {
          good_range <- 1:nrow(split_pre[[probe]])
          replacement <- split_post[[probe]][good_range, ]
          split_post[[probe]] <<- replacement
        } else {
          good_range <- 1:nrow(split_post[[probe]])
          replacement <- split_pre[[probe]][good_range, ]
          split_pre[[probe]] <<- replacement
        }
      }
    })
    pre$bg <- do.call(rbind, split_pre)
    post$bg <- do.call(rbind, split_post)

    if (units(pre$bg$o2_bg) != units(post$bg$o2_bg)) {
      stop("It seems the two background readings are not in the same unit! (",
        units(pre$bg$o2_bg), " != ", units(post$bg$o2_bg), ").")
    }
  }

  if (method == "parallel") {
    if (missing(ref_probe)) {
      stop("method = 'parallel' but argument ref_probe is missing.")
    }
    if (!any(ref_probe %in% input$cleaned$probe)) {
      stop("Could not find ref_probe in cleaned data.")
    }
  }

  # calculate bg objects to use downstream
  if (method == "pre") {
    my_bg <- pre$bg
    input$bg$pre <- pre
  }

  if (method == "post") {
    my_bg <- post$bg
    input$bg$post <- post
  }

  if (method == "average") {
    my_bg <- pre$bg
    my_bg$o2_bg <- (pre$bg$o2_bg + post$bg$o2_bg) / 2
    input$bg$pre <- pre
    input$bg$post <- post
  }

  if (method == "linear") {
    my_bg <- calc_linear_bg(input = input, pre = pre, post = post)
    input$bg$pre <- pre
    input$bg$post <- post
  }

  if (method == "parallel") {
    my_bg <- input$cleaned[input$cleaned$probe == ref_probe, ]
    my_bg$o2_bg <- my_bg$o2_delta
  }

  # make equivalent indexes 
  if (method %in% c("pre", "post", "average")) {
    input$cleaned$tmp_index <- paste(input$cleaned$probe, 
                                     input$cleaned$phase_time)
    
    my_bg$tmp_index <- paste(my_bg$probe, 
                             my_bg$phase_time)
  }

  if (method %in% c("linear", "exponential")) {
    input$cleaned$tmp_index <- paste(input$cleaned$cycle, 
                                     input$cleaned$probe,
                                     input$cleaned$phase_time)

    my_bg$tmp_index <- paste(my_bg$cycle, 
                             my_bg$probe, 
                             my_bg$phase_time)
  }

  if (method == "parallel") {
    my_bg$tmp_index <- paste(my_bg$cycle,
                             my_bg$phase_time)

    input$cleaned$tmp_index <- paste(input$cleaned$cycle,
                                     input$cleaned$phase_time)
  }

  # transfer bg readings
  if (method == "none") {
    input$cleaned$o2_bg <- 0
  } else {
    link <- match(input$cleaned$tmp_index, my_bg$tmp_index)
    input$cleaned$o2_bg <- my_bg$o2_bg[link]
    input$cleaned$tmp_index <- NULL
    
    if (any(is.na(input$cleaned$o2_bg)))
      warning("Some measurement phases are longer than the background.",
        " Discarding overextended points.", 
        immediate. = TRUE, call. = FALSE)

    input$cleaned <- input$cleaned[!is.na(input$cleaned$o2_bg), ]
  }

  #-----------------------------------------------------------------------------
  input$cleaned$o2_bg_delta <- input$cleaned$o2_bg - input$cleaned$o2_bg[1]
  input$cleaned$o2_cor <- input$cleaned$o2 - input$cleaned$o2_bg

  aux <- split(input$cleaned, paste0(input$cleaned$probe, input$cleaned$phase))
  aux <- lapply(aux, function(x) {
    x$o2_cordelta <- x$o2_cor - x$o2_cor[1]
    return(x)
  })

  output <- as.data.frame(data.table::rbindlist(aux))
  output <- transfer_attributes(input$cleaned, output)
  attributes(output)$correction_method <- method

  input$cleaned <- output
  
  return(input)
}



#' Internal function
#' 
#' Calculates the linear progression between pre- and post-background readings,
#' used by \code{\link{subtract_bg}}
#' 
#' @inheritParams subtract_bg
#' 
#' @keywords internal
#' 
calc_linear_bg <- function(input, pre, post) {
  cycles <- max(input$cleaned$cycle)

  my_bg <- lapply(unique(input$cleaned$probe), function(probe) {
    # cat(probe, "\n")
    pre_line <- pre$bg$o2_bg[pre$bg$probe == probe]
    post_line <- post$bg$o2_bg[post$bg$probe == probe]

    if (length(pre_line) != length(post_line)) {
      warning("The pre-bg and post-bg in probe ", probe, 
              " have different lengths! Truncating longer vector.", 
              immediate. = TRUE, call. = FALSE)
      
      if (length(pre_line) < length(post_line))
        post_line <- post_line[1:length(pre_line)]
      else
        pre_line <- pre_line[1:length(post_line)]
    }

    difference <- post_line - pre_line

    if (all(as.numeric(difference) == 0)) {
      stop("The pre-bg and post-bg are exactly the same!",
           " Cannot calculate linear progression.", call. = FALSE)
    }

    # linear increments.
    increments <- difference / (cycles - 1)

    my_list <- lapply(1:length(pre_line), function(i) {
      seq(from = pre_line[i], to = post_line[i], by = increments[i])
    })

    my_matrix <- as.data.frame(do.call(rbind, my_list))

    return(my_matrix)
  })
  names(my_bg) <- unique(input$cleaned$probe)

  my_bg_simplified <- lapply(names(my_bg), function(probe) {
    x <- suppressMessages(reshape2::melt(my_bg[[probe]]))
    colnames(x) <- c("cycle", "o2_bg")
    x$phase_time <- 1:nrow(my_bg[[probe]])
    x$cycle <- as.numeric(sub("V", "", x$cycle))
    x$probe <- probe

    units(x$phase_time) <- units(pre$bg$phase_time)
    units(x$o2_bg) <- units(pre$bg$o2_bg)
    return(x)
  })

  my_bg_df <- do.call(rbind, my_bg_simplified)
  return(my_bg_df)
}
