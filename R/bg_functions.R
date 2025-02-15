#' Calculate background values
#'
#' @param input The output of \code{\link{process_experiment}}.
#' @param method How should multiple background cycles for a single probe
#'  be handled? One of:
#'    'mean'  : All the measured cycles are used for the calculations
#'    'first' : Only the first cycle is used for the calculations
#'    'last'  : Only the last cycle is used for the calculations
#'
#' @return The input list, with a new "bg" object.
#'
#' @export
#'
calc_bg <- function(input, method = c('mean', 'first', 'last')){

  method <- match.arg(method)

  # extract some variables we'll be using for convenience
  trimmed <- input$trimmed
  cycles <- unique(trimmed$cycle)

  # discard unnecessary cycles if relevant
  if (method == 'first') {
    trimmed <- trimmed[trimmed$cycle == cycles[1], ]
  }

  if (method == 'last') {
    trimmed <- trimmed[trimmed$cycle == cycles[length(cycles)], ]
  }

  # split by probe to start working
  probe_lists <- split(trimmed, trimmed$probe)

  bg_lists <- lapply(probe_lists, function(probe) {

    aux <- with(probe,
                aggregate(as.numeric(o2_delta),
                          list(phase_time = as.numeric(phase_time)),
                          mean))
    colnames(aux)[2] <- "o2_delta"
    bg_lm <- lm(o2_delta ~ phase_time, data = aux)

    # force the intercept to be 0
    bg_lm$coefficients[1] <- 0

    output <- data.frame(slope = bg_lm$coefficients[2],
                         R2 = summary(bg_lm)$adj.r.squared)
    # done
    return(output)
  })

  output <- as.data.frame(data.table::rbindlist(bg_lists, idcol = 'probe'))

  units(output$slope) <- units(trimmed$o2_delta[1]/trimmed$phase_time[1])

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
#' @return the updated input object
#'
#' @export
#'
replace_bg <- function(input, replace, with) {
  bg <- input$bg

  if (is.null(bg)) {
    stop("Could not find a bg object in the input")
  }

  if (any(!(replace %in% bg$probe))) {
    stop("Could not find some of the specified probes to replace in bg")
  }

  if (length(with) != 1) {
    stop("Please chose only one probe to use as replacement in 'with'.")
  }

  if (!(with %in% bg$probe)) {
    stop("Could not find the replacement probe in bg")
  }

  # add a space to note that it has been replaced
  if (!("note" %in% colnames(bg))) {
    bg$note <- ""
  }

  for (i in replace) {
    if (sum(bg$probe == i) > sum(bg$probe == with)) {
      stop("the cycle for the replacement probe is shorted than the cycle",
           " for the_probe to be replaced.")
    }

    bg$slope[bg$probe == i] <- bg$slope[bg$probe == with]
    bg$R2[bg$probe == i] <- bg$R2[bg$probe == with]
    bg$note[bg$probe == i] <- paste("Copied from ", with)
  }

  input$bg <- bg
  return(input)
}

#' Subtract background from slopes
#'
#' @param input A list containing trimmed oxygen data.
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
  if (check) {
    if (missing(pre)) {
      stop("method = '", method, "' but argument pre is missing.")
    }
    probe_check <- unique(input$trimmed$probe) %in% unique(pre$bg$probe)
    if (any(!probe_check)) {
      stop("Could not find probe(s) ",
           paste(unique(input$trimmed$probe)[!probe_check], collapse = ", "),
           " in the pre background data.")
    }
  }

  check <- method %in% c("post", "average", "linear", "exponential")
  if (check) {
    if (missing(post)) {
      stop("method = '", method, "' but argument post is missing.")
    }
    probe_check <- unique(input$trimmed$probe) %in% unique(post$bg$probe)
    if (any(!probe_check)) {
      stop("Could not find probe(s) ",
           paste(unique(input$trimmed$probe)[!probe_check], collapse = ", "),
           " in the post background data.")
    }
  }

  if (method %in% c("average", "linear", "exponential")) {
    if (units(pre$bg$slope) != units(post$bg$slope)) {
      stop("It seems the two background readings are not in the same unit! (",
           units(pre$bg$slope), " != ", units(post$bg$slope), ").")
    }
  }

  if (method == "parallel") {
    if (missing(ref_probe)) {
      stop("method = 'parallel' but argument ref_probe is missing.")
    }
    if (!any(ref_probe %in% input$trimmed$probe)) {
      stop("Could not find ref_probe in trimmed data.")
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
    my_bg$slope <- (pre$bg$slope + post$bg$slope) / 2
    input$bg$pre <- pre
    input$bg$post <- post
  }

  if (method == "linear") {
    my_bg <- calc_linear_bg(input = input, pre = pre, post = post)
    input$bg$pre <- pre
    input$bg$post <- post
  }

  if (method == "parallel") {
    stop("this method needs to be updated for new slope bg output.")
    # my_bg <- input$slopes[input$slopes$probe == ref_probe, ]
    # my_bg$o2_bg_delta <- my_bg$o2_delta
  }

  # make equivalent indexes
  if (method %in% c("pre", "post", "average")) {
    input$slopes$tmp_index <- input$slopes$probe
    my_bg$tmp_index <- my_bg$probe
  }

  if (method %in% c("linear", "exponential")) {
    input$slopes$tmp_index <- paste(input$slopes$probe,
                                     input$slopes$cycle)

    my_bg$tmp_index <- paste(my_bg$probe,
                             my_bg$cycle)
  }

  if (method == "parallel") {
    stop("this method needs to be updated for new slope bg output.")
    # my_bg$tmp_index <- paste(my_bg$cycle,
    #                          my_bg$phase_time)

    # input$slopes$tmp_index <- paste(input$slopes$cycle,
    #                                  input$slopes$phase_time)
  }

  # transfer bg readings
  if (method == "none") {
    input$slopes$slope_bg <- 0
    units(input$slopes$slope_bg) <- units(input$slopes$slope)
  } else {
    link <- match(input$slopes$tmp_index, my_bg$tmp_index)
    input$slopes$slope_bg <- my_bg$slope[link]
    input$slopes$tmp_index <- NULL
  }
  #-----------------------------------------------------------------------------
  input$slopes$slope_cor <- with(input$slopes, slope - slope_bg)
  input$slopes$bg_pct_of_cor <- with(input$slopes, slope_bg / slope_cor)
  # units is now "1"; changing to percent automatically multiplies by 100
  units(input$slopes$bg_pct_of_cor) <- "percent"

  attributes(input$slopes)$correction_method <- method

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
  cycles <- max(input$slopes$cycle)

  my_bg <- lapply(unique(input$trimmed$probe), function(probe) {
    # cat(probe, "\n")
    pre_slope <- pre$bg$slope[pre$bg$probe == probe]
    post_slope <- post$bg$slope[post$bg$probe == probe]

    slope_diff <- post_slope - pre_slope

    if (as.numeric(slope_diff) == 0) {
      stop("The pre-bg and post-bg are exactly the same!",
           " Cannot calculate linear progression.", call. = FALSE)
    }

    # linear slope_incr.
    slope_incr <- slope_diff / (cycles - 1)

    slope_bg <- seq(from = pre_slope,
                   to = post_slope,
                   by = slope_incr)

    output <- data.frame(probe = probe,
                         cycle = 1:cycles,
                         slope_bg = slope_bg)

    return(output)
  })
  output <- do.call(rbind, my_bg)
  return(output)
}
