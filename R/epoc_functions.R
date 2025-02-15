#' Calculate area under the curve
#'
#' @param mmr post-chasing resp object
#' @param smr pre-chasing resp object
#' @param smr_col name of the SMR column to use (depends on the methods used
#'   to calculate the SMR in the SMR object)
#' @param smr_buffer Defaults to 0.1. Defines the zero for the auc calculations.
#'   0.1 means the AUC calculation will be 0 once the MO2 measured post-chase
#'   returns to 1.1 times the respective SMR of the animal.
#'
#' @return an updated mmr object, now containing 1) an updated mr object which
#'   shows the delta MO2, and 2) a new auc object. The auc object contains both
#'   the auc between each pair of points, and also the cumulative auc.
#'
calc_auc <- function(mmr, smr, smr_col, smr_buffer = 0.1) {
  if (is.null(mmr$mmr)) {
    stop("Couldn't find object 'mmr' inside mmr.",
         " Have you run process_experiment?", call. = FALSE)
  }
  if (is.null(smr$smr)) {
    stop("Couldn't find object 'smr' inside smr.",
         " Have you run process_experiment?", call. = FALSE)
  }
  if (length(smr_col) > 1) {
    stop("Select only one smr column", call. = FALSE)
  }
  if (!(smr_col %in% colnames(smr$smr))) {
    stop("Couldn't find column '", smr_col,
         "' in smr$smr.", call. = FALSE)
  }
  if (is.null(mmr$phases)) {
    stop("Couldn't find object 'phases' inside mmr.",
         " Is resp object corrupted?", call. = FALSE)
  }

  by_probe <- split(mmr$mr, mmr$mr$probe)

  # update mr table to include mr_over_smr and pct_over_smr
  by_probe <- lapply(by_probe, function(the_mr) {
    the_probe <- the_mr$probe[1]
    the_smr <- smr$smr[smr$smr$probe == the_probe, smr_col]
    the_phases <- mmr$phases[[which(names(mmr$phases) == the_probe)]]
    this_cycle <- the_phases$cycle == the_mr$cycle[1]
    this_phase <- the_phases$phase == "M"
    this_start <- the_phases$start[this_cycle & this_phase]

    the_mr$time_since_chase <- as_units(the_mr$date_time - this_start)
    units(the_mr$time_since_chase) <- "h"

    the_mr$mr_over_smr <- the_mr$mr_g - the_smr
    the_mr$fold_over_smr <- the_mr$mr_g / the_smr
    return(the_mr)
   })
  # splitting the lapply in two allows saving the updated mr tables

  # then use the new cols to calculate auc.
  recipient <- lapply(by_probe, function(the_mr) {
    the_probe <- the_mr$probe[1]
    the_smr <- smr$smr[smr$smr$probe == the_probe, smr_col]

    the_auc <- auc(y = the_mr$mr_over_smr,
                   x = the_mr$time_since_chase,
                   zero = the_smr * smr_buffer)
    colnames(the_auc)[1:2] <- c("time_delta", "mr_delta")
    the_auc$probe <- the_probe
    return(the_auc)
  })

  mmr$mr <- do.call(rbind, by_probe)
  mmr$auc <- do.call(rbind, recipient)

  # simplify row names
  rownames(mmr$mr) <- 1:nrow(mmr$mr)

  return(mmr)
}

#' Extracts the cumulative AUC
#'
#' Until the stepwise AUC reaches 0, or max_duration is reached.
#'
#' @param input a resp object with an "auc" sub-object.
#'   The output of \code{\link{calc_auc}}
#' @param max_duration maximum time allowed for the auc to reach 0 (in hours).
#'   if AUC does not reach 0 before this threshold, the cumulative AUC for the
#'   whole duration is used.
#'
#' @return the input object, with an 'epoc' sub-object
#'
extract_epoc <- function(input, max_duration = Inf) {
  if (is.null(input$auc)) {
    stop("Couldn't find object 'auc' inside input.",
         " Have you run calc_auc?", call. = FALSE)
  }
  units(max_duration) <- units(input$auc$time_delta)

  by_probe <- split(input$auc, input$auc$probe)

  recipient <- lapply(by_probe, function(x) {
    x <- x[x$time_delta <= max_duration, ]
    epoc_ended <- min(which(as.numeric(x$mr_delta) < 0), nrow(x))
    output <- data.frame(probe = x$probe[1],
                         epoc = x$cumauc[epoc_ended],
                         epoc_duration = x$time_delta[epoc_ended])
    this_probe <- input$mr$probe == x$probe[1]
    this_time <- input$mr$time_since_chase == x$time_delta[epoc_ended]
    output$fold_over_smr <- input$mr$fold_over_smr[this_probe & this_time]
    return(output)
  })
  the_epoc <- do.call(rbind, recipient)

  input$epoc <- the_epoc
  return(input)
}
