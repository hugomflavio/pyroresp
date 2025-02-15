#' Find cycle with highest mr for each probe
#'
#' @param mr The metabolic rate data frame. The output of \code{\link{calc_mr}}.
#'
#' @return A data frame with the mmr cycle and and MO2 for each probe.
#'
#' @export
#'
extract_mmr <- function(mr){
  by_probe <- split(mr, mr$probe)

  recipient <- lapply(by_probe, function(the_probe) {
      index <- head(order(the_probe$mr_g, decreasing = TRUE), 1)
      output <- the_probe[index, ]
    })

  mmr <- as.data.frame(data.table::rbindlist(recipient))

  # Keep only needed columns
  mmr <- mmr[, c("probe", "id", "mass", "date_time", "cycle", "mr_g")]

  return(mmr)
}
