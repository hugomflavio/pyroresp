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
    target <- ifelse ("animal_mass" %in% colnames(the_probe),
                      "mr_g", "mr_abs") {
    the_probe$mmr <- the_probe[, target]
    index <- order(the_probe$mmr, decreasing = TRUE)[1]
    output <- the_probe[index, ]
  })

  mmr <- as.data.frame(data.table::rbindlist(recipient))

  # Keep only needed columns
    mmr <- mmr[, c("probe", "date_time", "cycle", "mmr")]

  return(mmr)
}
