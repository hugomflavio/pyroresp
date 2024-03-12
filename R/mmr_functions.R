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
      index <- head(order(the_probe$mr_cor, decreasing = TRUE), 1)
      output <- the_probe[index, ]
    })

  mmr <- as.data.frame(data.table::rbindlist(recipient))
  
  # Keep only needed columns
  mmr <- mmr[, c("probe", "id", "mass", "date_time", "phase", "mr_cor")]

  return(mmr)
}


# #' Remove part of the mr trace
# #' 
# #' Useful when there are artifacts along the cycle that could disrupt fractioned
# #'  mr calculation. Refinds max mmr after removing the desired segment.
# #' 
# #' @param rolling_mmr the output of compile_rolling_mmr
# #' @param probe which probe to act upon
# #' @param from second from which to start exclusion
# #' @param to second where to end exclusion
# #' 
# #' @export
# #' 
# exclude_rolling_mmr_segment <- function(rolling_mmr, probe, from, to) {
#   stop("Don't run this until the code has been checked")
#   to_exclude <- rolling_mmr$detailed_mmr$probe == probe & 
#                 rolling_mmr$detailed_mmr$phase_time >= from & 
#                 rolling_mmr$detailed_mmr$phase_time <= to

#   rolling_mmr$detailed_mmr[to_exclude, 5:9] <- NA

#   aux <- rolling_mmr$detailed_mmr[rolling_mmr$detailed_mmr$probe == probe, ]
#   index <- which.max(aux$mr_cor)
  
#   prb <- rolling_mmr$max_mmr$probe == probe
#   rolling_mmr$max_mmr[prb, ] <- to_check[index, colnames(rolling_mmr$max_mmr)]

#   return(rolling_mmr)
# }



