#' Calculate SMR using various methods
#'
#' @param mr A data frame containing the metabolic rates.
#'  The output of \code{\link{calc_mr}}.
#' @param G A vector of the number of clusters to use for MLND SMR estimations.
#'  Defaults to 1:4. Set to NULL to skip MLND method.
#' @param q A vector of the quantiles to use for the quantile SMR estimations.
#'  Defaults to 0.2 and 0.25. Set to NULL to skip the quantile method.
#' @param p A vector of the percentages to use for low percent SMR estimations.
#'  Percent methods discard the five lowest measurements before calculating the
#'  average of the lowest x percentage. Defaults to 0.1. Set to NULL to skip
#'  the percent method.
#' @param n A vector of the number of measurements to use for the low SMR
#'  estimations. This method grabs the lowest n measurements and averages them.
#'  No values are discarded using this method. Defaults to 10. Set to NULL to
#'  skip the low n method.
#'
#' @importFrom mclust mclustBIC
#'
#' @references {Chabot, D., Steffensen, J. F., & Farrell, A. P. (2016).
#' The determination of standard metabolic rate in fishes.
#' Journal of Fish Biology, 88(1), 81-121.}
#'
calc_smr <- function(mr, G = 1:4, q = c(0.2, 0.25), p = 0.1, n = 10){
  by_probe <- split(mr, mr$probe)

  # flag to issue warning only once
  issue_low_n_warning <- TRUE

  recipient <- lapply(by_probe, function(the_probe) {

    if ("mass" %in% colnames(the_probe)) {
      output <- the_probe[1, c("probe", "id", "mass")]
      the_probe$input_mr <- the_probe$mr_g
    } else {
      output <- the_probe[1, c("probe", "id")]
      the_probe$input_mr <- the_probe$mr_abs
    }

    # MLND
    if (!is.null(G)) {
      calc_mlnd <- TRUE
      if (nrow(the_probe) < 10) {
        message("Skipping MLND method for probe ", the_probe$probe[1], 
                " (less than 10 values).")
        calc_mlnd <- FALSE
      }
      if (calc_mlnd && length(unique(the_probe$input_mr)) == 1) {
        message("Skipping SMR calculations for probe ", the_probe$probe[1], 
                " (only one unique value).")
        calc_mlnd <- FALSE
      }
      if (calc_mlnd) {
        the_mclust <- mclust::Mclust(the_probe$input_mr, G = G, verbose = FALSE)
        if (is.null(the_mclust)) {
          warning("mclust method failed for probe ", the_probe$probe[1],
                  ". Skipping MLND method. Low sample size?",
                  immediate. = TRUE, call. = FALSE)
          output$mlnd_mr <- NA
          units(output$mlnd_mr) <- units(the_probe$input_mr)
          output$mlnd_cv <- NA
        } else {
          cl <- the_mclust$classification
          clusters <- as.data.frame(table(cl))
          clusters$cl <- as.numeric(levels(clusters$cl))
          valid <- clusters$Freq >= 0.1 * sum(clusters$Freq)
          the_cl <- min(clusters$cl[valid])
          mlnd_index <- the_mclust$classification == the_cl
          mlnd <- the_mclust$parameters$mean[the_cl]

          output$mlnd_mr <- unname(mlnd)
          units(output$mlnd_mr) <- units(the_probe$input_mr)

          output$mlnd_cv <- sd(the_probe$input_mr[mlnd_index])/mlnd * 100
          attributes(output)$mlnd_phases <- the_probe$phase[mlnd_index]
        }
      }
    }

    # QUANTILES
    if (!is.null(q)) {
      for (i in q) {
        output[, paste0("q", i, "_mr")] <- quantile(the_probe$input_mr, i)
      }
    }

    if (!is.null(p)) {
      aux <- the_probe[order(the_probe$input_mr), ]
      if (nrow(aux) > 5) { # drop the lowest five values to avoid outlier effects
        aux <- aux[-(1:5), ]
      } else {
        aux <- aux[1, , drop = FALSE]
      }

      if (issue_low_n_warning && nrow(aux) < 35) {
        warning("The low percentage method does not work well",
                " with a low number of measurements (i.e. less than 35).",
                immediate. = TRUE, call. = FALSE)
        issue_low_n_warning <<- FALSE
        # we are inside an lapply here, hence the
        # assignment outside of the current environment
      }

      for (i in p) {
        aux <- aux[1:ceiling(nrow(aux) * i), ]
        output[, paste0("p", i, "_mr")] <- mean(aux$input_mr)
        attributes(output)[paste0("p", i, "_phase")] <- aux$phase
      }
    }

    if (!is.null(n)) {
      aux <- the_probe[order(the_probe$input_mr), ]
      for (i in n) {
        new_col <- paste0("low", i, "_mr")
        new_attr <- paste0("low", i, "_phases")
        if (nrow(aux) < n) {
          warning(paste0("Could not find enough measurements for probe ",
                    aux$probe[1], " to calculate low", i,
                    " smr method. Skipping."),
                  immediate. = TRUE, call. = FALSE)
          output[, new_col] <- NA
          units(output[, new_col]) <- units(the_probe$input_mr)
        } else {
          aux <- aux[1:i, ]
          output[, new_col] <- mean(aux$input_mr)
          attributes(output)[[new_attr]] <- aux$phase
        }
      }
    }

    return(output)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  return(output)
}
