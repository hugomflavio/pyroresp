#' Compile a list of points following Killen et al. 2021
#'
#' For more information on these points, see Killen et al. 2021.
#'
#' @param input A resp object that went through \code{\link{process_mr}}.
#'
#' @references Killen et al. (2021). Guidelines for reporting
#'   methods to estimate metabolic rates by aquatic intermittent-flow
#'   respirometry. Journal of Experimental Biology, 224(18), jeb242522.
#'   <https://doi.org/10.1242/jeb.242522>
#'
#' @export
#'
#' @return a list of summary information
#'
killen_summary <- function(input) {
  # 04
  pt04 <- with(input$probe_info,
               (volume - conv_w_to_ml(mass)) / conv_w_to_ml(mass))
  # 15
  recipient <- lapply(input$pyro$source_data, function(i) {
    prb <- data.frame(probe = paste0(attributes(i)$device, attributes(i)$ch))
    cal <- t(as.data.frame(attributes(i)$cal_settings))
    return(cbind(prb, cal))
  })

  pt15 <- do.call(rbind, recipient)
  pt15 <- type.convert(pt15, as.is = TRUE)
  date_cols <- grepl("lastCal", colnames(pt15))
  if (any(date_cols)) {
    exp_day <- as.Date(min(input$pyro$compiled_data$date_time))
    for(i in which(date_cols)) {
      pt15[, i] <- as.Date(pt15[, i])
      n <- sub("lastCal", "", colnames(pt15)[i])
      link <- paste0("cal", n, "_to_exp_start")
      pt15[, link] <- difftime(exp_day, pt15[, i], units = "days")
    }
  }
  rownames(pt15) <- 1:nrow(pt15)
  
  # 17
  sem_aux <- sem(input$trimmed$temp, na.rm = TRUE)
  units(sem_aux) <- units(input$trimmed$temp)
  pt17 <- c(mean = mean(input$trimmed$temp, na.rm = TRUE),
            sem = sem_aux,
            min = min(input$trimmed$temp, na.rm = TRUE),
            max = max(input$trimmed$temp, na.rm = TRUE))

  # 22

  pt22aux <- aggregate(input$trimmed$airsat,
                       list(input$trimmed$phase),
                       min, na.rm = TRUE)
  pt22 <- c(mean = mean(pt22aux$x),
            lowest = min(pt22aux$x))

  # 30
  pt30 <- sapply(input$bg, function(i) {
            length(unique(i$trimmed$phase))
          })

  # 36
  pt36 <- gsub("_mr_g", "",
               colnames(input$smr)[grepl("_mr_g", colnames(input$smr))])

  # 39
  by_probe <- table(input$slopes$probe)
  by_probe <- as.data.frame(by_probe)
  colnames(by_probe) <- c("probe", "n_all")
  aux <- as.data.frame(table(input$good_slopes$probe))
  colnames(aux) <- c("probe", "n_good")
  by_probe <- merge(by_probe, aux, by = "probe", all = TRUE)
  by_probe$n_good[is.na(by_probe$n_good)] <- 0
  by_probe$n_discarded <- by_probe$n_all - by_probe$n_good
  by_probe$pct_discarded <- by_probe$n_discarded / by_probe$n_all * 100
  units(by_probe$pct_discarded) <- "percent"

  run_total <- as.data.frame(apply(by_probe[, -1], 2, sum, simplify = FALSE))
  run_total$pct_discarded <- run_total$n_discarded / run_total$n_all * 100
  units(run_total$pct_discarded) <- "percent"

  pt39 <- list(by_probe = by_probe,
               run_total = run_total)

  # bring it all together
  output <- list(
    pt01 = list(description = "mass of the animals",
                value = input$probe_info$mass),
    pt02 = list(description = "volume of empty respirometer",
                value = input$probe_info$volume),
    pt04 = list(description = "ratio of net resp volume",
                value = pt04),
    pt14 = list(description = "wait time (in number of data points)",
                value = attributes(input$trimmed)$wait),
    pt15 = list(description = "O2 probe calibration details",
                value = pt15),
    pt17 = list(description = "temperature details",
                value = round(pt17, 3)),
    pt22 = list(description = "lowest water oxygen levels reached",
                value = round(pt22, 3)),
    pt24 = list(description = paste0("number of animals in this run",
                                     " (assuming one animal per chamber)"),
                value = nrow(input$probe_info)),
    pt30 = list(description = "number of bg cycles used",
                value = pt30),
    pt31 = list(description = "method used for bg correction over time",
                value = attributes(input$slopes)$correction_method),
    pt35 = list(description = "duration over which MO2 was recorded",
                value = difftime(max(input$trimmed$date_time),
                                 min(input$trimmed$date_time))),
    pt36 = list(description = "method used for smr calculation",
                value = pt36),
    pt37 = list(description = "number of MO2 values used for smr calculation",
                value = table(input$good_slopes$probe)),
    pt39 = list(description = "r2 threshold for slopes",
                value = attributes(input$good_slopes)$r2_threshold),
    pt40 = list(description = "proportion of slopes below r2 threshold",
                value = pt39),
    pt51 = list(description = "How were oxygen uptake rates calculated?",
                value = paste0("Oxygen uptake rates were calculated using ",
                               "the R package pyroresp (version ",
                               packageVersion("pyroresp"),
                               "; Fl\u00e1vio, 2025).")),
    pt52 = list(description = paste0("Was the volume of the animal susbtracted",
                                     " from the respirometer volume when",
                                     " calculating oxygen uptake rates?"),
                value = "Yes")
    )
  return(output)
}
