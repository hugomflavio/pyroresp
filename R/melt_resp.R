#' Dummy documentation
#' 
#' @export
#' 
melt_resp <- function(input, info.data, O2.unit) {

  pre_output <- data.table::data.table(Date.Time = rep(input$Date.Time, each = sum(grepl("Ox", colnames(input)))),
                       Probe = rep(sub("Phase.", "", colnames(input)[grep("Phase.", colnames(input))]), length.out = nrow(input) * sum(grepl("Phase", colnames(input)))),
                       Phase = as.vector(t(input[, grepl("Phase", colnames(input))])),
                       Pressure = as.vector(t(input[, grepl("Pressure", colnames(input))])),
                       Temp = as.vector(t(input[, grepl("Temp", colnames(input))])),
                       O2.raw = as.vector(t(input[, grepl("Ox", colnames(input))])))

  if (any(grepl("pH", colnames(input))))
    pre_output$pH <- as.vector(t(input[, grepl("pH", colnames(input))]))

  pre_output$Cycle <- as.numeric(sub("F|M", "", as.character(pre_output$Phase)))

  # if fish information is provided
  if (!missing(info.data)) {
    # include fish information
    link <- match(pre_output$Probe, info.data$Probe)
    pre_output <- cbind(pre_output, info.data[link, c("ID", "Mass", "Volume")])

    # and trim away the cycles that happen before the first measurement

    by.chamber <- split(pre_output, pre_output$Probe)
    
    recipient <- lapply(names(by.chamber), function(the.chamber, info.data) {

      trimmed.db <- by.chamber[[the.chamber]]

      first.meas <- info.data$First.meas[info.data$Probe == the.chamber]
      trimmed.db <- trimmed.db[trimmed.db$Cycle >= first.meas, ]
  
      return(trimmed.db)
    }, info.data = info.data)

    output <- as.data.frame(data.table::rbindlist(recipient))
  } else {
    output <- as.data.frame(pre_output)
  }

  output$Phase <- as.character(output$Phase)

  colnames(output)[colnames(output) == 'O2.raw'] <- paste0('O2.', O2.unit)

  return(output)
}



