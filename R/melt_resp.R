#' Dummy documentation
#' 
#' @export
#' 
melt_resp <- function(input, info.data) {
  base.cols <- c("Date.Time", "Date", "Real.Time", "Phase.Time", "Phase", "Start.Meas", "End.Meas")

  if (missing(info.data)) {
    ox.temp.cols <- colnames(input)[grepl("^Temp|^Ox", colnames(input))]
    ox.cols <- colnames(input)[grepl("^Ox", colnames(input))]
    temp.cols <- colnames(input)[grepl("^Temp", colnames(input))]
  } else {
    ox.temp.cols <- as.vector(outer(c('Temp', 'Ox'), sub('CH', '', info.data$Chamber.No), paste, sep = "."))
    ox.cols <- paste0("Ox.", sub('CH', '', info.data$Chamber.No))
    temp.cols <- paste0("Temp.", sub('CH', '', info.data$Chamber.No))
  }

  meas.data <- input[, c(base.cols, ox.temp.cols)]

  temp.df <- reshape2::melt(meas.data, id.vars = base.cols, measure.vars = ox.cols, variable.name = "Chamber.No", value.name = "O2.raw")
  temp.df$Chamber.No <- paste0("CH", sub("Ox.", "", temp.df$Chamber.No))
  temp.df$Temp <- reshape2::melt(meas.data, id.vars = NULL, measure.vars = temp.cols)$value

  if (!missing(info.data)) {
    link <- match(temp.df$Chamber.No, info.data$Chamber.No)
    temp.df <- cbind(temp.df, info.data[link, c("ID", "Mass", "Volume")])
  }
  
	by.chamber <- split(temp.df, temp.df$Chamber.No)
  
  recipient <- lapply(names(by.chamber), function(the.chamber, info.data) {

    trimmed.db <- by.chamber[[the.chamber]]

    if (!missing(info.data)) {
      cycles <- as.numeric(sub("F|M", "", as.character(trimmed.db$Phase)))
      first.meas <- info.data$First.meas[info.data$Chamber.No == the.chamber]
      trimmed.db <- trimmed.db[cycles >= first.meas, ]
    }

    by.phase <- split(trimmed.db, trimmed.db$Phase)
  	
    recipient <- lapply(by.phase, function(the.phase) {
  	  the.phase$O2.delta.raw <- the.phase$O2.raw - the.phase$O2.raw[1]
      return(the.phase)
  	})
    
    output <- data.table::rbindlist(recipient)
  }, info.data = info.data)

  output <- as.data.frame(data.table::rbindlist(recipient))

  output$Phase <- as.character(output$Phase)

  if (!missing(info.data))
    output$O2.unit <- info.data$O2.unit[1]

  return(output)
}
