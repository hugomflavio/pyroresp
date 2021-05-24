melt_resp <- function(input, info.data) {
  base.cols <- c("Date.Time", "Date", "Real.Time", "Phase.Time", "Phase", "Start.Meas", "End.Meas")

  meas.data <- input[, c(base.cols, colnames(input)[grepl("^Temp|^Ox", colnames(input))])]
  row.names(meas.data) <- NULL

  temp.df <- reshape2::melt(meas.data, id.vars = base.cols, measure.vars = colnames(meas.data)[grepl("^Ox", colnames(meas.data))], variable.name = "Chamber.No", value.name = "O2.raw")
  temp.df$Chamber.No <- paste0("CH", sub("Ox.", "", temp.df$Chamber.No))
  temp.df$Temp <- reshape2::melt(meas.data, id.vars = NULL, measure.vars = colnames(meas.data)[grepl("^Temp", colnames(meas.data))])$value

  if (!missing(info.data))
    temp.df <- cbind(temp.df, info.data[as.numeric(sub("CH", "", temp.df$Chamber.No)), c("ID", "Mass", "Volume")])
  
	aux <- split(temp.df, paste0(temp.df$Chamber.No, temp.df$Phase))
	aux <- lapply(aux, function(x) {
	  x$O2.delta.raw <- x$O2.raw - x$O2.raw[1]
	return(x)
	})

  output <- as.data.frame(data.table::rbindlist(aux))

  if (!missing(info.data))
    output$O2.unit <- info.data$O2.unit[1]

  return(output)
}
