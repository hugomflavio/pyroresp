#' Prepare imported measurements for further analyses
#' 
#' @param input the data frame containing imported oxygen measurements
#' @param meas.to.wait integer: the number of first rows for each measurement phase (M) which should be reassigned to the wait phase (W). The parameter should be used when the wait phase (W) is absent (e.g. in 'Q-box Aqua' logger software) or not long enough to eliminate non-linear change in DO concentration over time from the measurement phase (M) after shutting off water supply from the ambient water source.
#' 
#' @return A data.frame containing measurements valid for further analyses
#' 
#' @export
#' 
clean.meas <- function(input, meas.to.wait = 0){

  MR.data.all <- input  
  MR.data.all$Date <- as.Date(MR.data.all$Date.Time)
  MR.data.all$Real.Time <- chron::times(strftime(MR.data.all$Date.Time, "%H:%M:%S"))
  
  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  # Removing Non-Measurement Data
  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  
  MR.data.all <- MR.data.all[grepl("^M", MR.data.all$Phase), ]
  
  phase.order <- as.numeric(gsub("[M]", "", unique(MR.data.all$Phase)))
  
  MR.data.all$Phase <- factor(MR.data.all$Phase, levels = unique(MR.data.all$Phase)[phase.order])

  rm(phase.order)

  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  # Removing the final measurement Phase (tail error)
  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  rows.per.phase <- table(MR.data.all$Phase)

  if (length(rows.per.phase) > 1)
    mean.rows.per.phase <- mean(rows.per.phase[-length(rows.per.phase)])
  else
    mean.rows.per.phase <- rows.per.phase
  
  if (abs(tail(rows.per.phase, 1) - mean.rows.per.phase) > 1) {
    MR.data.all <- MR.data.all[MR.data.all$Phase != tail(levels(MR.data.all$Phase), 1), ]
    MR.data.all$Phase <- droplevels(MR.data.all$Phase)
  }

  rm(rows.per.phase, mean.rows.per.phase)

  # cut off first n rows from 'M' phase
  if(meas.to.wait != 0){
    idx <- unlist(tapply(1:nrow(MR.data.all), MR.data.all$Phase, tail, -(meas.to.wait)), use.names=FALSE)
    MR.data.all <- MR.data.all[idx, ]
    rm(idx)
  }

  row.names(MR.data.all) <- 1:nrow(MR.data.all)

  aux <- split(MR.data.all, MR.data.all$Phase)

  aux <- aux[sapply(aux, nrow) > 0]
  
  aux <- lapply(aux, function(x) {
    x$Start.Meas <- x$Real.Time[1]
    x$End.Meas <- x$Real.Time[nrow(x)]
    x$Phase.Time <- as.numeric(difftime(x$Date.Time, x$Date.Time[1], units = 's'))

    return(x)
  })
  MR.data.all <- as.data.frame(data.table::rbindlist(aux))
  
  rm(aux)

  return(MR.data.all)
}
