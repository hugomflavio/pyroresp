#' Prepare imported measurements for further analyses
#' 
#' @param input the data frame containing imported oxygen measurements
#' @param wait integer: the number of first rows for each measurement phase (M) which should be reassigned to the wait phase (W). The parameter should be used when the wait phase (W) is absent (e.g. in 'Q-box Aqua' logger software) or not long enough to eliminate non-linear change in DO concentration over time from the measurement phase (M) after shutting off water supply from the ambient water source.
#' 
#' @return A data.frame containing measurements valid for further analyses
#' 
#' @export
#' 
clean.meas <- function(input, wait = 0, auto.cut.last = FALSE){

  input$Date <- as.Date(input$Date.Time)
  input$Real.Time <- chron::times(strftime(input$Date.Time, "%H:%M:%S"))

  # Removing Non-Measurement Data  
  input <- input[grepl("^M", input$Phase), ]
  
  # this is used just to ensure that the phases maintain their order,
  # even it they don't start at one or are not sorted at the start.
  phase.order <- as.numeric(gsub("[M]", "", unique(input$Phase)))
  phase.order <- phase.order - min(phase.order) + 1
  
  input$Phase <- factor(input$Phase, levels = unique(input$Phase)[phase.order])


  #the rest has to be done on a chamber by chamber basis.

  by.chamber <- split(input, input$Chamber.No)

  recipient <- lapply(names(by.chamber), function(the.chamber) {

    trimmed.db <- by.chamber[[the.chamber]]

    # Removing the final measurement Phase if necessary or forced (tail error)
    rows.per.phase <- table(trimmed.db$Phase)

    if (length(rows.per.phase) > 1)
      mean.rows.per.phase <- mean(rows.per.phase[-length(rows.per.phase)])
    else
      mean.rows.per.phase <- rows.per.phase

    if (tail(rows.per.phase, 1) < wait | auto.cut.last) {
      trimmed.db <- trimmed.db[trimmed.db$Phase != tail(levels(trimmed.db$Phase), 1), ]
      trimmed.db$Phase <- droplevels(trimmed.db$Phase)
    }


    # cut off first n rows from 'M' phase
    if(wait != 0){
      # the code below grabs the 1:nrow vector, breaks it out by phase, and then uses tail() with
      # a negative n to grab all numbers but the first 30 that show up for each phase.
      index <- unlist(tapply(1:nrow(trimmed.db), trimmed.db$Phase, tail, -(wait)), use.names = FALSE)
      trimmed.db <- trimmed.db[index, ]
    }

    # reset rownames
    row.names(trimmed.db) <- 1:nrow(trimmed.db)

    # and now, by phase, calculate the passing time and store the start and end points of that measurement
    aux <- split(trimmed.db, trimmed.db$Phase)

    aux <- aux[sapply(aux, nrow) > 0]

    aux <- lapply(aux, function(x) {
      x$Start.Meas <- x$Real.Time[1]
      x$End.Meas <- x$Real.Time[nrow(x)]
      x$Phase.Time <- as.numeric(difftime(x$Date.Time, x$Date.Time[1], units = 's'))

      return(x)
    })
    trimmed.db <- as.data.frame(data.table::rbindlist(aux))
 
    return(trimmed.db)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))

  return(output)
}
