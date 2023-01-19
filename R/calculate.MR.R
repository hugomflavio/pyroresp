#' Calculation of Metabolic Rate
#'
#' The function is used to calculate and plot background respiration, absolute and mass-specific metabolic rates.
#'
#' @usage
#' calculate.MR(slope.data, density = 1000,
#'              plot.BR = TRUE,
#'              plot.MR.abs = TRUE,
#'              plot.MR.mass = TRUE)
#'
#' @param slope.data  a data frame obtained by using the function \code{\link{extract.slope}}
#' @param density  numeric: the density of an animal body (\eqn{kg/m^{3}})
#' @param plot.BR  logical: if TRUE, the graph of background respiration rate is plotted
#' @param plot.MR.abs  logical: if TRUE, the graph of absolute metabolic rate is plotted
#' @param plot.MR.mass  logical: if TRUE, the graph of mass-specific metabolic rate is plotted
#'
#' @return The function returns a data frame with calculated background respiration, absolute and mass-specific metabolic rates. The data frame is used in the function \code{\link{export.MR}}.
#'
#' @importFrom lattice xyplot
#' @importFrom grDevices dev.new
#' @importFrom graphics abline legend par plot
#' @importFrom stats coef lm predict.lm
#' @importFrom utils head read.table tail write.table
#'
#' @examples
#' # if the data have been already loaded to R,
#' # skip the first two lines of the code:
#' data(SMR.slope)
#' data(AMR.slope)
#'
#' SMR <- calculate.MR(SMR.slope,
#'                     density = 1000,
#'                     plot.BR = TRUE,
#'                     plot.MR.abs = TRUE,
#'                     plot.MR.mass = TRUE)
#'
#' AMR <- calculate.MR(AMR.slope,
#'                     density = 1000,
#'                     plot.BR = TRUE,
#'                     plot.MR.abs = TRUE,
#'                     plot.MR.mass = TRUE)
#'
#' @export

calculate.MR  <- function(slope.data, density = 1000, plot.BR = FALSE,
                          plot.MR.abs = FALSE, plot.MR.mass = FALSE){

  BW = slope.data$Mass/1000 # convert Mass from 'g' to 'kg'
  V = slope.data$Volume/1000 - (BW/(density/1000)) #convert Volume from 'mL' to 'L' and Density from 'kg/M^3' to 'kg/L'

  slope.data$MR.abs.with.BR = -(slope.data$Slope.with.BR*V*3600) #3600 sec = 1 hour

  slope.data$BR = (slope.data$Slope.with.BR - slope.data$Slope)/slope.data$Slope.with.BR*100 # in %
  slope.data$MR.abs = -(slope.data$Slope*V*3600)
  slope.data$MR.mass = -(slope.data$Slope*V*3600/BW)

  a <- xyplot(BR~Temp|ID, data=slope.data, as.table = TRUE,
              xlab = bquote("Temperature (" ~ C^o ~ ")"), ylab = "Background respiration (%)",
              main = "Percentage rate of background respiration")

  b <- xyplot(MR.abs~Temp|ID, data=slope.data, as.table = TRUE,
              xlab = bquote("Temperature (" ~ C^o ~ ")"), ylab = bquote("Absolute MR (" ~ .(slope.data$DO.unit[1]) ~ h^-1 ~ ")"),
              main = "Absolute metabolic rate")

  d <- xyplot(MR.mass~Temp|ID, data=slope.data, as.table = TRUE,
              xlab = bquote("Temperature (" ~ C^o ~ ")"), ylab = bquote("Mass-specific MR (" ~ .(slope.data$DO.unit[1]) ~ kg^-1 ~ h^-1 ~ ")"),
              main = "Mass-specific metabolic rate")

  if (plot.BR){
    par(mfrow = c(2, 1), ask = TRUE)
    print(a)
  }

  if (plot.MR.abs){
    par(mfrow = c(2, 1), ask = TRUE)
    print(b)
  }

  if (plot.MR.mass){
    par(mfrow = c(2, 1), ask = TRUE)
    print(d)
  }

  # a <- slope.data$DO.unit[1]
  # MR.data <- slope.data[,-12]
  # MR.data$DO.unit <- a[1]
  MR.data <- slope.data

  attributes(MR.data)$selection.method <- attributes(slope.data)$selection.method
  return(MR.data)
}




#' fraction MR
#' @export
fraction_mr <- function(input, chamber, phase, smoothing) {
  x <- input$corrected[input$corrected$Phase == phase & input$corrected$Chamber.No == chamber, ]
  head(x)

  recipient <- lapply(smoothing:nrow(x), function(i) {
    aux_m <- lm(O2.delta.corrected ~ Phase.Time,
          data = x[(i - smoothing + 1):i, ])
    output <- data.frame(
      Chamber.No = chamber,
      Phase = phase,
      Sec = i,
      Smoothing = smoothing,
      Slope = coef(aux_m)[2],
      R2 = summary(aux_m)$r.squared)
    return(output)
  })

  output <- data.table::rbindlist(recipient)

  BW <- x$Mass[1]/1000
  V <- x$Volume[1]/1000 - BW

  output$MR.abs <- -(output$Slope * V * 3600)
  output$MR.mass <- -(output$Slope * V * 3600 / BW)

  return(output)
}


#' dummy doc
#' 
#' @export
#' 
calculate.bg <- function(input, O2_col, method = c('mean', 'first', 'last'), force.linear = TRUE, smoothing = 30){

  method <- match.arg(method)

  n.phases <- unique(input$Phase)

  if (method == 'first')
    input <- input[input$Phase == n.phases[1], ]
    
  if (method == 'last')
    input <- input[input$Phase == n.phases[length(n.phases)], ]

  chamber.lists <- split(input, input$Chamber.No)

  bg.lists <- lapply(names(chamber.lists), function(chamber) {
    # cat(chamber, '\n')
    sub.data <- input[input$Chamber.No == chamber, ]

    if (force.linear) {

      bg.lm <- eval(parse(text = paste("lm(", O2_col, "~ Phase.Time, data = sub.data)")))
      bg.lm$coefficients[1] <- 0 

      output <- data.frame(Phase.Time = 1:max(sub.data$Phase.Time))
      output$O2.background <- as.vector(predict(bg.lm, output, type = "response", se.fit = FALSE))
    }
    else {
      output <- aggregate(sub.data[, O2_col], by = list(sub.data$Phase.Time), mean)
      colnames(output) <- c('Phase.Time', 'O2.background')
      
      if (smoothing > 1) {
        x <- stats::filter(output$O2.background, rep(1/smoothing, smoothing), sides = 2)

        if (smoothing%%2 == 0) {
          for (i in 1:(smoothing/2-1)) {
            r <- round(i/2)
            x[i] <- mean(output$O2.background[(i-r):(i+r)])
          }

          last.value <- length(x) - round(smoothing/2)
          smooth.values <- (last.value - round(smoothing/2)):last.value
          smooth.slope <- mean(x[smooth.values-1]-x[smooth.values])

          for (i in (length(x)-(smoothing/2)):length(x)) {
            x[i] <- x[i-1] - smooth.slope
          }
        }
        output$O2.background <- x
      }
    }  
    return(output)
  })
  names(bg.lists) <- names(chamber.lists)

  output <- as.data.frame(data.table::rbindlist(bg.lists, idcol = 'Chamber.No'))

  return(output)
}


#' Override background correction factor for any given
#' chambers using the values from another chamber.
#' 
#' @export
#' 
override_bg <- function(bg, replace, with) {
  if (!is.data.frame(bg))
    stop("bg must be a data.frame")

  if (all(!grepl("Chamber.No", colnames(bg))))
    stop("bg must contain a 'Chamber.No' column.")
  
  if (all(!grepl("O2.background", colnames(bg))))
    stop("bg must contain a 'O2.background' column.")
  
  if (any(!(replace %in% bg$Chamber.No)))
    stop("Could not find some of the specified chambers to replace in bg")

  if (length(with) != 1)
    stop("Please chose only one chamber to use as replacement in 'with'.")
  
  if (!(with %in% bg$Chamber.No))
    stop("Could not find the replacement chamber in bg")

  for (i in replace) {
    if (sum(bg$Chamber.No == i) > sum(bg$Chamber.No == with))
      stop("the cycle for the replacement chamber is shorted than the cycle for the chamber to be replaced.")

    bg$O2.background[bg$Chamber.No == i] <- bg$O2.background[bg$Chamber.No == with][1:sum(bg$Chamber.No == i)]
    # the additional subset at the end ensures 
  }

  return(bg)
}



#' Dummy documentation
#' 
#' @export
#' 
calc_delta <- function(input, O2_col) {

  input[, paste0('O2.delta', sub('O2', '', O2_col))] <- NULL
  
  by.chamber <- split(input, input$Chamber.No)
  
  recipient <- lapply(names(by.chamber), function(the.chamber, info.data) {

    trimmed.db <- by.chamber[[the.chamber]]

    by.phase <- split(trimmed.db, trimmed.db$Phase)
    
    recipient <- lapply(by.phase, function(the.phase) {
      the.phase$O2.delta <- the.phase[, O2_col] - the.phase[1, O2_col]
      return(the.phase)
    })
    
    output <- data.table::rbindlist(recipient)

  })

  output <- as.data.frame(data.table::rbindlist(recipient))
  colnames(output)[ncol(output)] <- paste0('O2.delta', sub('O2', '', O2_col))

  return(output)
}