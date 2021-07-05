calc.slope <- function(input, length = Inf) {
  # The operation is done by phase and by chamber, so the dataset is broken twice below
  by.chamber <- split(input, input$Chamber.No) # first by chamber

  recipient <- lapply(by.chamber, function(the.chamber) {
    by.phase <- split(the.chamber, the.chamber$Phase) # now by phase

    recipient <- lapply(by.phase, function(the.phase) {
      trimmed.phase <- the.phase[the.phase$Phase.Time <= length, ]

      model.with.BR <- lm(O2.raw ~ Phase.Time, data = trimmed.phase)

      model.without.BR <- lm(O2.corrected ~ Phase.Time, data = trimmed.phase)

      output <- data.frame(Chamber.No = trimmed.phase$Chamber.No[1],
                           ID = trimmed.phase$ID[1],
                           Mass = trimmed.phase$Mass[1],
                           Volume = trimmed.phase$Volume[1],
                           Date.Time = trimmed.phase$Date.Time[nrow(trimmed.phase)],
                           Phase = trimmed.phase$Phase[1],
                           Temp = mean(trimmed.phase$Temp),
                           Slope.with.BR = coef(model.with.BR)[2],
                           Slope = coef(model.without.BR)[2],
                           SE = summary(model.without.BR)$coef[4],
                           R2 = summary(model.without.BR)$r.squared)

      return(output)
    })
    # start rebinding back to a dataframe
    output <- data.table::rbindlist(recipient)
    return(output)
  })

  output <- as.data.frame(data.table::rbindlist(recipient))
  return(output)
}


#' Extraction of Slope(s)
#'
#' The function extracts the slopes of the linear regression of corrected \eqn{O_{2}} concentration over time with defined parameters (see Arguments).
#'
#' @usage
#' extract.slope(clean.data,
#'               method = c("all", "min", "max",
#'                          "lower.tail", "upper.tail",
#'                          "calcSMR.mlnd", "calcSMR.quant",
#'                          "calcSMR.low10", "calcSMR.low10pc"),
#'               r2=0.95, length = 999999, n.slope = 1000,
#'               percent = 10, p = 0.25, G = 1:4)
#'
#' @param clean.data  a data frame obtained by using the function \code{\link{correct.meas}}
#' @param method  string: the method of extracting slopes:
#'  \itemize{
#'   \item 'all' all slopes
#'   \item 'min' extracts lowest absolute slopes, specify the number of extracted\cr slopes (parameter: n.slope)
#'   \item 'max' extracts highest absolute slopes, specify the number of extracted slopes (parameter: n.slope)
#'   \item 'lower.tail' extracts slopes from a lower tail of absolute slope distribution, specify percentage of a lower tail (parameter: percent)
#'   \item 'upper.tail' extracts slopes from an upper tail of absolute slope distribution, specify percentage of an upper tail (parameter: percent)
#'   \item 'calcSMR.mlnd' calculates the mean of the lowest normal distribution\cr (MLND) using the parameter G (see Appendix S1 in Chabot et al, 2016)
#'   \item 'calcSMR.quant' calculates quantile value of slope distribution using the parameter p (see Appendix S1 in Chabot et al, 2016)
#'   \item 'calcSMR.low10' calculates the mean of the 10 lowest absolute slopes (see Appendix S1 in Chabot et al, 2016)
#'   \item 'calcSMR.low10pc' calculates the mean of the lowest 10% of absolute slopes, after the 5 from 10 lowest have been removed as outliers (see Herrmann & Enders, 2000; Appendix S1 in Chabot et al, 2016)
#'   }
#' @param r2  numeric: minimal coefficient of determination (\eqn{r^{2}}) for extracted slopes. Coefficient of determination is used as a threshold of quality to be determined by the user (by default \eqn{r^{2}} = 0.95)
#' @param length  integer: length of a measurement period for slope calculations (in seconds; by default - full length)
#' @param n.slope  integer: the number of extracted slopes, only one slope is calculated for each measurement phase (used in the methods "min" and "max"; by default - all slopes)
#' @param percent integer: percentage of lower or upper tail (used in the methods "lower.tail" and "upper.tail", respectively; by default percent = 10)
#' @param p integer: p-value of quantile used in the method "calcSMR.quant" (by default p = 0.25)
#' @param G integer: G value is used in the method "calcSMR.mlnd" (by default G = 1:4)
#'
#' @return The function returns a data frame with the information about extracted slopes. The data frame is used in the functions \code{\link{QC.slope}} and \code{\link{calculate.MR}}.
#'
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom chron chron times
#' @importFrom grDevices dev.new
#' @importFrom stats coef lm predict.lm quantile time
#' @importFrom utils head read.table tail write.table
#'
#' @references {
#' \enumerate{
#'   \item Chabot, D., Steffensen, J. F., & Farrell, A. P. (2016). The determination of standard metabolic rate in fishes. Journal of Fish Biology, 88(1), 81-121.
#'   \item Herrmann, J. P., & Enders, E. C. (2000). Effect of body size on the standard metabolism of horse mackerel. Journal of Fish Biology, 57(3), 746-760.
#'   }
#' }
#'
#' @examples
#' # if the data have been already loaded to R,
#' # skip the first two lines of the code:
#' data(SMR.clean)
#' data(AMR.clean)
#'
#' SMR.slope <- extract.slope(SMR.clean,
#'                            method = "min",
#'                            n.slope = 3,
#'                            r2=0.95,
#'                            length = 1200)
#'
#' AMR.slope <- extract.slope(AMR.clean,
#'                            method = "all",
#'                            r2=0.95,
#'                            length = 300)
#'
#' @export

extract.slope <- function(slopes, method = c("all", "min", "max", "lower.tail", "upper.tail",
                                             "calcSMR.mlnd", "calcSMR.quant", "calcSMR.low10", "calcSMR.low10pc"),
                          r2 = 0.95, n.slope = 1000, percent = 10, p = 0.25, G = 1:4){

  method <- match.arg(method)

  slopes <- slopes[slopes$R2 >= r2, ]
  
  by.chamber <- split(slopes, slopes$Chamber.No)

  recipient <- lapply(by.chamber, function(the.chamber) {
    if(method == "all"){
      s <- order(the.chamber$Slope)
      output <- the.chamber[s,]
    }
    if (method == "min"){
      s <- head(order(the.chamber$Slope, decreasing = TRUE), n.slope)
      output <- the.chamber[s,]
    }
    if (method == "max"){
      s <- head(order(the.chamber$Slope, decreasing = FALSE), n.slope)
      output <- the.chamber[s,]
    }
    if (method == "lower.tail"){
      rate <- as.numeric(quantile(abs(the.chamber$Slope), percent/100))
      output <- the.chamber[abs(the.chamber$Slope) <= rate,]
    }
    if (method == "upper.tail"){
      rate <- as.numeric(quantile(the.chamber$Slope, percent/100))
      output <- the.chamber[the.chamber$Slope <= rate,]
    }
    else if (method == "calcSMR.mlnd"){
      the.Mclust <- Mclust(the.chamber$Slope, G = G)
      cl <- the.Mclust$classification
      cl2 <- as.data.frame(table(cl))
      cl2$cl <- as.numeric(levels(cl2$cl))
      valid <- cl2$Freq >= 0.1 * length(time)
      the.cl <- max(cl2$cl[valid])
      left.distr.R2 <- the.chamber$R2[the.Mclust$classification == the.cl]
      left.distr.Temp <- the.chamber$Temp[the.Mclust$classification == the.cl]
      left.distr.Slope.with.BR <- the.chamber$Slope.with.BR[the.Mclust$classification == the.cl]
      left.distr.SE <- the.chamber$SE[the.Mclust$classification == the.cl]
      left.distr <- the.chamber$Slope[the.Mclust$classification == the.cl]
      mlnd = the.Mclust$parameters$mean[the.cl]
      s <- as.numeric(mlnd)

      output<-the.chamber[1,]
      output$Date.Time <- NA
      output$Phase <- "M"
      output$Temp <- mean(left.distr.Temp)
      output$Slope.with.BR <- mean(left.distr.Slope.with.BR)
      output$Slope <- s
      output$SE <- mean(left.distr.SE)
      output$R2 <- mean(left.distr.R2)
    }
    if (method == "calcSMR.quant"){
      quant <- quantile(abs(the.chamber$Slope), p)
      rate <- as.numeric(-abs(quant))
      s <- which.min(abs(the.chamber$Slope - rate))
      output <- the.chamber[s,]
      output$Slope <- rate
    }
    if (method == "calcSMR.low10"){
      u <- sort(the.chamber$Slope, decreasing = TRUE)
      max.low10 <- max(u[1:10])
      min.low10 <- min(u[1:10])
      low10 <- mean(u[1:10])
      s <- as.numeric(low10)

      output <- the.chamber[1,]
      output$Date.Time <- NA
      output$Phase <- NA
      output$Temp <- mean(the.chamber$Temp[the.chamber$Slope >= min.low10])
      output$Slope.with.BR <- mean(the.chamber$Slope.with.BR[the.chamber$Slope >= min.low10])
      output$Slope <- s
      output$SE <- NA
      output$R2 <- mean(the.chamber$R2[the.chamber$Slope >= min.low10])
    }
    if (method == "calcSMR.low10pc"){
      u <- sort(the.chamber$Slope, decreasing = TRUE)
      max.low10pc <- max(u[6:10])
      min.low10pc <- min(u[6:10])
      low10pc <- mean(u[6:(5 + round(0.1 * (length(u) - 5)))])
      s <- as.numeric(low10pc)
      u.BR <- the.chamber$Slope.with.BR[the.chamber$Slope >= min.low10pc]
      low10pc.BR <- mean(u.BR[6:(5 + round(0.1*(length(u.BR)-5)))])
      s.BR <- as.numeric(low10pc.BR)
      output <- the.chamber[1,]
      output$Date.Time <- NA
      output$Phase <- NA
      output$Temp <- mean(the.chamber$Temp[the.chamber$Slope <= max.low10pc & the.chamber$Slope >= min.low10pc])
      output$Slope.with.BR <- s.BR
      output$Slope <- s
      output$SE <- NA
      output$R2 <- mean(the.chamber$R2[the.chamber$Slope <= max.low10pc & the.chamber$Slope >= min.low10pc ])
    }
    return(output)
  })

  good.slopes <- as.data.frame(data.table::rbindlist(recipient))

  good.slopes$Volume <- as.numeric(as.character(good.slopes$Volume))
  good.slopes$Mass <- as.numeric(as.character(good.slopes$Mass))
  good.slopes$DO.unit <- slopes$DO.unit[1]

  return(good.slopes)
}
