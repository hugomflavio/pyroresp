calc.mean.slope <- function(x) {
  by.chamber <- split(x, x$Chamber.No)

  recipient <- lapply(by.chamber, function(the.chamber) {

    by.phase <- split(the.chamber, the.chamber$Phase)

    # Some phases may have an extra second, so we grab the smallest phase, to avoid NAs
    rows.per.phase <- min(sapply(by.phase, function(x) max(x$Phase.Time)))

    # extract times to confirm they match
    Phase.Time.list <- lapply(by.phase, function(x) x$Phase.Time[1:rows.per.phase])
    Phase.Time.df <- as.data.frame(Phase.Time.list)
    Phase.Time.df$Check <- apply(Phase.Time.df, 1, function(r) length(unique(r)) == 1)
    if (!all(Phase.Time.df$Check))
      stop('There is a mismatch on the phase times of the pre/post data.', call. = FALSE)

    # now extract the O2 values and calculate the mean
    O2.delta.raw.list <- lapply(by.phase, function(x) x$O2.delta.raw[1:rows.per.phase])
    O2.delta.raw.df <- as.data.frame(O2.delta.raw.list)
    O2.delta.raw.df$Mean <- apply(O2.delta.raw.df, 1, mean)

    # bind and return
    output <- cbind(Phase.Time = Phase.Time.df[, 1], O2.delta.raw.df)
    return(output)
  })

  # now combine all chambers as in the original input
  output <- as.data.frame(data.table::rbindlist(recipient, idcol = 'Chamber.No'))

  return(output)
}


#' Correction of Metabolic Rate Measurements
#'
#' The function is used to correct metabolic rate measurements for background respiration. To this end, oxygen consumption is estimated as the slope of the linear regression of measured \eqn{O_{2}} concentration over time, and is extracted for background respiration test and for each measurement phase. The correction is based on subtraction of oxygen consumption obtained during background respiration test from oxygen consumption obtained during metabolic rate measurements.
#'
#' @usage
#' correct.meas(pre.data, post.data, meas.data,
#'              method = c("pre.test", "post.test", "average",
#'                         "linear", "exponential", "parallel"),
#'              empty.chamber = c("CH1", "CH2", "CH3", "CH4",
#'                                "CH5", "CH6", "CH7", "CH8"))
#'
#' @param pre.data  a data frame obtained by using the function \code{\link{import.test}} for a blank test before actual metabolic rate measurements
#' @param post.data  a data frame obtained by using the function \code{\link{import.test}} for a blank test after actual metabolic rate measurements
#' @param meas.data  a data frame obtained by using the function \code{\link{import.meas}} for actual metabolic rate measurements
#' @param method  string: the name of the method used for background respiration correction:
#' \itemize{
#' \item  "pre.test" - subtracts oxygen consumption of pre.data from oxygen consumptions of meas.data
#' \item  "post.test" - subtracts oxygen consumption of post.data from oxygen consumptions of meas.data
#' \item  "average" - subtracts an averaged oxygen consumption of pre.data and\cr post.data from oxygen consumptions of meas.data
#' \item  "linear" - subtracts a vector of progressively changing microbial consumptions from oxygen consumptions of meas.data. The values of oxygen consumption are linearly predicted from two reference points: oxygen consumption of pre.data and oxygen consumption of post.data.
#' \item  "exponential" - subtracts a vector of progressively changing microbial consumptions from oxygen consumptions of meas.data. The values of oxygen consumption are exponentially predicted from two reference points: oxygen consumption of pre.data and oxygen consumption of post.data.
#' \item  "parallel" - subtracts oxygen consumption in an empty chamber from oxygen consumptions of meas.data for each chamber
#' \item  "none" - does not perform oxygen consumption subtraction. Not recommended for anything other than checking test data.
#' }
#' @param empty.chamber  string: the name of an empty chamber used only for the method 'parallel'
#'
#' @importFrom chron chron dates times
#' @importFrom grDevices dev.new
#' @importFrom graphics abline legend par plot
#' @importFrom stats coef lm predict.lm
#' @importFrom utils head read.table tail write.table
#'
#' @return  The function returns a data frame containing data of metabolic rate measurements corrected for background respiration. The data frame is used in the functions \code{\link{QC.meas}}, \code{\link{QC.activity}},\cr \code{\link{extract.slope}} and \code{\link{QC.slope}}.
#'
#' @examples
#' # if the data have been already loaded to R,
#' # skip the first five lines of the code:
#' data(info)
#' data(pre)
#' data(post)
#' data(AMR.raw)
#' \dontrun{
#' data(SMR.raw)
#' SMR.clean <- correct.meas(pre.data = pre,
#'                           meas.data = SMR.raw,
#'                           method = "pre.test")
#' }
#'
#' AMR.clean <- correct.meas(post.data = post,
#'                           meas.data = AMR.raw,
#'                           method = "post.test")
#'
#' @references {Svendsen, M. B. S., Bushnell, P. G., & Steffensen, J. F. (2016). Design and setup of intermittent-flow respirometry system for aquatic organisms. Journal of Fish Biology, 88(1), 26-50.}
#'
#' @export
#'
correct.meas <- function (pre.bg, post.bg, meas.data,
                        method = c("pre.test", "post.test", "average",
                                   "linear", "exponential", "parallel", "none"),
                        empty.chamber = c("CH1", "CH2", "CH3", "CH4",
                                          "CH5", "CH6", "CH7", "CH8")){

  method <- match.arg(method)

  M.total <- max(as.numeric(gsub("(F|M)", "", unique(meas.data$Phase))))

  if(method == "pre.test"){

    O2.background <- sapply(1:nrow(meas.data), function(i) {
      row <- which(pre.bg$Chamber.No == meas.data$Chamber.No[i] & pre.bg$Phase.Time == meas.data$Phase.Time[i])
      
      if (length(row) == 0)
        return(NA)
      else
        return(pre.bg$O2.background[row])
    })

    if (any(is.na(O2.background)))
      warning('Some measurement phases are longer than the background. Discarding overextended points', immediate. = TRUE, call. = FALSE)

    meas.data$O2.background <- O2.background
    meas.data <- meas.data[!is.na(O2.background), ]
  }

  if(method == "post.test"){

    O2.background <- sapply(1:nrow(meas.data), function(i) {
      row <- which(post.bg$Chamber.No == meas.data$Chamber.No[i] & post.bg$Phase.Time == meas.data$Phase.Time[i])
      
      if (length(row) == 0)
        return(NA)
      else
        return(post.bg$O2.background[row])
    })

    if (any(is.na(O2.background)))
      warning('Some measurement phases are longer than the background. Discarding overextended points', immediate. = TRUE, call. = FALSE)

    meas.data$O2.background <- O2.background
    meas.data <- meas.data[!is.na(O2.background), ]
  }

  if(method == "average"){
    stop("average has not been updated yet")
    pre.link <- pre.bg$Chamber.No == meas.data$Chamber.No & pre.bg$Phase.Time == meas.data$Phase.time
    post.link <- post.bg$Chamber.No == meas.data$Chamber.No & post.bg$Phase.Time == meas.data$Phase.time

    recipient <- data.frame(pre = pre.bg[pre.link], post = post.bg[post.link])

    recipient$mean <- apply(recipient, 1, function(x) mean(x, na.rm = TRUE))

    meas.data$O2.background <- recipient$mean
  }


  if(method == "linear"){
    stop('this method needs to be updated')

    all.meas <- unique(post.data$Phase)

    # The operation is done by phase and by chamber, so the dataset is broken twice below
    by.chamber <- split(meas.data, meas.data$Chamber.No) # first by chamber

    recipient <- lapply(names(by.chamber), function(the.chamber) { # use the chamber names so the respective pre and post data can be extracted

      pre.lm <- lm(O2.delta.raw ~ Phase.Time, data = pre.data[pre.data$Chamber.No == the.chamber, ])
      post.lm <- lm(O2.delta.raw ~ Phase.Time, data = post.data[post.data$Chamber.No == the.chamber, ])

      by.phase <- split(by.chamber[[the.chamber]], by.chamber[[the.chamber]]$Phase) # now by phase

      recipient <- lapply(names(by.phase), function(the.phase) { # use the phase names so the phase i can be extracted
        phase.i <- as.numeric(gsub("(F|M)", "", the.phase))

        phase.lm <- pre.lm

        # adjust lm weights based on phase distance to the pre and post data
        phase.lm$coefficients[1] <- 0
        phase.lm$coefficients[2] <- (1 - phase.i / (M.total + 1)) * pre.lm$coefficients[2] + phase.i / (M.total + 1) * post.lm$coefficients[2]

        # store values, done.
        by.phase[[the.phase]]$O2.background <- as.vector(predict.lm(phase.lm, by.phase[[the.phase]], type = "response", se.fit = FALSE))

        return(by.phase[[the.phase]])
      })
      # start rebinding back to a dataframe
      output <- data.table::rbindlist(recipient)
      return(output)
    })

    meas.data <- as.data.frame(data.table::rbindlist(recipient))
  }

  if(method == "exponential"){
    stop("exponential has not been updated yet")
    M.phase <- levels(meas.data$Phase)

    b<-NULL
    for (i in M.phase){
      a <- as.numeric(substr(i, 2, 3))
      b <-append(b, a)
    }

    y<-NULL
    for (i in b){
      temp.lm1<-lm(O2.delta.raw ~ Time, data = subset(pre.data, Chamber.No=="CH1"))
      temp.lm2<-lm(O2.delta.raw ~ Time, data = subset(post.data, Chamber.No=="CH1"))
      exp.coef<-sign(temp.lm2$coefficients[2] / temp.lm1$coefficients[2]) * abs(temp.lm2$coefficients[2] / temp.lm1$coefficients[2])^(1 / (M.total + 1))
      pro = temp.lm1$coefficients[2] * exp.coef^i
      temp.lm1$coefficients[1] <- 0
      temp.lm1$coefficients[2] <- pro
      lm.M <- assign(paste("lm.M", i, sep = ""), temp.lm1)
      x<-as.vector(predict.lm(lm.M, temp.df[temp.df$Phase == paste("M", i, sep=""), ], type = "response", se.fit = FALSE))
      y<-append(y, x)
    }
    any(y>0)
    temp.df$O2.background<-y
    rm(list = ls(pattern = "lm.M"))
    rm(list = c("temp.lm1", "temp.lm2", "a", "b", "i", "x", "y", "pro", "exp.coef"))
  }

  if(method == "none")
    meas.data$O2.background <- 0

  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  meas.data$O2.corrected <- meas.data$O2.raw - meas.data$O2.background

  aux <- split(meas.data, paste0(meas.data$Chamber.No, meas.data$Phase))
  aux <- lapply(aux, function(x) {
    x$O2.delta.corrected <- x$O2.corrected - x$O2.corrected[1]
    return(x)
  })

  output <- as.data.frame(data.table::rbindlist(aux))
  
  attributes(output)$correction_method <- method
  
  return(output)
}




calculate.bg <- function(input, method = c('mean', 'first', 'last'), force.linear = TRUE, smoothing = 30){

  method <- match.arg(method)

  n.phases <- unique(input$Phase)

  if (method == 'first')
    input <- input[input$Phase == n.phases[1], ]
    
  if (method == 'last')
    input <- input[input$Phase == n.phases[length(n.phases)], ]

  chamber.lists <- split(input, input$Chamber.No)

  bg.lists <- lapply(names(chamber.lists), function(chamber) {
    # cat(i, '\n')
    sub.data <- input[input$Chamber.No == chamber, ]

    if (force.linear) {
      bg.lm <- lm(O2.delta.raw ~ Phase.Time, data = sub.data)
      bg.lm$coefficients[1] <- 0 

      output <- data.frame(Phase.Time = 1:max(sub.data$Phase.Time))
      output$O2.background <- as.vector(predict.lm(bg.lm, output, type = "response", se.fit = FALSE))
    }
    else {
      output <- aggregate(sub.data$O2.delta.raw, by = list(sub.data$Phase.Time), mean)
      colnames(output) <- c('Phase.Time', 'O2.background')
      
      if (smoothing > 1) {
        x <- filter(output$O2.background, rep(1/smoothing, smoothing), sides = 2)

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






plot_bg <- function(obs, bg, mean.lwd = 1.5) {
  p <- ggplot(data = obs, aes(x = Phase.Time))
  p <- p + geom_line(aes(y = O2.delta.raw, group = Phase, colour = Phase))
  p <- p + geom_line(data = bg, aes(y = O2.background), col = 'red', lwd = mean.lwd)
  p <- p + facet_grid(. ~ Chamber.No)
  p
}


