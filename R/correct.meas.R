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
correct.meas <- function (pre.bg, post.bg, meas.data, O2_col = 'O2.raw',
                        method = c("pre.test", "post.test", "average",
                                   "linear", "exponential", "parallel", "none"),
                        empty.chamber = c("CH1", "CH2", "CH3", "CH4",
                                          "CH5", "CH6", "CH7", "CH8")){

  method <- match.arg(method)

  M.total <- max(meas.data$Cycle)

  if(method == "pre.test"){
    my_bg <- pre.bg

    meas.data$temporary_index <- paste(meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Probe, my_bg$Phase.Time)
  }

  if(method == "post.test"){
    my_bg <- post.bg

    meas.data$temporary_index <- paste(meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Probe, my_bg$Phase.Time)
  }

  if(method == "average"){
    stop("average has not been updated yet")
    pre.link <- pre.bg$Probe == meas.data$Probe & pre.bg$Phase.Time == meas.data$Phase.time
    post.link <- post.bg$Probe == meas.data$Probe & post.bg$Phase.Time == meas.data$Phase.time

    recipient <- data.frame(pre = pre.bg[pre.link], post = post.bg[post.link])

    recipient$mean <- apply(recipient, 1, function(x) mean(x, na.rm = TRUE))

    meas.data$O2.background <- recipient$mean
  }


  if(method == "linear"){
    my_bg <- calculate_linear_bg_progression(pre.bg = pre.bg, post.bg = post.bg, meas.data = meas.data)

    meas.data$temporary_index <- paste(meas.data$Cycle, meas.data$Probe, meas.data$Phase.Time)

    my_bg$temporary_index <- paste(my_bg$Cycle, my_bg$Probe, my_bg$Phase.Time)

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
      temp.lm1<-lm(O2.background ~ Phase.Time, data = subset(pre.data, Probe=="CH1"))
      temp.lm2<-lm(O2.background ~ Phase.Time, data = subset(post.data, Probe=="CH1"))
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

  if(method == "none") {
    meas.data$O2.background <- 0
  } else {
    link <- match(meas.data$temporary_index, my_bg$temporary_index)

    meas.data$O2.background <- my_bg$O2.background[link]

    meas.data$temporary_index <- NULL
    
    if (any(is.na(meas.data$O2.background)))
      warning('Some measurement phases are longer than the background. Discarding overextended points', immediate. = TRUE, call. = FALSE)

    meas.data <- meas.data[!is.na(meas.data$O2.background), ]
  }

  #--------------------------------------------------------------------------------------------------------------------------------------------------#
  meas.data$O2.corrected <- meas.data[, O2_col] - meas.data$O2.background

  aux <- split(meas.data, paste0(meas.data$Probe, meas.data$Phase))
  aux <- lapply(aux, function(x) {
    x$O2.delta.corrected <- x$O2.corrected - x$O2.corrected[1]
    return(x)
  })

  output <- as.data.frame(data.table::rbindlist(aux))
  
  attributes(output)$correction_method <- method
  
  return(output)
}







calculate_linear_bg_progression <- function(pre.bg, post.bg, meas.data) {
  cycles <- max(meas.data$Cycle)

  my_bg <- lapply(unique(meas.data$Probe), function(chamber) {
    # cat(chamber, "\n")
    pre_line <- pre.bg$O2.background[pre.bg$Probe == chamber]

    post_line <- post.bg$O2.background[post.bg$Probe == chamber]

    if (length(pre_line) != length(post_line)) {
      warning("The pre-bg and post-bg in chamber ", chamber, " have different lengths! Truncating longer vector.", immediate. = TRUE, call. = FALSE)
      
      if (length(pre_line) < length(post_line))
        post_line <- post_line[1:length(pre_line)]
      else
        pre_line <- pre_line[1:length(post_line)]
    }

    difference <- post_line - pre_line

    if (all(difference == 0)) {
      stop("The pre-bg and post-bg are exactly the same! Cannot calculate linear progression.", call. = FALSE)
    }

    increments <- difference / (cycles - 1)

    my_list <- lapply(1:length(pre_line), function(i) {
      seq(from = pre_line[i], to = post_line[i], by = increments[i])
    })

    my_matrix <- as.data.frame(do.call(rbind, my_list))

    return(my_matrix)
  })
  names(my_bg) <- unique(meas.data$Probe)



  my_bg_simplified <- lapply(names(my_bg), function(chamber) {
    x <- suppressMessages(reshape2::melt(my_bg[[chamber]]))
    colnames(x) <- c("Cycle", "O2.background")
    x$Phase.Time <- 1:nrow(my_bg[[chamber]])
    x$Cycle <- as.numeric(sub("V", "", x$Cycle))
    x$Probe <- chamber
    return(x)
  })

  my_bg_df <- do.call(rbind, my_bg_simplified)
  return(my_bg_df)
}