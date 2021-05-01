#' Plot raw measurements of oxygen and temperature
#' 
#' @param input A dataframe of measurements imported using one of FishResp's loading functions
#' @param phases A vector of phases to plot (leave empty to plot all phases)
#' @param temperature Logical: Should the temperature recordings be plotted?
#' @param oxygen.label The label of the Oxygen axis. Defaults to "Oxygen"
#' 
#' @return A ggplot object
#' 
#' @export
#' 
plot_meas <- function(input, phases, temperature = TRUE, oxygen.label = "Oxygen") {
	substrRight <- function(x, n){
	  substr(x, nchar(x)-n+1, nchar(x))
	}

	if (!missing(phases)) {
		if (!is.numeric(phases))
			stop("phases must be a numeric vector")
		phases <- as.vector(outer(c("M", "F"), phases, paste0))
		phase.string <- paste0("^", paste(phases, collapse = "$|^"), "$")
		input <- input[grepl(phase.string, input$Phase), ]
		input$Phase <- droplevels(input$Phase)
		rm(phases, phase.string)
	}

	if (any(grepl("F", input$Phase))) {
		paint_flush <- TRUE
		aux <- split(input, input$Phase)
		aux <- aux[grepl("F", names(aux))]
		aux <- lapply(aux, function(x) {
			data.frame(xmin = x$Date.Time[1], xmax = x$Date.Time[nrow(x)], ymin = -Inf, ymax = Inf)
		})
		flush.times <- as.data.frame(data.table::rbindlist(aux))
		rm(aux)
	} else {
		paint_flush <- FALSE
	}

	if (temperature) {
		to.melt.ox   <- input[, grepl("^Date.Time$|^Phase$|^Ox.[0-9]", colnames(input))]
		to.melt.temp <- input[, grepl("^Date.Time$|^Phase$|^Temp.[0-9]$", colnames(input))]
		
		aux.ox <- reshape2::melt(to.melt.ox, id.vars = c("Phase", "Date.Time"))
		colnames(aux.ox)[grepl("value", colnames(aux.ox))] <- "Oxygen"

		aux.temp <- reshape2::melt(to.melt.temp, id.vars = c("Phase", "Date.Time"))
		colnames(aux.temp)[grepl("value", colnames(aux.temp))] <- "Temperature"
		
		plotdata <- cbind(aux.ox, Temperature = aux.temp$Temperature)

	}	else {
		to.melt <- input[, grepl("^Date.Time$|^Phase$|^Ox.[0-9]", colnames(input))]
		colnames(aux.ox)[grepl("value", colnames(aux.ox))] <- "Oxygen"
	}

	plotdata$Chamber <- as.numeric(substrRight(as.character(plotdata$variable), 1))

	Date.Time <- Temperature <- Oxygen <- Chamber <- NULL
	xmin <- xmax <- ymin <- ymax <- NULL

	p <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Date.Time))

	if (paint_flush)
		p <- p + ggplot2::geom_rect(data = flush.times, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1)

	p <- p + ggplot2::geom_line(ggplot2::aes(y = Oxygen, col = oxygen.label))

	if (temperature) {
		ox.mean <- mean(plotdata$Oxygen)
		temp.mean <- mean(plotdata$Temperature)
		mean.dif <- temp.mean - ox.mean
		aux <- range(plotdata$Oxygen)
		ox.range <- aux[2] - aux[1]
		aux <- range(plotdata$Temperature)
		temp.range <- aux[2] - aux[1]
		compress.factor <- ox.range/temp.range

		data.link <- function(x, a = temp.mean, b = ox.mean, c.factor = compress.factor) {
			b + ((x - a) * c.factor) # = y
		}

		axis.link <- function(y, a = temp.mean, b = ox.mean, c.factor = compress.factor) {
		  (y + a * c.factor - b) / c.factor # = x
		}

		p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(Temperature), col = "Temperature"))
 		p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "Temperature"))
		p <- p + ggplot2::scale_colour_manual(values = c("blue", "red"))
	} else {
		p <- p + ggplot2::scale_colour_manual(values = "blue") 
	}

	p <- p + ggplot2::labs(y = oxygen.label, x = "Time")

	p <- p + ggplot2::facet_wrap(Chamber ~ ., ncol = 1)

	p <- p + ggplot2::theme_bw()

	p

	return(p)
}

