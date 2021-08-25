#' Dummy documentation
#' 
#' @export
#' 
plot_mr <- function(MR, SMR, MMR, chambers, target_col = 'MR.mass', ylabel = expression(paste("MR (", mg~O[2]~kg^-1~h^-1, ")"))) {
	MR$y_column <- MR[, target_col]

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)

		MR <- MR[MR$Chamber.No %in% chambers, ]
	}

	p <- ggplot2::ggplot(data = MR, ggplot2::aes(x = Date.Time, y = y_column))
	p <- p + ggplot2::geom_line()
	p <- p + ggplot2::geom_point()
	if (!missing(SMR)) {
		SMR$y_column <- SMR[, target_col]
		
		if (!missing(chambers))
			SMR <- SMR[SMR$Chamber.No %in% chambers, ]

		p <- p + ggplot2::geom_hline(data = SMR, ggplot2::aes(yintercept = y_column), col = "red")
	}
	if (!missing(MMR)) {
		MMR$y_column <- MMR[, target_col]		
		
		if (!missing(chambers)) {			
			MMR <- MMR[MMR$Chamber.No %in% chambers, ]
		}

		p <- p + ggplot2::geom_point(data = MMR, ggplot2::aes(x = Date.Time, y = y_column), col = "red", size = 2)
	}
	p <- p + ggplot2::facet_wrap(Chamber.No ~ ., ncol = 1)
	p <- p + ggplot2::labs(y = ylabel, x = '')
	p <- p + ggplot2::theme_bw()
	return(p)	
}


#' Dummy documentation
#' 
#' @export
#' 
plot_slopes <- function(input, chambers, r2 = TRUE) {
	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)
		input <- input[input$Chamber.No %in% chambers, ]
	}

	p <- ggplot2::ggplot(data = input, aes(x = Date.Time))
	p <- p + ggplot2::geom_line(aes(y = Slope, col = 'Slope'))
	p <- p + ggplot2::geom_point(aes(y = Slope, col = 'Slope'))


	if (r2) {
		slope.mean <- mean(input$Slope)
		r2.mean <- 0.5
		mean.dif <- r2.mean - slope.mean
		aux <- range(input$Slope)
		slope.range <- aux[2] - aux[1]
		r2.range <- 1

		if (slope.range == 0 | r2.range == 0)
			compress.factor <- 1
		else
			compress.factor <- slope.range/r2.range

		data.link <- function(x, a = r2.mean, b = slope.mean, c.factor = compress.factor) {
			b + ((x - a) * c.factor) # = y
		}

		axis.link <- function(y, a = r2.mean, b = slope.mean, c.factor = compress.factor) {
		  (y + a * c.factor - b) / c.factor # = x
		}

		p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(R2), col = "R2"))
		p <- p + ggplot2::geom_point(ggplot2::aes(y = data.link(R2), col = "R2"))
 		p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "R2"))
		p <- p + ggplot2::scale_colour_manual(values = c("Grey", "Black"))
		p <- p + labs(colour = 'Values:')
	} else {
		p <- p + ggplot2::theme(legend.position = "none") 
	}

	p <- p + ggplot2::facet_wrap(Chamber.No ~ ., ncol = 1)
	p <- p + ggplot2::theme_bw()
}


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
plot_meas <- function(input, cycles, chambers, temperature = FALSE, oxygen.label = "Oxygen") {
	substrRight <- function(x, n){
	  substr(x, nchar(x)-n+1, nchar(x))
	}

	if (!missing(cycles)) {
		if (!is.numeric(cycles))
			stop("cycles must be a numeric vector")

		n.cycles <- max(as.numeric(gsub("(F|M)", "", unique(input$Phase))))

		if (max(cycles) > n.cycles)
			stop("Requested cycles go over available cycles (", n.cycles, ").", call. = FALSE)
		
		phases <- as.vector(outer(c("M", "F"), cycles, paste0))
		phase.string <- paste0("^", paste(phases, collapse = "$|^"), "$")
		input <- input[grepl(phase.string, input$Phase), ]
		input$Phase <- droplevels(input$Phase)
		rm(phases, phase.string)
	}

	if (!missing(chambers)) {
		if (!is.numeric(chambers))
			chambers <- as.numeric(sub('CH', '', chambers))
		cols.to.keep <- c("Date.Time", "Phase", as.vector(outer(c("Ox.", "Temp."), chambers, paste0)))
		input <- input[,cols.to.keep]
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
		plotdata <- reshape2::melt(to.melt, id.vars = c("Phase", "Date.Time"))
		colnames(plotdata)[grepl("value", colnames(plotdata))] <- "Oxygen"
	}

	plotdata$Chamber <- paste0('CH', substrRight(as.character(plotdata$variable), 1))

	Date.Time <- Temperature <- Oxygen <- Chamber <- NULL
	xmin <- xmax <- ymin <- ymax <- NULL

	p <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Date.Time))

	if (paint_flush)
		p <- p + ggplot2::geom_rect(data = flush.times, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1)

	p <- p + ggplot2::geom_line(ggplot2::aes(y = Oxygen, col = oxygen.label, group = Phase))
	p <- p + ggplot2::theme_bw()

	if (temperature) {
		ox.mean <- mean(plotdata$Oxygen)
		temp.mean <- mean(plotdata$Temperature)
		mean.dif <- temp.mean - ox.mean
		aux <- range(plotdata$Oxygen)
		ox.range <- aux[2] - aux[1]
		aux <- range(plotdata$Temperature)
		temp.range <- aux[2] - aux[1]

		if (ox.range == 0 | temp.range == 0)
			compress.factor <- 1
		else
			compress.factor <- ox.range/temp.range

		data.link <- function(x, a = temp.mean, b = ox.mean, c.factor = compress.factor) {
			b + ((x - a) * c.factor) # = y
		}

		axis.link <- function(y, a = temp.mean, b = ox.mean, c.factor = compress.factor) {
		  (y + a * c.factor - b) / c.factor # = x
		}

		p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(Temperature), col = "Temperature", group = Phase))
 		p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "Temperature"))
		p <- p + ggplot2::scale_colour_manual(values = c("royalblue", "red"))
	} else {
		p <- p + ggplot2::scale_colour_manual(values = "blue") 
	}
	p <- p + ggplot2::theme(legend.position = "none") 

	p <- p + ggplot2::labs(y = oxygen.label, x = "Time")
	p <- p + ggplot2::facet_wrap(Chamber ~ ., ncol = 1)

	return(p)
}



#' Dummy documentation
#' 
#' @export
#' 
plot_deltas <- function(input, cycles, chambers) {
	if (!missing(cycles)) {
		if (!is.numeric(cycles))
			stop("cycles must be a numeric vector")

		n.cycles <- max(as.numeric(gsub("(F|M)", "", unique(input$Phase))))

		if (max(cycles) > n.cycles)
			stop("Requested cycles go over available cycles (", n.cycles, ").", call. = FALSE)
		
		phases <- as.vector(outer(c("M", "F"), cycles, paste0))
		phase.string <- paste0("^", paste(phases, collapse = "$|^"), "$")
		input <- input[grepl(phase.string, input$Phase), ]
		rm(phases, phase.string)
	}

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)
		input <- input[input$Chamber.No %in% chambers, ]
	}


	p <- ggplot(data = input, aes(x = Date.Time, Group = Phase))
	p <- p + geom_path(aes(y = O2.background, col = 'Background'))
	p <- p + geom_path(data = input, aes(y = O2.delta.raw, col = 'Raw'))
	p <- p + geom_path(data = input, aes(y = O2.delta.corrected, col = 'Corrected'))
	p <- p + theme_bw()
	p <- p + ggplot2::scale_colour_manual(values = c("Grey", "royalblue", 'black'))
	# p <- p + ggplot2::scale_linetype_manual(values = c("dashed", "solid", 'solid'))

	p <- p + facet_wrap(.~Chamber.No)
	p <- p + labs(x = '', y = expression(Delta~O[2]), title = paste('Correction method used:', attributes(input)$correction_method),
				  colour = 'Values:', linetype = 'Values:')
	p
}