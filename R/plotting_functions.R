#' Dummy documentation
#' 
#' @export
#' 
plot_mr <- function(MR, SMR = NULL, MMR = NULL, cycles, chambers, target_col = 'MR.mass', ylabel = expression(paste("MR (", mg~O[2]~kg^-1~h^-1, ")"))) {
	Date.Time <- NULL
	y_column <- NULL

	MR$y_column <- MR[, target_col]

	if (!missing(cycles)) {
		if (!is.numeric(cycles))
			stop("cycles must be a numeric vector")

		n.cycles <- max(as.numeric(gsub("(F|M)", "", unique(MR$Phase))))

		if (max(cycles) > n.cycles)
			stop("Requested cycles go over available cycles (", n.cycles, ").", call. = FALSE)
		
		phases <- as.vector(outer(c("M", "F"), cycles, paste0))
		phase.string <- paste0("^", paste(phases, collapse = "$|^"), "$")
		MR <- MR[grepl(phase.string, MR$Phase), ]
		# MR$Phase <- droplevels(MR$Phase)
		rm(phases, phase.string)
	}

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)

		MR <- MR[MR$Probe %in% chambers, ]
	}

	p <- ggplot2::ggplot()
	p <- p + ggplot2::geom_line(data = MR, ggplot2::aes(x = Date.Time, y = y_column))
	p <- p + ggplot2::geom_point(data = MR, ggplot2::aes(x = Date.Time, y = y_column))

	if (!is.null(SMR)) {
		SMR$y_column <- SMR[, target_col]
		
		if (!missing(chambers)) 
			SMR <- SMR[SMR$Probe %in% chambers, ]

		p <- p + ggplot2::geom_hline(data = SMR, ggplot2::aes(yintercept = y_column), col = "red")
	}
	
	if (!is.null(MMR)) {
		MMR$y_column <- MMR[, target_col]		
		
		if (!missing(chambers)) {			
			MMR <- MMR[MMR$Probe %in% chambers, ]
		}

		p <- p + ggplot2::geom_point(data = MMR, ggplot2::aes(x = Date.Time, y = y_column), col = "red", size = 2)
	}
	
	p <- p + ggplot2::facet_wrap(Probe ~ ., ncol = 1)
	p <- p + ggplot2::labs(y = ylabel, x = '')
	p <- p + ggplot2::theme_bw()
	return(p)	
}


#' Dummy documentation
#' 
#' @export
#' 
plot_slopes <- function(input, cycles, chambers, r2 = TRUE) {
	Date.Time <- NULL
	Slope <- NULL
	R2 <- NULL

	if (!missing(cycles)) {
		if (!is.numeric(cycles))
			stop("cycles must be a numeric vector")

		n.cycles <- max(as.numeric(gsub("(F|M)", "", unique(input$Phase))))

		if (max(cycles) > n.cycles)
			stop("Requested cycles go over available cycles (", n.cycles, ").", call. = FALSE)
		
		phases <- as.vector(outer(c("M", "F"), cycles, paste0))
		phase.string <- paste0("^", paste(phases, collapse = "$|^"), "$")
		input <- input[grepl(phase.string, input$Phase), ]
		# input$Phase <- droplevels(input$Phase)
		rm(phases, phase.string)
	}

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)
		input <- input[input$Probe %in% chambers, ]
	}

	p <- ggplot2::ggplot(data = input, ggplot2::aes(x = Date.Time))
	p <- p + ggplot2::geom_line(ggplot2::aes(y = Slope, col = 'Slope'))
	p <- p + ggplot2::geom_point(ggplot2::aes(y = Slope, col = 'Slope'))


	if (r2) {
		slope.mean <- mean(input$Slope, na.rm = TRUE)
		r2.mean <- 0.5
		mean.dif <- r2.mean - slope.mean
		aux <- range(input$Slope, na.rm = TRUE)
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
		p <- p + ggplot2::labs(colour = 'Values:')
	} else {
		p <- p + ggplot2::theme(legend.position = "none") 
	}

	p <- p + ggplot2::facet_wrap(Probe ~ ., ncol = 1)
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
plot_meas <- function(input, cycles, chambers, temperature = FALSE, oxygen.col = "O2.hPa", oxygen.label = "Oxygen (hPa)") {
	Phase <- NULL
	Temp <- NULL

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
		# input$Phase <- droplevels(input$Phase)
		rm(phases, phase.string)
	}

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)

		input <- input[input$Probe %in% chambers, ]
	}

	if (any(grepl("F", input$Phase))) {
		paint_flush <- TRUE

		aux <- split(input, input$Probe)
		
		recipient <- lapply(names(aux), function(the.chamber) {
			# instead of using split I need to manually break these
			# because F0 can appear multiple times. split would
			# group them together and mess up the plot
			breaks <- c(1, cumsum(rle(aux[[the.chamber]]$Phase)$lengths))
			aux2 <- lapply(2:length(breaks), function(i) {
				input[breaks[i-1]:breaks[i], ]
			})
			# --
			names(aux2) <- rle(aux[[the.chamber]]$Phase)$values
			aux2 <- aux2[grepl("F", names(aux2))]
			# now extract start and end of phase
			aux2 <- lapply(aux2, function(the.phase) {
				data.frame(Probe = the.chamber, xmin = the.phase$Date.Time[1], xmax = the.phase$Date.Time[nrow(the.phase)], ymin = -Inf, ymax = Inf)
			})
			return(as.data.frame(data.table::rbindlist(aux2)))
		})
		rm(aux)

		flush.times <- as.data.frame(data.table::rbindlist(recipient))
	} else {
		paint_flush <- FALSE
	}

	input$Oxygen <- input[, oxygen.col] #to eliminate the units part

	Date.Time <- Temperature <- Oxygen <- Chamber <- NULL
	xmin <- xmax <- ymin <- ymax <- NULL

	p <- ggplot2::ggplot(data = input, ggplot2::aes(x = Date.Time))

	if (paint_flush)
		p <- p + ggplot2::geom_rect(data = flush.times, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.1)

	p <- p + ggplot2::geom_line(ggplot2::aes(y = Oxygen, col = oxygen.label, group = Phase))
	p <- p + ggplot2::theme_bw()

	if (temperature) {
		aux <- range(input$Oxygen, na.rm = TRUE)
		ox.range <- aux[2] - aux[1]
		ox.mid <- mean(aux)

		aux <- range(input$Temp, na.rm = TRUE)
		temp.range <- aux[2] - aux[1]
		temp.mid <- mean(aux)

		mid.dif <- temp.mid - ox.mid

		if (ox.range == 0 | temp.range == 0) {
			compress.factor <- 1
		} else {
			compress.factor <- ox.range/temp.range
		}

		data.link <- function(x, a = temp.mid, b = ox.mid, c.factor = compress.factor) {
			b + ((x - a) * c.factor) # = y
		}

		axis.link <- function(y, a = temp.mid, b = ox.mid, c.factor = compress.factor) {
		  (y + a * c.factor - b) / c.factor # = x
		}

		p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(Temp), col = "Temperature", group = Phase))
 		p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "Temperature"))
		p <- p + ggplot2::scale_colour_manual(values = c("royalblue", "red"))
	} else {
		p <- p + ggplot2::scale_colour_manual(values = "blue") 
	}
	p <- p + ggplot2::theme(legend.position = "none") 

	p <- p + ggplot2::labs(y = oxygen.label, x = "Time")
	p <- p + ggplot2::facet_wrap(Probe ~ ., ncol = 1)

	return(p)
}



#' Dummy documentation
#' 
#' @export
#' 
plot_deltas <- function(input, cycles, chambers, raw_delta_col = 'O2.delta.raw') {
	Date.Time <- NULL
	Phase <- NULL
	O2.background <- NULL
	plot_this <- NULL
	O2.delta.corrected <- NULL

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
		input <- input[input$Probe %in% chambers, ]
	}

	input$plot_this <- input[, raw_delta_col]
	p <- ggplot2::ggplot(data = input, ggplot2::aes(x = Date.Time, Group = Phase))
	p <- p + ggplot2::geom_path(ggplot2::aes(y = O2.background, col = 'Background'))
	p <- p + ggplot2::geom_path(data = input, ggplot2::aes(y = plot_this, col = 'Raw'))
	p <- p + ggplot2::geom_path(data = input, ggplot2::aes(y = O2.delta.corrected, col = 'Corrected'))
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::scale_colour_manual(values = c("Grey", "royalblue", 'black'))
	# p <- p + ggplot2::scale_linetype_manual(values = c("dashed", "solid", 'solid'))

	p <- p + ggplot2::facet_wrap(.~Probe)
	p <- p + ggplot2::labs(x = '', y = expression(Delta~O[2]~"(umol/L)"), title = paste('Correction method used:', attributes(input)$correction_method),
				  colour = 'Values:', linetype = 'Values:')
	p
}


#' dummy documentation
#' 
#' plots the whole experiment data, from background to MO2
#' 
#' @export
#' 
plot_experiment <- function(pre, post, mr, cycles, chamber, smr = FALSE, mmr = FALSE, raw_delta_col = 'O2.delta.raw',
							title = 'Experiment measurements',
							mr.col = 'MR.mass', mr.label = 'MO2', verbose = FALSE) {

	if (verbose)
		message('Plotting Pre-background')

	if (!missing(pre))
		p_pre <- plot_bg(pre$meas, O2_col = 'O2.delta.umol.l', pre$bg, chambers = chamber) + ggplot2::labs(title = 'Pre-background')
	else
		p_pre <- patchwork::wrap_elements(grid::textGrob('Pre-background was not recorded'))
	
	if (verbose)
		message('Plotting Post-background')

	if (!missing(post))
		p_post <- plot_bg(post$meas, O2_col = 'O2.delta.umol.l', post$bg, chambers = chamber) + ggplot2::labs(title = 'Post-background')
	else
		p_post <- patchwork::wrap_elements(grid::textGrob('Post-background was not recorded'))

	if (verbose)
		message('Plotting measurements')

	p1 <- plot_meas(mr$meas_raw, oxygen.col = "O2.umol.l", oxygen.label = "O2 (umol/L)", cycles = cycles, chambers = chamber, temperature = TRUE) 
	p1 <- p1 + ggplot2::labs(title = title, x = '')
	
	if (verbose)
		message('Plotting deltas')

	B <- plot_deltas(mr$corrected, cycles = cycles, chambers = chamber, raw_delta_col = raw_delta_col) + mimic_x_datetime(p1)

	if (verbose)
		message('Plotting plotting slopes')

	p2 <- plot_slopes(mr$all.slopes, cycles = cycles, chambers = chamber) + mimic_x_datetime(p1) + ggplot2::xlab('')

	if (smr)
		the_smr <- mr$smr
	else
		the_smr <- NULL

	if (mmr)
		the_mmr <- mr$mmr
	else
		the_mmr <- NULL

	if (verbose)
		message('Plotting metabolic rate')

	p3 <- plot_mr(MR = mr$mr, SMR = the_smr, MMR = the_mmr, chambers = chamber, 
		  target_col = mr.col, ylabel = mr.label)
	p3 <- p3 + mimic_x_datetime(p1)

	to.print <- p_pre + p_post + p1 + B + p2 + p3 + patchwork::plot_layout(design = 'AB\nCC\nDD\nEE\nFF')

	return(to.print)
}

mimic_x_datetime <- function(input) {
	ggplot2::xlim(as.POSIXct(ggplot2::layer_scales(input)$x$range$range, origin = "1970-01-01 00:00.00"))	
}

mimic_y_datetime <- function(input) {
	ggplot2::ylim(as.POSIXct(ggplot2::layer_scales(input)$y$range$range, origin = "1970-01-01 00:00.00"))	
}

mimic_x <- function(input) {
	ggplot2::xlim(ggplot2::layer_scales(input)$x$range$range)
}

mimic_y <- function(input) {
	ggplot2::ylim(ggplot2::layer_scales(input)$y$range$range)
}






#' dummy doc
#' 
#' @export
#' 
plot_bg <- function(obs, bg, chambers, O2_col, mean.lwd = 1.5) {
	Phase.Time <- NULL
	O2.background <- NULL

	if (!missing(chambers)) {
		if (is.numeric(chambers))
			chambers <- paste0('CH', chambers)

		obs <- obs[obs$Probe %in% chambers, ]
		bg <- bg[bg$Probe %in% chambers, ]
	}

	obs$Phase <- sub('M', '', obs$Phase)
	
	p <- ggplot2::ggplot(data = obs, ggplot2::aes(x = Phase.Time))
	p <- p + ggplot2::geom_line(ggplot2::aes_string(y = O2_col, group = "Phase", colour = "Phase"))
	p <- p + ggplot2::geom_line(data = bg, ggplot2::aes(y = O2.background), col = 'red', lwd = mean.lwd)
	p <- p + ggplot2::facet_wrap(. ~ Probe, ncol = 4)
	p <- p + ggplot2::theme_bw()
	# p <- p + theme(legend.position = 'none')
	p
}



#' Plot rolling mmr
#' 
#' @param rolling_mmr the output of compile_rolling_mmr (which can have been edited by  exclude_rolling_mmr_segment)
#' 
#' @export
#' 
plot_rolling_mmr <- function(rolling_mmr) {
  Probe <- NULL
  R2 <- NULL
  Sec <- NULL
  MR.mass.umol.g <- NULL

  p <- ggplot2::ggplot(data = rolling_mmr$detailed_mmr, ggplot2::aes(x = Sec, y = MR.mass.umol.g, group = Probe))
  p <- p + ggplot2::geom_line()
  p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), col = "grey")
  p <- p + ggplot2::geom_point(data = rolling_mmr$max_mmr, col = 'red')
  p <- p + ggplot2::ylab(expression(paste("MR (", mu, mol~O[2]~g^-1~h^-1, ")"))) + ggplot2::xlab('Time (seconds)')

  mmr.mean <- mean(rolling_mmr$detailed_mmr$MR.mass.umol.g, na.rm = TRUE)
  r2.mean <- 0.5
  mean.dif <- r2.mean - mmr.mean
  aux <- range(rolling_mmr$detailed_mmr$MR.mass.umol.g, na.rm = TRUE)
  mmr.range <- aux[2] - aux[1]
  r2.range <- 1

  if (mmr.range == 0 | r2.range == 0) {
    compress.factor <- 1
  } else {
    compress.factor <- mmr.range/r2.range
  }

  data.link <- function(x, a = r2.mean, b = mmr.mean, c.factor = compress.factor) {
    b + ((x - a) * c.factor) # = y
  }

  axis.link <- function(y, a = r2.mean, b = mmr.mean, c.factor = compress.factor) {
    (y + a * c.factor - b) / c.factor # = x
  }

  p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(R2)), col = "grey")
  p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "R2"))

  p <- p + ggplot2::facet_wrap(Probe ~ .)
  p <- p + ggplot2::theme_bw()
  return(p)
}


plot_o2 <- function(folder) {
	Timestamp <- NULL
	Oxygen <- NULL
	Temperature <- NULL

	files <- list.files(paste0(folder, '/ChannelData/'))

	O2.file.link <- grepl('Oxygen.txt', files)

	if (all(!O2.file.link)) {
		stop('No oxygen files found')
	}

	O2.files <- files[O2.file.link]

	aux <- lapply(O2.files, function(i) {
		x <- load_pyro_o2_file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
		x[, c('Date.Time', 'Oxygen.Main', 'Sample.CompT')]
	})
	names(aux) <- stringr::str_extract(O2.files,'Ch.[0-9]')

	plotdata <- data.table::rbindlist(aux, idcol = 'Chamber')
	colnames(plotdata)[2:4] <- c('Timestamp', 'Oxygen', 'Temperature')

	p <- ggplot2::ggplot(data = plotdata, ggplot2::aes(x = Timestamp, y = Oxygen, col = 'Oxygen'))
	p <- p + ggplot2::geom_line()

	aux <- range(plotdata$Oxygen)
	ox.range <- aux[2] - aux[1]
	ox.mid <- mean(aux)

	aux <- range(plotdata$Temperature)
	temp.range <- aux[2] - aux[1]
	temp.mid <- mean(aux)

	mid.dif <- temp.mid - ox.mid

	if (ox.range == 0 | temp.range == 0) {
		compress.factor <- 1
	} else {
		compress.factor <- ox.range/temp.range
	}

	data.link <- function(x, a = temp.mid, b = ox.mid, c.factor = compress.factor) {
		b + ((x - a) * c.factor) # = y
	}

	axis.link <- function(y, a = temp.mid, b = ox.mid, c.factor = compress.factor) {
	  (y + a * c.factor - b) / c.factor # = x
	}

	p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(Temperature), col = "Temperature"))
	p <- p + ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = axis.link, name = "Temperature"))
	p <- p + ggplot2::scale_colour_manual(values = c("royalblue", "red"))

	p <- p + ggplot2::facet_wrap(Chamber ~ ., ncol = 1)
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::theme(legend.position = "none") 

	return(p)
}