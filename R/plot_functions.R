#' Plot background readings
#' 
#' @param input An experiment list with calculated background.
#' 	The output of \code{\link{calc_bg}}.
#' @param linewidth The width of the averaged/linear bg.
#' @inheritParams plot_meas
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_bg <- function(input, probes, linewidth = 1.5) {
	# ggplot variables
	o2_delta <- cycle <- o2_bg_delta <- phase_time <- NULL

	if(is.null(input$bg))
		stop("Could not find a bg object in the input")

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, input$cleaned$probe, "probes")
		input$cleaned <- input$cleaned[input$cleaned$probe %in% probes, ]
		input$bg <- input$bg[input$bg$probe %in% probes, ]
	}

	input$cleaned$cycle <- as.factor(input$cleaned$cycle)
	input$probe <- factor(input$idprobe, levels = input$probe_info$probe)
	
	p <- ggplot2::ggplot(data = input$cleaned, ggplot2::aes(x = phase_time))

	p <- p + ggplot2::geom_line(ggplot2::aes(y = o2_delta, 
								group = cycle,
								colour = cycle))
	p <- p + ggplot2::geom_line(data = input$bg, 
								ggplot2::aes(y = o2_bg_delta),
								col = 'red',
								linewidth = linewidth)
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::labs(x = "Phase time", y = "Delta~O[2]")
	p <- p + ggplot2::facet_wrap(. ~ probe, ncol = 4)
	
	p
}

#' Plot raw measurements of oxygen and temp
#' 
#' @param input An experiment list with melted oxygen and temperature data.
#' 	The output of \code{\link{process_experiment}} or
#' 	\code{\link{melt_resp}}.
#' @param cycles A numeric vector of which cycles to plot
#' @param probes A string of which probes to plot
#' @param show_temp Should the temperature recordings be plotted?
#' @param verbose should argument warnings be displayed?
#'
#' @return A ggplot object
#' 
#' @export
#' 
plot_meas <- function(input, cycles, probes,
					  show_temp = FALSE, verbose = TRUE) {
	# ggplot variables
	phase <- date_time <- temp <- o2 <- probe <- NULL
	xmin <- xmax <- ymin <- ymax <- NULL

	meas <- input$melted
	meas$idprobe <- paste0(meas$id, " (", meas$probe, ")")
	idprobe_levels <-paste0(input$probe_info$id, 
							" (", input$probe_info$probe, ")")
	meas$idprobe <- factor(meas$idprobe, levels = idprobe_levels)
	
	# strip units to avoid conflicts
	o2_unit <- units(meas$o2)
	meas$o2 <- as.numeric(meas$o2)
	temp_unit <- units(meas$temp)
	meas$temp <- as.numeric(meas$temp)

	if (!missing(cycles)) {
		cycles <- check_arg_in_data(cycles, meas$cycle,
									"cycles", verbose = verbose)
		meas <- meas[meas$cycle %in% cycles, ]
	}

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, meas$probe,
									"probes", verbose = verbose)
		meas <- meas[meas$probe %in% probes, ]
		meas$idprobe <- droplevels(meas$idprobe)
	}

	if (any(grepl("F", meas$phase))) {
		paint_flush <- TRUE

		aux <- split(meas, meas$idprobe)
		
		recipient <- lapply(names(aux), function(the_idprobe) {
			# instead of using split, I need to manually break these
			# because F0 can appear multiple times. split would
			# group them together and mess up the plot
			breaks <- c(1, cumsum(rle(aux[[the_idprobe]]$phase)$lengths))
			aux2 <- lapply(2:length(breaks), function(i) {
				meas[breaks[i-1]:breaks[i], ]
			})
			# --
			names(aux2) <- rle(aux[[the_idprobe]]$phase)$values
			aux2 <- aux2[grepl("F", names(aux2))]
			# now extract start and end of phase
			aux2 <- lapply(aux2, function(the_phase) {
				data.frame(idprobe = the_idprobe,
						   xmin = the_phase$date_time[1], 
						   xmax = the_phase$date_time[nrow(the_phase)], 
						   ymin = -Inf, 
						   ymax = Inf)
			})
			return(as.data.frame(data.table::rbindlist(aux2)))
		})
		rm(aux)

		flush_times <- as.data.frame(data.table::rbindlist(recipient))
	} else {
		paint_flush <- FALSE
	}

	p <- ggplot2::ggplot(data = meas)

	if (paint_flush) {
		p <- p + ggplot2::geom_rect(data = flush_times, 
									ggplot2::aes(xmin = xmin, xmax = xmax, 
												 ymin = ymin, ymax = ymax),
									fill = "black", alpha = 0.1)
	}

	p <- p + ggplot2::geom_line(ggplot2::aes(x = date_time, 
											 y = o2, group = phase,
											 colour = "o2"))
	p <- p + ggplot2::theme_bw()

	if (show_temp) {
		aux <- range(meas$o2, na.rm = TRUE)
		ox_range <- aux[2] - aux[1]
		ox_mid <- mean(aux)

		aux <-  range(meas$temp, na.rm = TRUE)
		temp_range <- aux[2] - aux[1]
		temp_mid <- mean(aux)

		mid_dif <- temp_mid - ox_mid

		if (ox_range == 0 | temp_range == 0) {
			compress_factor <- 1
		} else {
			compress_factor <- ox_range/temp_range
		}

		data.link <- function(x, a = temp_mid, 
							  b = ox_mid, c_factor = compress_factor) {
			b + ((x - a) * c_factor) # = y
		}

		axis_link <- function(y, a = temp_mid, 
							  b = ox_mid, c_factor = compress_factor) {
		  (y + a * c_factor - b) / c_factor # = x
		}

		p <- p + ggplot2::geom_line(ggplot2::aes(x = date_time,
												 y = data.link(temp), 
												 group = phase,
												 colour = "temp"))
 		p <- p + ggplot2::scale_y_continuous(
			sec.axis = ggplot2::sec_axis(trans = axis_link, 
		 	name = units::make_unit_label("Temperature", temp_unit))
								)
		p <- p + ggplot2::scale_colour_manual(values = c("royalblue", "red"))
	} else {
		p <- p + ggplot2::scale_colour_manual(values = "royalblue") 
	}

	p <- p + ggplot2::theme(legend.position = "none") 
	p <- p + ggplot2::labs(y = units::make_unit_label("O[2]", o2_unit),
						   x = "Time")
	p <- p + ggplot2::facet_wrap(idprobe ~ ., ncol = 1)

	return(p)
}

#' Plot the calculated oxygen deltas
#' 
#' @param input An experiment list with calculated deltas.
#' 	The output of \code{\link{process_experiment}} or
#' 	\code{\link{calc_delta}}.
#' @inheritParams plot_meas
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_deltas <- function(input, cycles, probes, verbose = TRUE) {
	# ggplot variables
	o2_delta <- o2_cordelta <- o2_bg_delta <- date_time <- NULL
	phase <- plot_this <- NULL

	cleaned <- input$cleaned
	cleaned$idprobe <- paste0(cleaned$id, " (", cleaned$probe, ")")
	idprobe_levels <-paste0(input$probe_info$id, 
							" (", input$probe_info$probe, ")")
	cleaned$idprobe <- factor(cleaned$idprobe, levels = idprobe_levels)

	if (!missing(cycles)) {
		cycles <- check_arg_in_data(cycles, cleaned$cycle,
									"cycles", verbose = verbose)
		cleaned <- cleaned[cleaned$cycle %in% cycles, ]
	}

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, cleaned$probe,
									"probes", verbose = verbose)
		cleaned <- cleaned[cleaned$probe %in% probes, ]
		cleaned$idprobe <- droplevels(cleaned$idprobe)
	}

	p <- ggplot2::ggplot(data = cleaned, 
		ggplot2::aes(x = date_time, Group = phase))
	p <- p + ggplot2::geom_path(ggplot2::aes(y = o2_bg_delta, col = 'Background'))
	p <- p + ggplot2::geom_path(ggplot2::aes(y = o2_delta, 
										     col = 'Raw'))
	p <- p + ggplot2::geom_path(ggplot2::aes(y = o2_cordelta,
											 col = 'Corrected'))
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::scale_colour_manual(values = c("Grey", "royalblue", 
													 "black"))
	p <- p + ggplot2::facet_wrap(.~idprobe)
	p <- p + ggplot2::labs(x = "", y = expression(Delta~O[2]), 
						   title = paste('Correction method used:',
						 				 attributes(input$cleaned)$correction_method),
						   colour = 'Values:', linetype = 'Values:')
	
	return(p)
}

#' Plot the calculated slopes
#' 
#' @param input An experiment list with calculated slopes.
#' 	The output of \code{\link{process_experiment}} or
#' 	\code{\link{calc_slopes}}.
#' @param r2 show the respective R2 for each slope?
#' @inheritParams plot_mr
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_slopes <- function(input, cycles, probes, r2 = TRUE, verbose = TRUE) {
	# ggplot variables
	date_time <- slope_cor <- valid <- NULL

	slopes <- input$all_slopes
	slopes$idprobe <- paste0(slopes$id, " (", slopes$probe, ")")
	idprobe_levels <-paste0(input$probe_info$id, 
							" (", input$probe_info$probe, ")")
	slopes$idprobe <- factor(slopes$idprobe, levels = idprobe_levels)
	
	# remove units to make the double y-axis easy
	slope_unit <- units(slopes$slope_cor)
	slopes$slope_cor <- as.numeric(slopes$slope_cor)
	r2_threshold <- attributes(input$good_slopes)$r2_threshold
	slopes$valid <- slopes$r2 > r2_threshold
	slopes$valid[is.na(slopes$valid)] <- FALSE

	if (!missing(cycles)) {
		cycles <- check_arg_in_data(cycles, slopes$cycle,
									"cycles", verbose = verbose)
		slopes <- slopes[slopes$cycle %in% cycles, ]
	}

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, slopes$probe,
									"probes", verbose = verbose)
		slopes <- slopes[slopes$probe %in% probes, ]
		slopes$idprobe <- droplevels(slopes$idprobe)
	}

	p <- ggplot2::ggplot(data = slopes, ggplot2::aes(x = date_time))
	p <- p + ggplot2::geom_line(ggplot2::aes(y = slope_cor, col = 'slope_cor'))
	p <- p + ggplot2::geom_point(ggplot2::aes(y = slope_cor, col = 'slope_cor'))


	if (r2) {
		aux <- range(slopes$slope_cor, na.rm = TRUE)
		# artificially raise centre so R2 tends to stay about the slopes
		slope_centre <- aux[1] + (((aux[2] - aux[1]) / 2) * 1.3)
		slope_range <- aux[2] - aux[1]
		r2_centre <- 0.5
		r2_range <- 1

		if (slope_range == 0 | r2_range == 0)
			compress_factor <- 1
		else
			compress_factor <- slope_range/r2_range

		data.link <- function(x, a = r2_centre, b = slope_centre,
							  c_factor = compress_factor) {
			b + ((x - a) * c_factor) # = y
		}

		axis_link <- function(y, a = r2_centre, b = slope_centre,
							  c_factor = compress_factor) {
		  (y + a * c_factor - b) / c_factor # = x
		}

		p <- p + ggplot2::geom_hline(yintercept = data.link(r2_threshold))
		p <- p + ggplot2::geom_line(ggplot2::aes(y = data.link(r2)), 
									col = "grey")
		p <- p + ggplot2::geom_point(ggplot2::aes(y = data.link(r2),
									col = valid))
 		p <- p + ggplot2::scale_y_continuous(sec.axis = 
 			ggplot2::sec_axis(trans = axis_link, name = "R2")
 		)
		if (any(!slopes$valid)) {
			p <- p + ggplot2::scale_colour_manual(
				values = c("Red", "Black", "Grey"))
		} else {
			p <- p + ggplot2::scale_colour_manual(values = c("Black", "Grey"))
		}
	}
	
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::theme(legend.position = "none") 
	p <- p + ggplot2::labs(y = units::make_unit_label("Slope", slope_unit), 
												 x = "Time")
	p <- p + ggplot2::facet_wrap(idprobe ~ ., ncol = 1)
	return(p)
}

#' Plot the metabolic rates
#' 
#' @param input An experiment list with calculated metabolic rates.
#' 	The output of \code{\link{process_experiment}} or
#' 	\code{\link{calc_mr}}.
#' @inheritParams plot_meas
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_mr <- function(input, cycles, probes, verbose = TRUE) {
	# ggplot variables
	date_time <- mr_cor <- value <- Method <- NULL

	mr <- input$mr
	mr$idprobe <- paste0(mr$id, " (", mr$probe, ")")
	idprobe_levels <-paste0(input$probe_info$id, 
							" (", input$probe_info$probe, ")")
	mr$idprobe <- factor(mr$idprobe, levels = idprobe_levels)

	if (!missing(cycles)) {
		cycles <- check_arg_in_data(cycles, mr$cycle, 
									"cycles", verbose = verbose)
		mr <- mr[mr$cycle %in% cycles, ]
	}

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, mr$probe,
									"probes", verbose = verbose)
		mr <- mr[mr$probe %in% probes, ]
		mr$idprobe <- droplevels(mr$idprobe)
	}

	p <- ggplot2::ggplot()
	p <- p + ggplot2::geom_line(data = mr, 
								ggplot2::aes(x = date_time, y = mr_cor))
	p <- p + ggplot2::geom_point(data = mr, 
								 ggplot2::aes(x = date_time, y = mr_cor))

	if (!is.null(input$smr)) {

		mr_cols <- grepl("mr_cor", colnames(input$smr))
		smr <- reshape2::melt(input$smr, 
							  id.vars = c("probe", "id"),
							  measure.vars = colnames(input$smr)[mr_cols])
		smr$Method <- sub("_mr_cor", "", smr$variable)
		smr$idprobe <- paste0(smr$id, " (", smr$probe, ")")
		smr$idprobe <- factor(smr$idprobe, levels = idprobe_levels)

		if (!missing(probes)) {
			smr <- smr[smr$probe %in% probes, ]
			smr$idprobe <- droplevels(smr$idprobe)
		}

		p <- p + ggplot2::geom_hline(data = smr, 
									 ggplot2::aes(yintercept = value, 
					 							  linetype = Method),
									 col = "red")
	}
	
	if (!is.null(input$mmr)) {
		mmr <- input$mmr
		mmr$idprobe <- paste0(mmr$id, " (", mmr$probe, ")")
		mmr$idprobe <- factor(mmr$idprobe, levels = idprobe_levels)

		if (!missing(probes)) {			
			mmr <- mmr[mmr$probe %in% probes, ]
			mmr$idprobe <- droplevels(mmr$idprobe)
		}

		p <- p + ggplot2::geom_point(data = mmr, 
									 ggplot2::aes(x = date_time, 
						 						  y = mr_cor),
									 col = "red", size = 2)
	}
	
	p <- p + ggplot2::theme_bw()
	p <- p + ggplot2::facet_wrap(idprobe ~ ., ncol = 1)
	p <- p + ggplot2::labs(y = "\u1E40[O[2]]", x = '', linetype = "SMR method")
	return(p)	
}

#' Plot the standard metabolic rates
#' 
#' @param input An experiment list with calculated SMR.
#' 	The output of \code{\link{process_experiment}} or
#' 	\code{\link{calc_smr}}.
#' @param probes A string of which probes to plot
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_smr <- function(input, probes) {
	# ggplot variables
	mr_cor <- value <- Method <- NULL

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, input$smr$probe, "probes")
		input$smr <- input$smr[input$smr$probe %in% probes, ]
		input$mr <- input$mr[input$mr$probe %in% probes, ]
	}

	mr_cols <- grepl("mr_cor", colnames(input$smr))
	aux <- reshape2::melt(input$smr, 
						  id.vars = c("probe", "id"),
						  measure.vars = colnames(input$smr)[mr_cols])
	aux$Method <- sub("_mr_cor", "", aux$variable)

	# change facet labels
	input$mr$idprobe <- paste0(input$mr$id, " (", input$mr$probe, ")")
	idprobe_levels <- paste0(input$probe_info$id, 
							 " (", input$probe_info$probe, ")")
	input$mr$idprobe <- factor(input$mr$idprobe, levels = idprobe_levels)

	aux$idprobe <- paste0(aux$id, " (", aux$probe, ")")
	aux$idprobe <- factor(aux$idprobe, levels = idprobe_levels)

	p <- ggplot2::ggplot(data = input$mr)
	p <- p + ggplot2::geom_histogram(ggplot2::aes(x = mr_cor), 
		fill = "white", colour = "grey", binwidth = 0.1)
	p <- p + ggplot2::geom_vline(data = aux, ggplot2::aes(xintercept = value,
														  linetype = Method))
	p <- p + ggplot2::labs(x = "\u1E40[O[2]]")
	p <- p + ggplot2::facet_grid(idprobe ~ .)
	p
}


#' Plot the standard metabolic rates
#' 
#' @param input An experiment list with calculated rolling MR values.
#' 	The output of \code{\link{roll_mr}}.
#' @inheritParams plot_meas
#' 
#' @return a ggplot object
#' 
#' @export
#' 
plot_rolling_mr <- function(input, probes, cycles) {
	# ggplot variables
	phase_time <- date_time <- mr_cor <- r2 <- NULL

	mr <- input$rolling_mr$values
	max_mr <- input$rolling_mr$max

	if (is.null(mr)) {
		stop("Could not find rolling mr values in input")
	}

	if (!missing(probes)) {
		probes <- check_arg_in_data(probes, mr$probe,
									"probes", verbose = TRUE)
		mr <- mr[mr$probe %in% probes, ]
		max_mr <- max_mr[max_mr$probe %in% probes, ]
	}

	if (!missing(cycles)) {
		cycles <- check_arg_in_data(cycles, mr$cycle,
									"cycles", verbose = TRUE)
		mr <- mr[mr$cycle %in% cycles, ]
		max_mr <- max_mr[max_mr$cycle %in% cycles, ]
	}

	# remove units
	o2_unit <- units(mr$mr_cor)
	mr$mr_cor <- as.numeric(mr$mr_cor)
	max_mr$mr_cor <- as.numeric(max_mr$mr_cor)

	# make facets
	aux <- unique(mr[, c("id", "probe", "cycle")])
	aux <- aux[order(aux$probe), ]

	ipp_levels <- paste0(aux$id, " (", 
						 aux$probe, ") - ", 
						 aux$cycle)

	mr$idprobecycle <- paste0(mr$id, " (", mr$probe, ") - ", mr$cycle)
	mr$idprobecycle <- factor(mr$idprobecycle, levels = ipp_levels)
	max_mr$idprobecycle <- paste0(max_mr$id,
								  " (", max_mr$probe, ") - ", max_mr$cycle)
	max_mr$idprobecycle <- factor(max_mr$idprobecycle, levels = ipp_levels)

	
	p <- ggplot2::ggplot(data = mr, ggplot2::aes(x = phase_time, y = mr_cor))
	p <- p + ggplot2::geom_line()
	p <- p + ggplot2::geom_point(data = max_mr, col = 'red')
	p <- p + ggplot2::theme_bw()


	rer_mean <- as.numeric(mean(mr$mr_cor, na.rm = TRUE))
	r2_mean <- 0.5
	mean_dif <- r2_mean - rer_mean
	aux <- as.numeric(range(mr$mr_cor, na.rm = TRUE))
	rer_range <- aux[2] - aux[1]
	r2_range <- 1

	if (rer_range == 0 | r2_range == 0) {
		compress_factor <- 1
	} else {
		compress_factor <- rer_range/r2_range
	}

	data_link <- function(x, a = r2_mean, b = rer_mean,
						  c_factor = compress_factor) {
		b + ((x - a) * c_factor) # = y
	}

	axis_link <- function(y, a = r2_mean, b = rer_mean, 
						  c_factor = compress_factor) {
	  (y + a * c_factor - b) / c_factor # = x
	}

	p <- p + ggplot2::geom_line(ggplot2::aes(y = data_link(r2)), col = "grey")
	p <- p + ggplot2::scale_y_continuous(
		sec.axis = ggplot2::sec_axis(trans = axis_link, name = "R2")
	)

	p <- p + ggplot2::labs(y = units::make_unit_label("\u1E40[O[2]]", o2_unit),
						   x = "Time within phase")
	p <- p + ggplot2::facet_wrap(idprobecycle ~ ., ncol = 1)
	p
}


#' Wrapper to plot the whole experiment data for one probe
#' 
#' @param input An experiment list with the experimental data.
#' The output of \code{\link{process_experiment}}.
#' @param probe The probe to plot (can only plot one probe at a time).
#' @inheritParams plot_meas
#' 
#' @export
#' 
plot_experiment <- function(input, cycles, probe, verbose = TRUE) {

	if (length(probe) > 1) {
		stop("Please select only one probe at a time.")
	}
	if (!is.null(input$bg$pre)) {
		# if (verbose) message('Plotting pre-background')
		p_pre <- plot_bg(input$bg$pre, probes = probe) 
		p_pre <- p_pre + ggplot2::labs(title = 'Pre-background')
	} else {
		p_pre <- patchwork::wrap_elements(
			grid::textGrob('Pre-bg was not included'))
	}

	if (!is.null(input$bg$post)) {
		# if (verbose) message('Plotting post-background')
		p_post <- plot_bg(input$bg$post, probes = probe) 
		p_post <- p_post + ggplot2::labs(title = 'Post-background')
	} else {
		p_post <- patchwork::wrap_elements(
			grid::textGrob('Post-bg was not included'))
	}

	# if (verbose) message('Plotting measurements')

	p_meas <- plot_meas(input, cycles = cycles, probes = probe,
						show_temp = TRUE, verbose = verbose)
	p_meas <- p_meas + ggplot2::labs(x = '')
	
	# if (verbose) message('Plotting deltas')

	p_delta <- plot_deltas(input, cycles = cycles,
			 			   probes = probe, verbose = FALSE)
	p_delta <- p_delta + mimic_x_datetime(p_meas)

	# if (verbose) message('Plotting plotting slopes')

	p_slopes <- plot_slopes(input, cycles = cycles,
						    probes = probe, verbose = FALSE) 
	p_slopes <- p_slopes + mimic_x_datetime(p_meas) + ggplot2::xlab('')

	# if (verbose) message('Plotting metabolic rate')

	probe_check <- any(input$mr$probe == probe)
	cycle_check <- (missing(cycles) || any(input$mr$cycle %in% cycles))
	combined_check <- (missing(cycles) || 
					  any(input$mr$probe == probe & input$mr$cycle %in% cycles))
	if (probe_check && cycle_check && combined_check) {
		p_mr <- plot_mr(input, cycles = cycles, probes = probe, verbose = FALSE)
		p_mr <- p_mr + mimic_x_datetime(p_meas)
	} else {
		p_mr <- patchwork::wrap_elements(
			grid::textGrob('No valid MO2 values found'))
	}
	p_final <- p_pre + p_post + p_meas + p_delta + p_slopes + 
		 	   p_mr + patchwork::plot_layout(design = 'AB\nCC\nDD\nEE\nFF')

	return(p_final)
}



mimic_x_datetime <- function(input) {
	ggplot2::xlim(as.POSIXct(ggplot2::layer_scales(input)$x$range$range,
	 			  origin = "1970-01-01 00:00.00"))	
}

mimic_y_datetime <- function(input) {
	ggplot2::ylim(as.POSIXct(ggplot2::layer_scales(input)$y$range$range,
	 			  origin = "1970-01-01 00:00.00"))	
}

mimic_x <- function(input) {
	ggplot2::xlim(ggplot2::layer_scales(input)$x$range$range)
}

mimic_y <- function(input) {
	ggplot2::ylim(ggplot2::layer_scales(input)$y$range$range)
}


