plot_o2 <- function(folder) {

	files <- list.files(paste0(folder, '/ChannelData/'))

	O2.file.link <- grepl('Oxygen.txt', files)

	if (all(!O2.file.link)) {
		stop('No oxygen files found')
	}

	O2.files <- files[O2.file.link]

	aux <- lapply(O2.files, function(i) {
		x <- load.pyroscience.o2.file(paste0(folder, '/ChannelData/', i), date.format = '%d-%m-%Y')
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