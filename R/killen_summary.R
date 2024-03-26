#' Compile a list of points following Killen et al. 2021
#' 
#' For more information on these points, see Killen et al. 2021:
#' 
#' @references Killen, S. S., Christensen, E. A. F., Cortese, D., Závorka, L.,
#' 	Norin, T., Cotgrove, L., Crespel, A., Munson, A., Nati, J. J. H., 
#' 	Papatheodoulou, M., & McKenzie, D. J. (2021). Guidelines for reporting 
#' 	methods to estimate metabolic rates by aquatic intermittent-flow 
#' 	respirometry. Journal of Experimental Biology, 224(18), jeb242522. 
#' 	<https://doi.org/10.1242/jeb.242522>
#' 
#' @export
#' 
#' @return a list of summary information
#' 
killen_summary <- function(input) {
	# 04
	pt04 <- with(input$probe_info, (volume - conv_w_to_ml(mass)) / conv_w_to_ml(mass))
	# 15
	recipient <- lapply(input$pyro$source_data, function(i) {
		prb <- data.frame(probe = paste0(attributes(i)$device, attributes(i)$ch))
		cal <- t(as.data.frame(attributes(i)$cal_settings))
		return(cbind(prb, cal))
	})

	pt15 <- do.call(rbind, recipient)
	pt15 <- type.convert(pt15, as.is = TRUE)
	date_cols <- grepl("lastCal", colnames(pt15))
	if (any(date_cols)) {
		exp_day <- as.Date(min(input$pyro$compiled_data$date_time))
		for(i in which(date_cols)) {
			pt15[, i] <- as.Date(pt15[, i])
			pt15[, paste0("cal", i, "_to_exp_start")] <- 
				difftime(exp_day, pt15[, i], units = "days")

		}
	}
	rownames(pt15) <- 1:nrow(pt15)
	#17
	sem_aux <- sem(input$cleaned$temp, na.rm = TRUE)
	units(sem_aux) <- units(input$cleaned$temp)
	pt17 <- c(mean = mean(input$cleaned$temp, na.rm = TRUE),
		     sem = sem_aux)

	# 30
	pt30 <- sapply(input$bg, function(i) {
				length(unique(i$cleaned$phase))
		   })

	# 36
	pt36 <- gsub("_mr_cor", "", colnames(input$input)[grepl("_mr_cor", colnames(input$input))])

	# 39
	by_probe <- table(input$all_slopes$probe)
	by_probe <- as.data.frame(by_probe)
	colnames(by_probe) <- c("probe", "n_all")
	by_probe$n_good <- as.data.frame(table(input$good_slopes$probe))$Freq
	by_probe$n_discarded <- by_probe$n_all - by_probe$n_good
	by_probe$pct_discarded <- by_probe$n_discarded / by_probe$n_all * 100
	units(by_probe$pct_discarded) <- "percent"

	run_total <- as.data.frame(apply(by_probe[, -1], 2, sum, simplify = FALSE))
	run_total$pct_discarded <- run_total$n_discarded / run_total$n_all * 100
	units(run_total$pct_discarded) <- "percent"

	pt39 <- list(by_probe = by_probe,
			   run_total = run_total)

	# bring it all together
	output <- list(
		pt01 = input$probe_info$mass,
	    pt02 = input$probe_info$volume,
		pt04 = pt04,
		pt14 = attributes(input$cleaned)$wait_time,
		pt15 = pt15,
		pt17 = pt17,
		pt22 = min(input$cleaned$airsat, na.rm = TRUE),
		pt24 = nrow(input$probe_info),
		pt30 = pt30,
		pt31 = attributes(input$cleaned)$correction_method,
		pt35 = difftime(max(input$cleaned$date_time), min(input$cleaned$date_time)),
		pt36 = pt36,
		pt37 = table(input$good_slopes$probe),
		pt38 = attributes(input$good_slopes)$r2_threshold,
		pt39 = pt39,
		pt51 = paste0("Oxygen uptake rates were calculated using ",
		   		    "the R package pyroresp (version ",
		   		    packageVersion("pyroresp"), "; Flávio, 2024)."),
		pt52 = "Yes"
		)
	return(output)
}
