#' Convert a human-friendly respirometry table to a computer-friendly format
#'
#' Also appends additional probe info if provided.
#'
#' @param input A resp list containing a phased object.
#'
#' @export
#'
melt_resp <- function(input) {

  if (is.null(input$phased)) {
    stop("Couldn't find object 'phased' inside input.",
         " Have you run assign_phases?", call. = FALSE)
  }
  phased <- input$phased

  ox_aux <- reshape2::melt(phased,
    id.vars = "date_time",
    measure.vars = colnames(phased)[grepl("ox_", colnames(phased))],
    value.name = "o2",
  )
  ox_aux$probe <- stringr::str_extract(ox_aux$variable, "(?<=ox_).*")
  ox_aux$variable <- NULL
  colnames(ox_aux)[1] <- "date_time"
  ox_aux <- ox_aux[, c("date_time", "probe", "o2")]

  pressure_aux <- reshape2::melt(phased,
    id.vars = "date_time",
    measure.vars = colnames(phased)[grepl("pressure_", colnames(phased))],
    value.name = "pressure"
  )
  pressure_aux$variable <- NULL
  pressure_aux$date_time <- NULL

  temp_aux <- reshape2::melt(phased,
    id.vars = "date_time",
    measure.vars = colnames(phased)[grepl("temp_", colnames(phased))],
    value.name = "temp"
  )
  temp_aux$variable <- NULL
  temp_aux$date_time <- NULL

  sal_aux <- reshape2::melt(phased,
    id.vars = "date_time",
    measure.vars = colnames(phased)[grepl("sal_", colnames(phased))],
    value.name = "sal"
  )
  sal_aux$variable <- NULL
  sal_aux$date_time <- NULL

  phase_aux <- reshape2::melt(phased,
    id.vars = "date_time",
    measure.vars = colnames(phased)[grepl("phase_", colnames(phased))],
    value.name = "phase"
  )
  phase_aux$variable <- NULL
  phase_aux$date_time <- NULL

  pre_output <- cbind(ox_aux, pressure_aux, temp_aux, sal_aux, phase_aux)

  if (any(grepl("ph_", colnames(phased)))) {
    pre_output$ph <- as.vector(t(phased[, grepl("ph_", colnames(phased))]))
  }

  pre_output$cycle <- as.numeric(sub("F|M", "", as.character(pre_output$phase)))

  # if fish information is provided
  if (!is.null(input$probe_info)) {
    # Include the extra information
    link <- match(pre_output$probe, input$probe_info$probe)
    output <- cbind(pre_output,
                    input$probe_info[link, c("id", "mass", "volume")])
  } else {
    output <- as.data.frame(pre_output)
  }

  output$phase <- as.character(output$phase)

  input$melted <- output
  return(input)
}
