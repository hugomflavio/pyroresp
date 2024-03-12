.onAttach <- function(libname, pkgname){
  inst.ver <- utils::packageVersion("pyroresp")
  aux <- as.matrix(data.frame(Package = "pyroresp", LibPath = NA, Version = as.character(utils::packageVersion("pyroresp")), Priority = NA, Built = NA))
  packageStartupMessage("Welcome to pyroresp (", inst.ver, ")!\nRun ?pyroresp for starting tips.")
  # new.ver <- tryCatch(old.packages(instPkgs = aux, repos = "https://cloud.r-project.org"), warning = function(w) NULL, error = function(e) NULL)
  # if (!is.null(new.ver)) { # nocov start
  #   packageStartupMessage(paste0("-------------------------------------------------------------\n!!! A NEW VERSION of pyroresp is available! (v.", inst.ver, " -> v.", new.ver[, "ReposVer"], ")\n!!! You should update pyroresp before continuing.\n!!! You can update your packages by running update.packages()\n-------------------------------------------------------------\n"))
  # } # nocov end
}

if (FALSE) {
  knitr::knit()
  rmarkdown::render()
}
# These dummy lines suppress the R check note regarding unused dependencies.
# knitr and rmarkdown are needed to compile the vignettes.
