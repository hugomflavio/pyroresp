# Startup message
#
# Startup text message linked to the publication describing the 'FishResp' package
#' @importFrom utils packageDescription

.onAttach<- function (lib, pkg){
  if(interactive()){
    pkg.version = packageDescription("FishResp", fields = "Version")
    startup.txt = paste("
THIS IS A FORKED VERSION OF THE FishResp PACKAGE!\n
Do not use this repository if you are looking for the real FishResp!\n
Go here instead: https://github.com/embedded-sergey/FishResp-Rpackage\n
")
    packageStartupMessage(startup.txt)
  }
}
