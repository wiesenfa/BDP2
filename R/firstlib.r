.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("\nBDP2 ",ver," loaded.\n"#, # Type 'help(\"futility-package\")' for an overview.
                       # "\nPlease cite as:\n   "
                         , domain = NULL,  appendLF = TRUE)
}

.onLoad <- function(...) {
}
