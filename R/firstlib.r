.onAttach <- function (lib, pkg) {
  ver <- read.dcf(file.path(lib,pkg,"DESCRIPTION"),"Version")
  ver <- as.character(ver)
  packageStartupMessage("\nBDP2 ",ver," loaded.\n", 
                        "Use BDP2workflow() for interactive shiny app, see 'help(\"BDP2-package\")' and vignette for more information.\n",
                      # "Please cite as:\n   ",format(citation("BDP2"), style = "text") ,
                         domain = NULL,  appendLF = TRUE)
}

.onLoad <- function(...) {
}
