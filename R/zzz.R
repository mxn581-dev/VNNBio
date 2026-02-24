.onAttach <- function(libname, pkgname) {
    ## Do NOT call juliaSetupOk() here -- it can open a socket connection
    ## that R CMD check flags as "connections left open".
    ## Julia availability is checked lazily when initJuliaBackend() is called.
    invisible(NULL)
}
