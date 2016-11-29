## Run Shiny App

runWebInterface <- function() {
    appDir <- system.file("webInterface", package = "supern")
    if (appDir == "") {
        stop("Could not find web app directory..", call. = FALSE)
    }
    
    shiny::runApp(appDir, display.mode = "normal")
}

