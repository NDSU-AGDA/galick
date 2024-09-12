library(shiny)

ui <- fluidPage(
  fileInput(
    "vcf",
    "Tu vcf"
  ),
  downloadButton("descargar")
)

server <- function(input, output, session) {
  
  output$descargar <- downloadHandler(
    filename = "imputed.vcf.gz",
    content = function(file) {
      if (is.null(input$vcf$datapath)) {
        showNotification("No file uploaded")
        return()
      }
      beagler::run_beagle(gt = input$vcf$datapath, out = file, n_threads = 1)
    }
  )

}

shinyApp(ui, server)
