source("./functions.r")
QC_normalized_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    downloadButton(ns('download_density'), 'Download plot'),
    plotOutput(ns("density_plot"),
                  brush = brushOpts(id = ns("dense_brush"),resetOnNew = TRUE),
                  dblclick = ns("dense_dblclick")),
    downloadButton(ns('download_ecdf'), 'Download plot'),
    plotOutput(ns("ecdf_plot"))
  )
}
QC_normalized<-function(input,output,session,normal,zoom)
{# Single zoomable plot (on left)
  dense_ranges <- reactiveValues(x = c(0,1000), y = c(0,0.035))
  
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Preparing plots", value = 0)
  
  # Number of times we'll go through the loop to update progress bar
  n <- 3
  
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Step 1: get density plot
  output$density_plot<-renderPlot({
  density_plot(normal,dense_ranges$x,dense_ranges$y,zoom)
  })
  
  observeEvent(input$dense_dblclick,{
    brush <- input$dense_brush
    if (!is.null(brush)) {
      
      dense_ranges$x <- c(brush$xmin, brush$xmax)
      dense_ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      dense_ranges$x<-c(0,1000)
      dense_ranges$y<-  c(0,0.035)
    }
  })
  
  #download density plot
  output$download_density <- downloadHandler(
    filename = paste(" Density plot for normalized data.png"),
    content = function(file) {
      png(file)
      density_plot(normal,dense_ranges$x,dense_ranges$y,zoom)
      dev.off()
    }) 
  
  # Increment the progress bar after density plot is obtained, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Step 2: Get ECDF plot
  output$ecdf_plot<-renderPlot(
  {
    ECDF_plot(normal)
  })
  #download ecdf plot
  output$download_ecdf <- downloadHandler(
    filename = paste("ECDF plot of normalized data.png"),
    content = function(file) {
      png(file)
      ECDF_plot(normal)
      dev.off()
    }) 
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
}