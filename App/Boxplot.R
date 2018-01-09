source("./functions.r")
Boxplot_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    downloadButton(ns('download_boxplot'), 'Download plot'),
    plotOutput(ns("boxplot")),
    fluidRow(column (7,uiOutput(ns("annotation"))))
  )
}
Boxplot_module<-function(input,output,session,data)
{
  output$annotation <-
    renderUI({
      #if (!is.null(input$file2)) {
        print("create Checkbox boxplot")
        checklist = list()
        for (i in seq_along(colnames(data[[2]]))) {
          checklist[[colnames(data[[2]])[[i]]]] = i
        }
        #print(checklist[[1]])
        radioButtons(session$ns("annotation"), "Choose the annotation", checklist)
     # }
    })
  
  output$boxplot<-renderPlot({
    print(head(data[[1]]))
             req(input$annotation)
             boxplot_output(data[[1]],data[[2]],input$annotation)
                 })
  
  output$download_boxplot <- downloadHandler(
    
    filename = function(){
      paste('Boxplot of Raw data .png')
    },
    content = function(file) {
      png(file)
      boxplot_output(data[[1]],data[[2]],input$annotation)
      dev.off()
      
    })  
  return(reactive({req(input$annotation)
    boxplot_output(data[[1]],data[[2]],input$annotation)}))
}