source("./functions.R")
#Module PCA
#UI component
#Server component: Calls Pcaplot function from fucntions.R

PCA_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    #User input: PCA type: 2D/3D
    radioButtons(ns("plotType"), "PCA Type:", choices = c("2D", "3D")),
    #User input: Annotation variables to color samples in PCA
    uiOutput(ns("annotation_pca")),
    #User input: Compute PCA for top 500 most variable genes/All genes
    selectInput(ns("top"), label = h5("Should PCA be computed for:"), 
                choices = list("Top 500 genes" = 1, "All" = 2),
                selected = 2),
    #Download button: PCA gets downloaded on user click
    fluidRow(column(3," "),column(6,downloadButton(ns('download_pca_svg'), 'Download plot as svg'))),
    #Output: Display PCA plot
    fluidRow(column(1," "),
             column(6,plotlyOutput(ns("pca"))))
  )
}
PCA<-function(input,output,session,data,resize_factor=NULL)
{
     #Creates a checkbox widget containing annotation variables to color samples in PCA
     output$annotation_pca<-renderUI({
                            #if (!is.null(input$file2)) {
                            print("create Checkbox")
                            checklist = list()
                            #print(head(data[[3]]))
                            for (i in seq_along(colnames(data[[2]]))) { 
                                  checklist[[colnames(data[[2]])[[i]]]] = i
                                 }
                             radioButtons(session$ns("annotation_pca"), "Choose the annotation", checklist)
    #}
                           })
  
     #Generate PCA plot
      output$pca<-renderPlotly({
  
                  pcaplot(data,input$annotation_pca,input$top,input$plotType,resize_factor)[[1]]
  
                  })
    #Download pca in svg format
      output$download_pca_svg <- downloadHandler(
  
                                 filename = function(){
                                   if (identical(input$plotType, "2D"))
                                   {
                                     if(input$top==2)  paste('2D PCA of Raw data all genes.svg')
                                     else  paste('2D PCA of Raw data top 500 genes.svg')
                                   }
                                   
                                   else if(identical(input$plotType, "3D"))
                                   {
                                     if(input$top==2)  paste('3D PCA of Raw data all genes.pdf')
                                     else  paste('3D PCA of Raw data top 500 genes.pdf')
                                   }
                                 },
                                 content = function(file) {
                                   
                                   p<-pcaplot(data,input$annotation_pca,input$top,input$plotType)[[2]]
                                   ggsave(file,p)
                                 })
      return(reactive({input$top}))
}
