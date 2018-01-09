gene_count_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    uiOutput(ns("plot_option_gene_n")),
    downloadButton(ns('download_gene_plot'), 'Download plot'),
    plotOutput(ns("gene_plot"))
  )
}

gene_count_module<-function(input,output,session,dataset,dds.fc,gene_name)
{
  output$plot_option_gene_n<-renderUI({
    
    radioButtons(session$ns("plotType_gene_n"), "use log scale?:", choices = c("Yes", "No"),selected = "Yes")
    
  })
  
  geneplot<-reactive({
    req(input$plotType_gene_n)
    gene_counts(dataset,dds.fc(),gene_name(),input$plotType_gene_n)
  })
  
  output$gene_plot<-renderPlot({
    geneplot()
    
  })
  #download boxplot for batch corrected table
  output$download_gene_plot <- downloadHandler(
    filename = paste(" Gene counts of ",gene_name()," normalized .svg"),
    content = function(file) {
      ggsave(file, geneplot())
    }) 
  
  
}
