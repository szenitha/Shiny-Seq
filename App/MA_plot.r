source("./functions.r")
MA_plot_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    uiOutput(ns("ma_comb")),
    sliderInput(ns("scale"),"y-limit",1,15,value = 4,step = 1),
    plotOutput(ns("MA_plot")),
    downloadButton(ns('download_ma_plot'), 'Download plot')
  )
}
MA_plot<-function(input,output,session,combination,DE_genes,p_values)
{
  combo<-combination()
  #combinations: Display comparisons A vs B , C Vs D etc
  output$ma_comb <- renderUI({
    if(!is.null(combo()))
    {
      print("inside ma plot module line 16")
      #num <- length(combination())
      num<- length(combo())
      comb<-lapply(1:num, function(i) {
        combo()[[i]]
        #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
        
      })
      
      checklist<-list()
      for (i in seq_along(comb)) {
        checklist[[comb[[i]]]] = i
      }
      selectInput(session$ns("ma_choice"),label = h5("Choose comparison") ,
                  choices = checklist,selected = 1)
    }
  })
  
  #ma plot for the comparison selected
  output$MA_plot<- renderPlot({
    
    print(p_values())
    print("inside map plot line 40")
    print(combo()[as.numeric(input$ma_choice)])
    ma_plot(input$ma_choice,combo(),DE_genes(),input$scale,p_values())
    
    
  })
  
  #download MA plot
  output$download_ma_plot <- downloadHandler(
    filename = function(){paste("MA plot",combination()[as.numeric(input$ma_choice)],".pdf")},
    content = function(file) {
      p<-ma_plot(input$ma_choice,combination(),DE_genes(),input$scale,p_values())
      ggsave(filename=file, plot=p)
      
    })
  # 
return(list(
  input_scale=reactive({input$scale}),
  p_values=reactive({p_values()})
))
}




