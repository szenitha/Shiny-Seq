source("./functions.r")
Volcano_plot_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    uiOutput(ns("vol_comb")),
    uiOutput(ns("vol_y")),
    uiOutput(ns("vol_x")),
    plotlyOutput(ns("volcano_plot"))
    
  )
}
Volcano_plot<-function(input,output,session,combination,DE_genes,p_values,hypothesis_choice)
{
  print("inside volcano plot line 14")
  #combinations for volcano plot (A vs B, C vs D)
  output$vol_comb <- renderUI({
    if(!is.null(combination()))
    {
      #num <- length(combination())
      num<- length(combination())
      comb<-lapply(1:num, function(i) {
        input$combination[i]
      })
      
      checklist<-list()
      for (i in seq_along(comb)) {
        checklist[[comb[[i]]]] = i
      }
      selectInput(session$ns("vol_choice"),label = h5("Choose comparison") ,
                  choices = checklist,selected = 1)
    }
  })
  
  #Slider to adjust Y axis of volcano plot
  output$vol_y<-renderUI({
    if(!is.null(input$vol_choice))
    {
      num<- length(combination())
      #Get the deseq2 dataset
      res<-DE_genes()[[as.numeric(input$vol_choice)]][5] [[1]]
      num<-(-log10(na.omit(res$padj)))
      idx<-which(is.finite(num))
      max<-max(round(num[idx],1))
      #max<-(-log10(min(na.omit(res$padj))))
      #max<-round(max,1)
      sliderInput(session$ns("scale_voly"),"y-limit",0,max+5,value = max,step = 0.5)
    }
  })
  #slider to adjust x axis of volcano plot
  
  output$vol_x<-renderUI({
    if(!is.null(input$vol_choice))
    {
      num<- length(combination())
      #Get the deseq2 dataset
      res<-DE_genes()[[as.numeric(input$vol_choice)]][5] [[1]]
      max<-max(na.omit(abs(res$log2FoldChange)))
      max<-round(max,1)
      print(max)
      sliderInput(session$ns("scale_volx"),"x-limit",0,max+5,value = max+5,step = 0.5)
    }
  })
  #volcano plot
  output$volcano_plot<- renderPlotly({
    print("inside volcano plot line 64")
    req(input$vol_choice)
    req(input$scale_volx)
    req(input$scale_voly)
    if(!is.null(input$vol_choice) & !is.null(input$scale_volx) & !is.null(input$scale_voly))
    {
      print("inside volcano plot line 69")
      volcano_plot(input$vol_choice,combination(),DE_genes(),input$scale_volx,input$scale_voly,p_values(),hypothesis_choice())[[2]]
     }
     else plotly_empty()
  })
  
return(list(
  input_scale_volx=reactive({input$scale_volx}),
  input_scale_voly=reactive({input$scale_voly}),
  hypothesis_choice=reactive({hypothesis_choice()})
))  
}