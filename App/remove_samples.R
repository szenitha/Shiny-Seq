source("./functions.r")
#######remove samples###
Remove_samples_UI<-function(id)
{
  ns<-NS(id)
  tagList(
          fluidRow(column(6,h3("PCA"))),
          fluidRow(column(10,plotlyOutput(ns("pca_rem")))),
         bscols(
        selectInput(ns("top_rem"), label = h5("Should PCA be computed for:"),
                   choices = list("Top 500 genes" = 1, "All" = 2),
                   selected = 2),
          column(6,uiOutput(ns("annotation_rem"))),
          column(6,textOutput(ns("box1")))),
        fluidRow(column(6,h3("Boxplot"))),
        fluidRow(column(10,plotOutput(ns("boxPlot_rem"),click =ns("clickBar"),dblclick=ns("deselect"))))#,
    # bscols(widths = 10,div(style="height: 95px;",width = '200px')),
    #bscols(column(3,actionButton(ns("del"), "Delete"),actionButton(ns("reset"), "Reset")))#,
    # bscols(
    #   column(10,DT::dataTableOutput(ns("outliers"))))
           
                                        
  )
  
}

# creates a checkbox widget to select the annotation for PCA(for normalized data)
Remove_samples<-function(input,output,session,dds.fc,rld,normal)#edata=NULL,pData=NULL,
{
  print("inside remove samples line 30 ")
  v_data<-NULL
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Processing outlier tab", value = 0)
  
  # Number of times we'll go through the loop to update progress bar
  n <- 6
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Step 1: Set up a list of annotation variables to color samples in PCA plot for visualization
output$annotation_rem <-
    renderUI({
      #if (!is.null(input$file2) && input$ok2>0) {
      #print("create Checkbox")
      checklist = list()
      for (i in seq_along(colnames(colData(dds.fc)))) {
        checklist[[colnames(colData(dds.fc))[[i]]]] = i
      }
      #radioButtons("an_pca_norm", "Choose the annotation", checklist)
      div(style="height: 95px;",radioButtons(session$ns("annotation_rem"), "Choose the annotation", checklist,inline=TRUE,width = '200px'))
      #}
    })

# Increment the progress bar once above step is successful, and update the detail text.
progress$inc(1/n, detail = paste("Doing part", 2,"/",n))

# Pause for 0.1 seconds to simulate a long computation.
Sys.sleep(0.1)

  #Step 2: Display PCA for unnormalized data
  
  output$pca_rem<-renderPlotly({
 
      req(input$annotation_rem)
    
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Making PCA plot", value = 0)
      
      # Number of times we'll go through the loop to update progress bar
      n1 <- 3
      
      #Get expression table and annotation table
      data<-rld#normal()
      
      rv <- rowVars(assay(data))
      genes<-order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]
      pheno<-colData(dds.fc)
      
      #Get which column u want to plot by
      var<-as.numeric(input$annotation_rem)
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n1, detail = paste("Doing part", 1,"/",n1))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      parameters<-list()
      if(input$top_rem==2) parameters<-pca_components(assay(data),var,pheno)#pca of all genes
      else parameters<-pca_components(assay(data)[genes, ],var,pheno)#pca of top 500 most variable genes
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n1, detail = paste("Doing part", 2,"/",n1))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      #Define colors
      library(RColorBrewer)
      n <- 60
      qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      
      #print(parameters)
      xlab = paste("PC1:",parameters[[4]],"%")
      ylab = paste("PC2:",parameters[[5]],"%")
      zlab = paste("PC3:",parameters[[6]],"%")
      df<-data.frame(PCA1=parameters[[1]],PCA2=parameters[[2]],PCA3=parameters[[3]],an=parameters[[7]])
      #print(head(df))
      p<-  plot_ly(df,x=~PCA1,y=~PCA2,z=~PCA3,key=rownames(df),source="B") %>%
        add_markers(type = "scatter3d",color=df$an,colors = unique(colors[df$an]))%>%
        layout(scene = list(xaxis = list(title = xlab),
                            yaxis = list(title = ylab),
                            zaxis = list(title = zlab)))#%>%
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n1, detail = paste("Doing part", 3,"/",n1)) 
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      p
    
  })
  # Increment the progress baronce above step is succesful, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
# Step 3:Display interactive boxplot
 output$boxPlot_rem <- renderPlot({

      req(input$annotation_rem)
      boxplot_output(normal,colData(dds.fc),input$annotation_rem)

  })
 # Increment the progress bar, and update the detail text.
 progress$inc(1/n, detail = paste("Doing part", 4,"/",n))
 
 # Pause for 0.1 seconds to simulate a long computation.
 Sys.sleep(0.1)
 
 #Step 4: display sample ID of selected boxplot
output$box1 <- renderText({
    if (is.null(input$clickBar)) return("")
    else {
      pheno<-colData(dds.fc)
      lvls <- pheno[,1]
      #print(lvls)
      #print(input$clickBar)
      # print(input$clickBar$x)
      # print(round(input$clickBar$x))
      name <- lvls[round(input$clickBar$x)]
      paste("You have selected",name)
    }

  })
# Increment the progress bar, and update the detail text.
progress$inc(1/n, detail = paste("Doing part", 5,"/",n))

# Pause for 0.1 seconds to simulate a long computation.
Sys.sleep(0.1)
# 
v1 <- reactiveValues(
    clickBar_list = NULL,
    clickpca_list=NULL)#Represents the first mouse click, if any

#   
#Step 6: Handle clicks on the plot 
#selection
observe({
    input$clickBar
    isolate({
      #print('hey')
      #print(round(input$clickBar$x))
      v1$clickBar_list = c(v1$clickBar_list, input$clickBar$x)
      #print('ok')
    })
  })
# deselection
observe({
    input$deselect
    isolate({
      #print('hey')
      #print(round(input$clickBar$x))
      if(!is.null(v1$clickBar_list))
      {
        print(input$deselect$x)
        idx<-which(v1$clickBar_list %in% input$deselect$x)
        v1$clickBar_list = v1$clickBar_list[-idx]
        #print('ok')
      }


    })
  })
# 
# Handling pca  plot clicks
observe({

    event_data("plotly_click",source="B")
    d<-event_data("plotly_click",source = "B")
    isolate({
      #print('hey')
      #print(round(input$clickBar$x))
      v1$clickpca_list = c(v1$clickpca_list, d$key)
      #print('ok')
    })
  })
return(list(d=reactive({v1$clickpca_list}),bar=reactive({v1$clickBar_list})))

# Increment the progress bar afterprocessing is completed, and update the detail text.
progress$inc(1/n, detail = paste("Doing part", 6,"/",n))

# Pause for 0.1 seconds to simulate a long computation.
Sys.sleep(0.1)
}           






