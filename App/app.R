source("./library.r")
source("./module1.r")
source("./PCA.r")
source("./Boxplot.r")
source("./module2.r")
source("./remove_samples.r")
source("./QC_normalized_data.r")
source("./module3.r")
source("./functions.r")
source("./Differential_expression.r")
source("./MA_plot.r")
source("./Volcano_plot.r")
source("./WGCNA.r")
source("./Kegg_module.r")
source("./Biological_Process_module.r")
source("./Hallmark_module.r")
source("./FC-FC_plot_module.r")
source("./Venn_diagram_module.r")
source("./Enriched_markers_module.R")
source("./Heatmap_module.R")
source("./ANOVA_module.r")
source("./Transcription_factor_prediction_module.r")
source("./powerpoint_module.r")
# Define UI for application that draws a histogram 
#shinyUI(fluidPage(
jsResetCode <- "shinyjs.reset = function() {history.go(0)}" # Define the js method that resets the page

ui<-navbarPage("Shiny-Seq",
               navbarMenu("Raw Data",
                       tabPanel("Raw Data",
                                 Module_Raw_data_UI("module"),
                                 useShinyjs(),# Set up shinyjs
                                 # Add a CSS class for red text colour
                                 inlineCSS(list(.blue = "background: lightblue")),
                                  uiOutput("button"),
                                 #Alter user if the number of samples in count/expression data
                                 #is not same as annotation table
                                fluidRow(column(8,
                                                bsAlert("alert_app"))),
                                 #Display annotation table
                                 downloadButton('downloadepData', 'Download full Data'),
                                 fluidRow(column (10,DT::dataTableOutput("pData"))),
                                 #Display expression table
                                 downloadButton('downloadeData', 'Download full Data'),
                                 fluidRow(column(12,DT::dataTableOutput("edata")))
                               ),

                        tabPanel("Boxplot",
                                 Boxplot_module_UI("module")
                          ),
                       tabPanel( "PCA",
                                 PCA_UI("module")
                       )
               ),
               navbarMenu("Normalization",
                          tabPanel("Normalized table",
                                    conditionalPanel(condition="input.ok1==0",
                          #                           #p(input.ok1),
                                                     p("please press start in the unnormalized table tab")),
                                    conditionalPanel(condition="input.ok1 >0",
                                                     fluidRow(column(4,uiOutput("condition"))),                 
                                                     fluidRow(column(5, uiOutput("design"))),
                                                     useShinyjs(), # Set up shinyjs
                                                    # Add a CSS class for red text colour
                                                     inlineCSS(list(.blue = "background: lightblue")),
                                                     uiOutput("button1"),
                                                    actionButton("rem_samp", "Detect outliers"),
                                                    actionButton("qc","Quality control"),
                                                    #actionButton("skip_bc","Skip Batch analysis"),
                                                    bsModal("modalqc","QC for normalized data","qc",size = "large",
                                                            conditionalPanel("input.ok2==0",
                                                                             p("Please press start button in the normalized tab and wait for a few minutes prior to clicking this button")),
                                                            conditionalPanel("input.ok2>0",
                                                                             QC_normalized_UI("module"))),
                                                    Module_Normalized_data_UI("module"),
                                                    
                                                    bsModal("modaloutlier", "Detect and remove outlier samples in addition to identifying batch effects", "rem_samp",size = "large",#uiOutput("plots")),
                                                            conditionalPanel("input.ok2==0",
                                                                             p("please press start button in the normalized tables tab and wait for a few minutes
                                                                               prior to clicking this button")
                                                                             ),
                                                            conditionalPanel("input.ok2>0",
                                                            Remove_samples_UI("module"),
                                                            bscols(column(3,actionButton("del", "Delete"),actionButton("reset", "Reset"))),
                                                            bscols(
                                                              column(10,DT::dataTableOutput("outliers")))
                                                            
                                                    ))

                         )
               ),
                          tabPanel("Boxplot",
                                   conditionalPanel(condition = "input.ok2==0",
                                                    p("Please press start button in normalized table tab")),
                                   conditionalPanel(condition="input.ok2 >0",
                                                    #p("hey"),
                                   Boxplot_module_UI("module1"))
                          ),
               tabPanel("Exploratory data anlysis",
                        tabsetPanel(
                           tabPanel("PCA",
                                    conditionalPanel(condition = "input.ok2==0",
                                                     p("Please press start button in normalized table tab")),
                                    conditionalPanel(condition = "input.ok2>0",
                                                     PCA_UI("module1"))),
                            tabPanel("Heatmap of samples",
                                     conditionalPanel(condition = "input.ok2==0",
                                                      p("Please press start button in normalized table tab")),
                                     conditionalPanel(condition = "input.ok2>0",
                                                     plotOutput("heatmap_sample")))
                        )
               )  
               ),
          tabPanel("Batch effect analysis",
                   tabsetPanel(
                     tabPanel("Compute batch effect",
                              actionButton("ok3","Start"),
                              conditionalPanel("input.ok3>0",
                                               Module_Batch_effect_UI("module"))
                              ),
                     
                     tabPanel("Source of variation",
                              conditionalPanel(condition = "input.ok2==0",
                                               p("Please press start button in normalized table tab")),
                              conditionalPanel(condition = "input.ok2>0",
                                               uiOutput("sov_data_input"),
                                               uiOutput("go_sov"),
                                               plotlyOutput("sov")
                                               )
                              )
                   )
                     ),
          tabPanel("WGCNA",
                   fluidRow(column(3,
                                   selectInput("option_wgcna", label = h5("Select choice"), 
                                               choices = list("","
                                                              Perform WGCNA independently" = 1,
                                                              "Perform as part of DE analysis" = 2),
                                               selected = NULL))),
                   uiOutput("file_input1"),
                   uiOutput("file_input2"),
                   fluidRow(column (10,DT::dataTableOutput("pheno_table"))),
                   uiOutput("go_wgcna"),
                   conditionalPanel("input.go_wgcna>0",
                                    WGCNA_module_UI("module"))
                  
                                 
                   ),
          navbarMenu("Downstream Analysis",
                       tabPanel("Differential expression analysis",
                                tabsetPanel(
                                  tabPanel("Get list of DE genes",
                                           conditionalPanel("input.ok3>0",
                                                            Module_Differential_Expression_UI("module")),
                                           conditionalPanel("input.ok3==0",
                                                            p("Please press start button in batch effect analysis tab"))
                                           ),
                                  tabPanel("Get list of Transcription factors",
                                           Enriched_markers_module_UI("module")
                                           
                                  ),
                                  tabPanel("Get list of marker genes",
                                           helpText("This tab enables you to search for any type of markers
                                                                 in the file that are enriched in the dataset uploaded
                                                    please make sure that the first two columns in the file
                                                    correspond to gene names present in mouse(column 1) and human(column 2"),
                                           fluidRow(
                                             column(8,
                                                    fileInput("add_file", label = h3('Choose file to upload'),
                                                              accept = c(
                                                                'text/csv',
                                                                'text/comma-separated-values',
                                                                'text/tab-separated-values',
                                                                'text/plain',
                                                                '.csv',
                                                                '.tsv')
                                                    ))),
                                           Enriched_markers_module_UI("module1")
                                           
                                  ),
                                  tabPanel("Transcription factor prediction",
                                           predicted_TF_module_UI("module")
                                  ),
                                  tabPanel("Visualization",
                                           tabsetPanel(
                                              tabPanel("MA plot",
                                                       MA_plot_UI("module")
                                                       ),
                                              tabPanel("Volcano plot",
                                                       Volcano_plot_UI("module")),
                                              tabPanel("Venn Diagram",
                                                       Venn_diagram_module_UI("module")),
                                              tabPanel("FC-FC plot",
                                                       FC_FC_plot_module_UI("module"))
                                           ))
                                
                                )
                       ),
                     tabPanel("Heatmap",
                              tabsetPanel(
                                tabPanel( "Heatmap ofDE genes",
                                          # actionButton("DE_go","GO"),
                                          # conditionalPanel("input.DE_go>0",
                                                           heatmap_module_UI("module")
                                                           #)
                                          
                                          ),
                                tabPanel("heatmap of TF",
                                         # actionButton("TF_go","GO"),
                                         # conditionalPanel("input.TF_go>0",
                                                          heatmap_module_UI("module1")
                                                         # )
                                         
                                         ),
                                tabPanel("heatmap of marker genes",
                                         heatmap_module_UI("module2")
                                         )
                                
                                )
                              ),
                     tabPanel("Enrichment Analysis",
                              tabsetPanel(
                                tabPanel("Enrichment analysis: KEGG",
                                         Kegg_module_UI("module")),
                                tabPanel("Enrichment analysis-GO:Biological processes",
                                         Biological_process_module_UI("module")),
                                tabPanel("Hallmark plot",
                                         Hallmark_module_UI("module"))
                              ))
                     ),
          tabPanel("ANOVA",
                   tabsetPanel(
                     tabPanel("ANOVA Table",
                              conditionalPanel("input.ok3>0",
                                               ANOVA_module_UI("module")    
                              )
                              ),
                     tabPanel("heatmap of 1000 most variable genes",
                              heatmap_module_UI("module3")
                     )
                   )
                   
                   ),
          tabPanel("Summary and results ", 
                   powerpoint_module_UI("module")
          )
          
                   
)
#################################Begin of code#################
# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 19MB.
options(shiny.maxRequestSize = 19*1024^2)

server<-function(input, output,session) {
  input_data<-callModule(Module_Raw_data_Input,"module")
  print("inside app- line 98 module 1 pass")
  
  #display output (count data and annotation file)
  # 
  #V contains two reactive values wherein the variable v$data represents the count/expression table
  # v$data is updated when either of the buttons "delete" or "reset" is clicked
   v <- reactiveValues(data=NULL,click1=NULL,raw_pca=NULL,norm_pca=NULL,heat_de=NULL,heat_tf=NULL,heat_anova=NULL,
                       dds=NULL,batch=NULL,anova=NULL,kegg=NULL,bp=NULL,hallmark=NULL,
                       de=NULL,TF=NULL,marker=NULL,data_wgcna=NULL,ma=NULL,vol=NULL,
                       pheno_wgcna=NULL,
                       wgcna_output=NULL,
                       wgcna_click=FALSE)
   #makeReactiveBinding("v")
#   
observeEvent(input_data$file2(),{
  
  if(!is.null(input_data$pData())) {
    #print('k')
    closeAlert(session,"exampleAlert_app")
    print("inside app line 111")
    #print(head(input_data$edata()))
    data<-input_data$edata()#edata()
    pheno<-input_data$pData()#pData()
    #set order of columns in expression data as same as order of sample  ID in pheno data
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Processing Data", value = 0)
    # Increment the progress bar, and update the detail text.
    progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    #print(as.vector(pheno[,1]))
    #print(pheno)
    data=data[,as.vector(pheno[,1])]
    print("line 291 app")
    print(head(data))
    v$data<-list(data,pheno,rlogTransformation(as.matrix(data))) #time consuming step
    # Increment the progress bar, and update the detail text.
    progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    #v$data<-list(data,pheno)
  }
  else {
    createAlert(session,"alert_app", "exampleAlert_app", title = "Oops!",
                content = paste0("Some samples in the expression/count table is absent in the annotation table!
                  Please give correct input"), append = FALSE)
  }

           
  })
observeEvent(input$del, {
  #v$data <- runif(100)
  #input$del
  print("inside app line 149")
  idx<-input$outliers_rows_selected
  data<-input_data$edata()
  pheno<-input_data$pData()
  #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]
  # print('data tab')
  #print(idx)
  #data<-input_data$edata()
  #   pheno<-input_data$pData()
  #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]
  # print('data tab')
 # print(idx)
  #set order of columns in expression data as same as order of sample  ID in pheno data
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Processing data", value = 0)
  # Increment the progress bar, and update the detail text.
  progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  # Pause for 0.1 seconds to simulate a long computation.
  #   Sys.sleep(0.1)
  if(!is.null(idx))
  {
    if(is.null(v$click1))
    {
      print("inside app line 181, delete button clicked hence deleteing selected samples")
      v$data<-list(data[,-idx],pheno[-idx,],rlogTransformation(as.matrix(data[,-idx])))
      print(ncol(v$data[[1]]))
      v$click1<-1
    }
    else
    {
      v$data[[1]]<-v$data[[1]][,-idx]
      v$data[[2]]<-v$data[[2]][-idx,]
      v$data[[3]]<-rlogTransformation(as.matrix(v$data[[1]][,-idx]))
    }
    v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1))
    normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
    callModule(Boxplot_module,"module1",normal)
    v$norm_pca<-callModule(PCA,"module1",normal)
    callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
    v_list<-callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                       v$dds$dds.fc()[[2]],
                       v$dds$normal())
    callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
    
  output$sov_data_input<-renderUI({
      print("create Checkbox")
      checklist = list()
      for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
        
        checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
      
      selectizeInput("sov_inp", "Choose the variables to be included
                     to compute the amount of contribution to variance",
                     checklist,selected = NULL,multiple=T)
      
    })
    observeEvent(input$sov_inp,
                 {
                   output$go_sov<-renderUI({
                     actionButton("go_sov","Go")
                   })
                   observeEvent(input$go_sov,
                                {
                                  output$sov<-renderPlotly({
                                    id<-as.numeric(input$sov_inp)
                                    print("inside app line 528")
                                    print(id)
                                    if(length(id)!=0){
                                    df<-source_of_variation_op(v$dds$normal(),
                                                               colData(v$dds$dds.fc()[[1]]),
                                                               id)
                                    if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                    print(df)
                                    p <- plot_ly(df, labels = ~Variation,
                                                 values = ~Percentage, type = 'pie',
                                                 textposition = 'inside',
                                                 textinfo = 'label+percent',
                                                 insidetextfont = list(color = '#FFFFFF')) %>%
                                      layout(title = 'Source of variation plot',
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                    }
                                    else plotly_empty()
                                  })
                                })
                 })
      
    #
    # Increment the progress bar, and update the detail text.
    progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
  }
})
#when reset button is clicked

observeEvent(input$reset, {
  #v$data <- rnorm(100)
  data<-input_data$edata()
  pheno<-input_data$pData()
  #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]
  #set order of columns in expression data as same as order of sample  ID in pheno data
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Processing data", value = 0)
  # Increment the progress bar, and update the detail text.
  progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  v$data<-list(data,pheno,rlogTransformation(as.matrix(data)))
  v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1))
  normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
  callModule(Boxplot_module,"module1",normal)
  v$norm_pca<-callModule(PCA,"module1",normal)
  callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
  callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                     v$dds$dds.fc()[[2]],
                     v$dds$normal())
  callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
  output$sov_data_input<-renderUI({
    print("create Checkbox")
    checklist = list()
    for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
      
      checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
    
    selectizeInput("sov_inp", "Choose the variables to be included
                   to compute the amount of contribution to variance",
                   checklist,selected = NULL,multiple=T)
    
  })
  observeEvent(input$sov_inp,
               {
                 output$go_sov<-renderUI({
                   actionButton("go_sov","Go")
                 })
                 observeEvent(input$go_sov,
                              {
                                output$sov<-renderPlotly({
                                  id<-as.numeric(input$sov_inp)
                                  print("inside app line 528")
                                  print(id)
                                  if(length(id)!=0){
                                    df<-source_of_variation_op(v$dds$normal(),
                                                               colData(v$dds$dds.fc()[[1]]),
                                                               id)
                                    if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                    print(df)
                                    p <- plot_ly(df, labels = ~Variation,
                                                 values = ~Percentage, type = 'pie',
                                                 textposition = 'inside',
                                                 textinfo = 'label+percent',
                                                 insidetextfont = list(color = '#FFFFFF')) %>%
                                      layout(title = 'Source of variation plot',
                                             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                    
                                  }
                                  else plotly_empty()
                                  
                                })
                              })
               })
  # Increment the progress bar, and update the detail text.
  progress$inc(2/2, detail = paste("Doing part", 2,"/",2))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
})  

  #Display expression table
  output$edata <- DT::renderDataTable({

    #print(head(v$data[[3]]))
    DT::datatable(v$data[[1]],class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 500,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list('copy',list(extend = 'collection',buttons = c('csv', 'excel'),text = 'Download only genes on clipboard'))#I('colvis')extend = 'collection',
                  )#list('copy',list(extend = 'collection',buttons = c('csv', 'excel'),text = 'Download'))
    )
  })
  
  output$downloadeData <- downloadHandler(
    filename = function() { 'exp.xlsx' },
    content = function(file) {
      #if(!is.null(data()[[1]]))
      #write.csv(v$data[[1]], file)
      write.xlsx2(v$data[[1]], file, sheetName = "Expression data",
                  col.names = TRUE, row.names = TRUE, append = FALSE)
    }
  )
  # #Display annotation table
  output$pData <- DT::renderDataTable({

    DT::datatable(v$data[[2]],class = 'cell-border stripe',
                  selection = list(target = 'column'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list('copy',list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download only genes on clipboard'))#I('colvis')
                  )#buttons = list('copy',list(extend = 'collection'
    )

  })
  # #outputOptions(output, 'pData', suspendWhenHidden=FALSE)
  #
  #dowmload annotation data
  output$downloadepData <- downloadHandler(
    filename = function() { 'annotation.xlsx' },
    content = function(file) {
      #write.csv(v$data[[2]], file)
      write.xlsx2(v$data[[2]], file, sheetName = "Annotation",
                  col.names = TRUE, row.names = TRUE, append = FALSE)
    }
  )

  #start button in raw data tab
  observeEvent(input_data$file2(),{
    #div(style="background-color:orange;",input$ok1)
    #toggleClass("ok1", "grey")
    output$button<-renderUI({
      actionButton("ok1", label = "Start")
    })

  })
 
  #start button turns blue on click
  #call boxplot and PCA for raw data input
  observeEvent(input$ok1,{
    
    toggleClass("ok1", "blue")
    callModule(Boxplot_module,"module",v$data)
    v$raw_pca<-callModule(PCA,"module",v$data)
 
    # creates a checkbox widget to select the important condition
    output$condition <-
      renderUI({
        #if (!is.null(input$file2)) {
        #print("create Checkbox condition")
        checklist = list()
        for (i in seq_along(colnames(v$data[[2]]))[-1]) {
          checklist[[colnames(v$data[[2]])[[i]]]] = i
        }
        radioButtons("conchoice", "Choose the conditon", checklist,selected = 2)
        #}
      })
    
    # creates a checkbox widget to select the variables to be included in design
    output$design <-
      renderUI({
        req(input$conchoice)
        checklist = list()
        for (i in seq_along(colnames(v$data[[2]]))[-1]) {
          # disply all variables except condition variable
          if(i!=as.numeric(input$conchoice)) checklist[[colnames(v$data[[2]])[[i]]]] = i
        }
        #print(checklist)
        if(length(checklist)!=0){
          checkboxGroupInput(session$ns("designchoice"), "Choose the variables to be included in design", checklist)
        }
      })
    output$button1<-renderUI({
      req(input$conchoice)
      actionButton("ok2", label = "Start")
    })
    #req(input$conchoice)
    
    
  })

observeEvent(input$ok2,
             {
               toggleClass("ok2", "blue")
               print("inside app line 324 calling normalized table module")
               
               v$dds<-callModule(Module_Normalized_data,"module",reactive(v$data),input$conchoice,input$designchoice,reactive(v$click1))
              # print("dat")
               normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
               #normalized table,annotation table,rlog transformed data for PCA
               callModule(Boxplot_module,"module1",normal)
               v$norm_pca<-callModule(PCA,"module1",normal)
               callModule(QC_normalized,"module",normal[[1]],zoom=TRUE)
               
               output$heatmap_sample<-renderPlot({
                 sampleDists <- dist(t(assay(v$dds$dds.fc()[[2]])))
                 var<-as.numeric(input$conchoice)
                 print("inside app line 368")
                 #print(colData(v$dds$dds.fc()[[1]])[,var])
                 heatmap_sample(sampleDists,v$dds$dds.fc()[[2]])
               })
               output$sov_data_input<-renderUI({
                     print("create Checkbox")
                     checklist = list()
                     for (i in seq_along(colnames(colData(v$dds$dds.fc()[[1]])))[-1]) { 
                       
                         checklist[[colnames(colData(v$dds$dds.fc()[[1]]))[[i]]]] = i}
                     
                     selectizeInput("sov_inp", "Choose the variables to be included
                                      to compute the amount of contribution to variance",
                                    checklist,selected = NULL,multiple=T)
                   
                 })
               observeEvent(input$sov_inp,
                            {
                              output$go_sov<-renderUI({
                                actionButton("go_sov","Go")
                              })
                              observeEvent(input$go_sov,
                                           {
                                             output$sov<-renderPlotly({
                                               id<-as.numeric(input$sov_inp)
                                               print("inside app line 528")
                                               print(id)
                                               if(length(id)!=0){
                                               df<-source_of_variation_op(v$dds$normal(),
                                                                          colData(v$dds$dds.fc()[[1]]),
                                                                                  id)
                                               if(sum(df$Percentage)<1) df<-rbind(df,c("Unknown",1-sum(df$Percentage)))
                                               print(df)
                                               p <- plot_ly(df, labels = ~Variation,
                                                            values = ~Percentage, type = 'pie',
                                                            textposition = 'inside',
                                                            textinfo = 'label+percent',
                                                            insidetextfont = list(color = '#FFFFFF')) %>%
                                                 layout(title = 'Source of variation plot',
                                                        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                                                        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
                                               }
                                               else plotly_empty()
                                             })
                                           })  
                            })
               
               
observeEvent(input$rem_samp,
             {
               v_list<-callModule(Remove_samples,"module",v$dds$dds.fc()[[1]],
                                                          v$dds$dds.fc()[[2]],
                                                          v$dds$normal())
                                                          
               print("inside app line 328 bar")
               # print(v_list$d())
               # print(v_list$bar())
               
               #display annotation table after removing outliers
               output$outliers <- DT::renderDataTable({
                 print("inside app line 334 outab")
                 dat <- as.data.frame(colData(v$dds$dds.fc()[[2]]))
                 # print("dat")
                 # print(head(dat))
                 # print(typeof(dat))
                 d <- v_list$d()
                 v1_clickBar_list<-v_list$bar()
                 idx<-NULL
                 if (!is.null(d) && is.null(v1_clickBar_list))#is.null(input$clickBar))
                 {
                   
                   idx<-which(dat[,1] %in% d)#d$key)
                 }
                 else if(is.null(d) && !is.null(v1_clickBar_list))
                 {
                   idx<-round(v1_clickBar_list)
                   #print(idx)
                 }
                 else if(!is.null(d) && !is.null(v1_clickBar_list))
                 {
                   id<-c(which(dat[,1] %in% d),round(v1_clickBar_list))
                   #print(id)
                   idx<-unique(id)
                 }
                 
                 library(DT)
                 DT::datatable(dat,class = 'cell-border stripe',
                               selection = list(mode='multiple',selected=idx,target='row'),
                               extensions = list('Scroller'=NULL,'Buttons'=NULL),
                               options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 350,scroller = TRUE,dom = 'Bfrtip',
                                              buttons = list()))#buttons = list('copy',list(extend = 'collection',
               })
               
               print("inside app line 372")
            
               
             })
observeEvent(input$ok3,
             {
               print("inside app line 406")
               print(head(v$dds$dds.fc()[[1]]))
               v$batch<-callModule(Module_Batch_effect,"module",input$conchoice,input$designchoice,reactive({v$dds$dds.fc()}),reactive({v$dds$normal()}))
               
               normal<-list(v$dds$normal(),colData(v$dds$dds.fc()[[1]]),assay(v$dds$dds.fc()[[2]]))
               
               observeEvent(input$go_wgcna,
                            {
                              #Estimate size factors  
                              dds.norm=estimateSizeFactors(v$dds$dds.fc()[[1]])
                              if(v$batch$batch_choice()>1)
                              {
                                #callModule(WGCNA_module,"module",infile=reactive({v$data_wgcna}),de=TRUE,condition)
                                condition<-colData(v$batch$batch_data_for_DESeq()[[1]])[,as.numeric(input$conchoice)]
                                v$wgcna_output<-callModule(WGCNA_module,"module",
                                                           infile=reactive({v$batch$batch_corrected_data()}),
                                                           perform_voom=FALSE,condition,
                                                           reactive({v$dds$dds.fc()[[1]]}),
                                                           reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                                                           reactive({input$conchoice}),
                                                           reactive({input_data$organism()}),
                                                           reactive({v$de$ok3()}),
                                                           reactive({v$de$combination}),
                                                           reactive({v$anova$anova_table()}),
                                                           reactive({v$batch$batch_choice()}),
                                                           reactive({v$batch$batch_corrected_data()}),
                                                           reactive({v$dds$normal()})
                                                           )
                                
                              }
                              else
                              {
                                #callModule(WGCNA_module,"module",infile=reactive({v$data_wgcna}),de=TRUE,condition)
                                condition<-colData(v$batch$batch_data_for_DESeq()[[1]])[,as.numeric(input$conchoice)]
                                
                                v$wgcna_output<-callModule(WGCNA_module,"module",
                                                           infile=reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                                                           perform_voom=TRUE,condition,
                                                           reactive({v$dds$dds.fc()[[1]]}),
                                                           reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                                                           reactive({input$conchoice}),
                                                           reactive({input_data$organism()}),
                                                           reactive({v$de$ok3()}),
                                                           reactive({v$de$combination}),
                                                           reactive({v$anova$anova_table()}),
                                                           reactive({v$batch$batch_choice()}),
                                                           NULL,
                                                           reactive({v$dds$normal()})
                                                           )
                                
                              }
                              
                              print("inside app line 484")
                            })
               
               v$de<-callModule(Module_Differential_Expression,"module",input$conchoice,
                                reactive({v$batch$batch_data_for_DESeq()}),
                                #reactive({v$wgcna_click}),
                                reactive({v$wgcna_output}),
                                reactive({v$dds$normal()}),
                                reactive({v$batch$batch_choice()}),
                                reactive({v$batch$batch_corrected_data()}),
                                reactive({v$anova$anova_table()})
                                )
               
               callModule(FC_FC_plot_module,"module",reactive({v$de$de_genes}),
                          reactive({v$de$combination})
                          )
               
               # 
               
               v$ma<-callModule(MA_plot,"module",reactive({v$de$combination}),
                          reactive({v$de$de_genes()}),
                          reactive({v$de$p_value()}))
               v$vol<-callModule(Volcano_plot,"module",
                          reactive({v$de$combination()}),
                          reactive({v$de$de_genes()}),
                          reactive({v$de$p_value()}),
                          reactive({v$de$hypothesis_choice()}))
               
               callModule(Venn_diagram_module,"module",reactive({v$de$de_genes}),
                                     reactive({v$de$combination}),
                                     reactive({v$wgcna_output}))
               
                             v$heat_de<-callModule(heatmap_module,"module","DE",NULL,
                                        reactive({v$batch}),
                                        reactive({v$de$de_genes}),
                                        reactive({v$de$combination}),
                                        reactive({v$wgcna_output})
                                        # reactive({v$batch$batch_choice}),
                                        # reactive({v$batch$batch_corrected_data})

                             )
                             #}
                           #}) 
            
               
               
               observeEvent(input$add_file,
                            {
                              filepath<-input$add_file$datapath
                              TF_list<-read.csv(filepath, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                              
                              v$marker<-callModule(Enriched_markers_module,"module1",reactive({TF_list}),
                                                 reactive({input_data$organism()}),
                                                 reactive({v$de$de_genes}),
                                                 reactive({v$de$combination}),
                                                 reactive({v$wgcna_output}),
                                                 input$conchoice,
                                                 reactive({v$dds$normal()}),
                                                 reactive({v$batch$batch_data_for_DESeq()}),
                                                 reactive({v$batch$batch_choice()}),
                                                 reactive({v$batch$batch_corrected_data()}),
                                                 reactive({v$anova$anova_table()}),
                                                 reactive({v$dds$dds.fc()[[1]]})
                                                 )
                              
                              callModule(heatmap_module,"module2","TF",NULL,
                                         reactive({v$batch}),
                                         reactive({v$marker$de_genes}),
                                         reactive({v$de$combination}),
                                         reactive({v$wgcna_output})
                                         # reactive({v$batch$batch_choice}),
                                         # reactive({v$batch$batch_corrected_data})
                                         
                              )
                              
                              
                            })
               
               TF_list<-read.csv("./www/Transcriptome_TFcat.txt", header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
               
          v$TF<-callModule(Enriched_markers_module,"module",
                          reactive({TF_list}),
                          reactive({input_data$organism()}),
                          reactive({v$de$de_genes}),
                          reactive({v$de$combination}),
                          reactive({v$wgcna_output}),
                          input$conchoice,
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_data_for_DESeq()}),
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_corrected_data()}),
                          reactive({v$anova$anova_table()}),
                          reactive({v$dds$dds.fc()[[1]]})
                          )
          
          v$heat_tf<-callModule(heatmap_module,"module1","TF",NULL,
                     reactive({v$batch}),
                     reactive({v$TF$de_genes}),
                     reactive({v$de$combination}),
                     reactive({v$wgcna_output})
                     # reactive({v$batch$batch_choice}),
                     # reactive({v$batch$batch_corrected_data})
                     
          )
          
           v$anova<-callModule(ANOVA_module,"module",
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$batch$batch_data_for_DESeq()[[2]]}),
                          reactive({v$de$combination}),
                          reactive({input$conchoice}),
                          reactive({v$de$de_genes}),
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_corrected_data()})
                          
                         
               )
           v$heat_anova<-callModule(heatmap_module,"module3","ANOVA",
                                 reactive({v$anova$dds()}),
                                 reactive({v$batch}),
                                 reactive({v$de$de_genes}),
                                 reactive({v$de$combination}),
                                 reactive({v$wgcna_output})
                                 # reactive({v$batch$batch_choice}),
                                 # reactive({v$batch$batch_corrected_data})
                                 
           )
               
               # DE_genes
               # combination
               # anova_table
               # dds.fc()[[1]]#batch design
               # conchoice
               # wgcna_output
               # organism
           
               callModule(predicted_TF_module,"module",
                          reactive({v$de$de_genes}),
                          reactive({v$batch$batch_data_for_DESeq()}),
                          reactive({v$anova$anova_table()}),
                          reactive({v$de$combination}),
                          reactive({input$conchoice}),
                          reactive({v$wgcna_output}),
                          reactive({input_data$organism()}),
                          reactive({v$dds$normal()}),
                          reactive({v$batch$batch_choice()}),
                          reactive({v$batch$batch_corrected_data()})
                          
               )
               v$kegg<-callModule(Kegg_module,"module",reactive({v$de$de_genes}),
                                               reactive({input_data$organism()}),
                                               reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                                               reactive({v$de$combination}),
                                               reactive({v$wgcna_output}),
                                               reactive({v$anova$anova_table()})
                          )
                                                                #dds.fc()
               v$bp<-callModule(Biological_process_module,"module",reactive({v$de$de_genes}),
                          reactive({input_data$organism()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$de$combination}),
                          reactive({v$wgcna_output}))
               
               v$hallmark<-callModule(Hallmark_module,"module",reactive({v$de$de_genes}),
                          reactive({input_data$organism()}),
                          reactive({v$batch$batch_data_for_DESeq()[[1]]}),
                          reactive({v$de$combination}),
                          reactive({v$wgcna_output}))
               
               
               callModule(powerpoint_module,"module",
                          reactive({v$dds$normal()}),reactive({v$dds$dds.fc()[[1]]}),reactive({input$conchoice}),
                          reactive({v$raw_pca()}),reactive({v$norm_pca()}),reactive({normal}),reactive(v$data),
                          reactive({v$batch$top_batch_pca()}),reactive({v$batch$batch_corrected_data()}),reactive({v$batch$batch_choice()}),
                          reactive({v$de$combination}),reactive({v$de$ok3()}),reactive({v$de$de_genes()}),reactive({v$TF$de_genes()}),
                          reactive({v$ma$input_scale()}),reactive({v$ma$p_values()}),
                          reactive({v$vol$input_scale_volx()}),reactive({v$vol$input_scale_voly()}),reactive({v$vol$hypothesis_choice()}),
                          reactive({v$heat_de$input_Distance()}),reactive({v$heat_de$input_Linkage()}),
                          reactive({v$heat_tf$input_Distance()}),reactive({v$heat_tf$input_Linkage()}),
                          reactive({v$anova$dds()}),reactive({v$heat_anova$input_Distance()}),reactive({v$heat_anova$input_Linkage()}),
                          reactive({v$kegg$Enriched_Kegg_table()}),reactive({v$kegg$Enriched_Kegg_obj()}),
                          reactive({v$bp$Enriched_BP_table()}),reactive({v$bp$Enriched_BP_obj()}),
                          reactive({v$hallmark$Enriched_hall_table()}),reactive({v$hallmark$Enriched_hall_obj()})
                          )
               print("inside app line 455")
              
               
               
             })
             })


    
output$file_input1<-renderUI({
        req(input$option_wgcna)
  if(as.numeric(input$option_wgcna)==1)
  {
        fileInput("input_file",label=h5("Please upload normalized/batch corrected data"),
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'))
  }
})
        
  observeEvent(input$input_file,
                     {
                       filepath<-input$input_file$datapath
                       if(!is.null(filepath))
                       {
                         input_exp<-read.csv(filepath, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                         data<-input_exp[,-1]
                         
                         #convert entries to integers
                         data=as.matrix(data)
                         storage.mode(data)="double"
                         data= data.frame(data)
                         rownames(data)=input_exp[,1]
                         v$data_wgcna<-data
                         
                         output$file_input2<-renderUI({
                           fileInput("input_file2",label=h5("Please upload annotation data"),
                                     accept = c('text/csv',
                                                'text/comma-separated-values',
                                                'text/tab-separated-values',
                                                'text/plain',
                                                '.csv',
                                                '.tsv'))
                         })
                         
                         observeEvent(input$input_file2,
                                      {
                                        filepath2<-input$input_file2$datapath
                                        if(!is.null(filepath2))
                                        {
                                          input_p<-read.csv(filepath2, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                                          pheno<-input_p#[,-1]
                                          print("line")
                                          print(pheno)
                                          v$pheno_wgcna<-get_pheno(v$data_wgcna,pheno)
                                          
                                        }
                                        
                                      })
                       }
                       
                     })
      
  
output$pheno_table <- DT::renderDataTable({
  print("line2")
  print(v$pheno_wgcna)
  DT::datatable(v$pheno_wgcna,class = 'cell-border stripe',
                selection = list(mode='single',target = 'column'),
                extensions = list('Scroller'=NULL,'Buttons'=NULL),
                options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 200,scroller = TRUE,dom = 'Bfrtip',
                               buttons = list('copy',list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download only genes on clipboard'))#I('colvis')
                )#buttons = list('copy',list(extend = 'collection'
  )
  
})
output$go_wgcna<-renderUI({
  req(input$option_wgcna)
  if(as.numeric(input$option_wgcna)==1)
  {
    req(input$pheno_table_columns_selected)
    print(as.numeric(input$pheno_table_columns_selected))
  if (as.numeric(input$pheno_table_columns_selected)>1) {actionButton("go_wgcna","GO")} #any variable other than sample id(assuming saple id to be the first column)
  }
  else if(as.numeric(input$option_wgcna)==2) actionButton("go_wgcna","GO")
  
})
observeEvent(input$pheno_table_columns_selected,{
  print('hey')
  print(input$pheno_table_columns_selected)
  col<-as.numeric(input$pheno_table_columns_selected)
  condition<-v$pheno_wgcna[,col]
  print("condition")
  #didsplay help text
  
  observeEvent(input$go_wgcna,
               {
                 callModule(WGCNA_module,"module",
                            infile=reactive({v$data_wgcna}),perform_voom=FALSE,condition,
                            NULL,
                            NULL,
                            NULL,
                            NULL,
                            reactive({0}),
                            NULL,
                            NULL
                            )
                 
               })
  
})




}
shinyApp(ui = ui, server = server)

