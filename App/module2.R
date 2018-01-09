source("./functions.r")
source("./gene_count_module.r")
#source("E:/Zenitha/App/module1.r")
Module_Normalized_data_UI<-function(id)
{
  ns<-NS(id)
  tagList(
          fluidRow(column(10,bsAlert("message"))),
          fluidRow(column(5,actionButton(ns("help"), "Help"))),
          fluidRow(column(8,
                            bsAlert("alert"),
                            plotlyOutput(ns("cut_off_plot"))
                          ),
                   column(4,textInput(ns("cutoff"), "cut-off", value = NULL, width = NULL, placeholder = NULL)),
                   
                   column(12,
                             # normalized data table
                             downloadButton(ns('download'), 'Download full Data'),
                             DT::dataTableOutput(ns("normal")),
                             gene_count_module_UI(ns("module"))
                   #           uiOutput("plot_option_gene_n"),
                   #           downloadButton('download_norm_genecount', 'Download plot'),
                   #           plotOutput("gene_plot_n"))
                  )),
          
          bsModal(ns("modalhelp"), "Help for error model matrix not full rank", ns("help"),size = "large",
                             helpText("There are two main reasons for this problem: either one or more
                                      columns in the model matrix are linear combinations of other columns, or there
                                      are levels of factors or combinations of levels of multiple factors which are missing
                                      samples. We address these two problems below"),
                             helpText("Case 1:"),
                             fluidRow(column(5,img(src="error1.png"))),
                             helpText(" In the above table two variables contain exactly the same information. The software cannot
                                      fit an effect for batch and condition, because
                                      they produce identical columns in the model matrix. This is also referred to as
                                      perfect confounding."),
                             p("Solution:"),
                             helpText("the batch effect cannot be fit and must be removed
                                      from the model formula. There is just no way to tell apart the condition effects
                                      and the batch effects. The options are either to assume there is no batch effect
                                      (which we know is highly unlikely given the literature on batch effects in sequencing
                                      datasets) or to repeat the experiment and properly balance the conditions across
                                      batches."),
                             helpText("Case 2:"),
                             fluidRow(column(5,img(src="error2.png"))),
                             helpText("In the above table the variables are not identical,
                                      but one variable can be formed by the combination of other factor levels. In the
                                      following example, the effect of batch 2 vs 1 cannot be fit because it is identical
                                      to a column in the model matrix which represents the condition C vs A effect."),
                             p("Solution:"),
                             helpText("the batch effect cannot be fit and must be removed
                                      from the model formula. There is just no way to tell apart the condition effects
                                      and the batch effects. The options are either to assume there is no batch effect
                                      (which we know is highly unlikely given the literature on batch effects in sequencing
                                      datasets) or to repeat the experiment and properly balance the conditions across
                                      batches."),
                             helpText("Case 3:"),
                             fluidRow(column(5,img(src="error3.png"))),
                             helpText("In the above case where we can in fact perform inference.In the above example the experiment
                                      comprises of grouped individuals, where the goal is to test the group-specific effect of
                                      a treatment, while controlling for individual effects."),
                             p("Solution:"),
                             helpText("unfortunately the tool cannot be currently used for a similar scenario
                                      as above.We are working on it.You will have to use DESeq2 to define the design."),
                             helpText("Case 4:"),
                             fluidRow(column(5,img(src="error4.png"))),
                             helpText("In the above case as highlighted a level is missing from a factor.The group 3 does not have level C"),
                             p("Solution:"),
                             helpText("Remove these samples with missing level.In the above case remove samples 13 to 16")
                           )#,
               
    )
  
}
#normalized tab
#create a gene plots module
Module_Normalized_data<-function(input,output,session,vdata,conchoice,designchoice,v_click1)
{
  data<-vdata()
  print(ncol(data[[1]]))
  trim <- reactive({
    #req(input$conchoice)
    tryCatch(
      
      #check if user has provided the choice of treatment/condition variable
      # if (!is.null(input$conchoice))#((input$ok1 > 0 )&& (!is.null(input$conchoice)))#
       {
         print("inside module 2 line 86")
        print("treatment/condition variable provided")
        #Get variable to consider as treatment/condition
        condition<-as.numeric(conchoice)
        #print(condition)
        
        #===============Main purpose of this reactive==========================#
        #Step 1: prepare for DESeq2
        #Inputs needed: expression/count table, annotation table and design formula
        
        #Step 2: prepare expression/count table for trimming (Quality control used to remove genes with low counts)
        #=======================================================================#
        
        #Input 1: Get expression data
        edata <-data[[1]]
        #print(head(edata))
        
        #Input 2: Get annotation data
        pheno<-data[[2]]
        #print(pheno)
        
        #Input 3: Get design formula
        design<-as.numeric(designchoice)
        #print(design)
        
        #Rename treatment/condition variable to condition
        #This is done to ensure one standard name to refer to , in the future computation steps
        colnames(pheno)[condition]<-"condition"
        
        #Since the design formula contains the name of the treatment/condition variable(based on input annotation table)
        #We rename the variable to "condirion" in the design formula
        #Updating design formula
        d<-''
        for (i in design)
        {
          #print(i)
          d<-paste(d,colnames(pheno)[i]," + ",sep=" ")
        }
        #print(d)
        d<-paste('~',' ',d,colnames(pheno)[condition],sep = " ")
        
        #print new formula
        #print(formula(d))
        
        #once we have obtained the input we store it in a container (a subclass)
        #DESeqDataSet is a subclass of RangedSummarizedExperiment, used to store the input values,
        #intermediate calculations and results of an analysis of differential expression.
        #check dese2 vignette for more information.
        
        #Prepare DESeq dataset using inputs
        dds<-DESeqDataSetFromMatrix(edata,colData=pheno,formula(d))
        
        #A check to see if the annotation table has been correctly stored and if the value can be retrieved
        print("inside module 2 line 139")
        print('Annotation table')
        print(colData(dds))
        
        #Normalization: DESeq2 performs normalization by first computing size factor for each sample
        
        #Estimate size factors  
        dds.norm=estimateSizeFactors(dds)
        
        #call function min_samples_three to get a list of treatment/condition groups that contain only one or two samples
        samp<-min_samples_three(dds.norm)
        
        #prior to proceeding with normalization, if treatment/condition groups have one or two samples
        #then user is alerted through a warning message (error message if a treatment/condition group has only one sample) 
        
        if(is.null(samp)) {
          #print('k')
          closeAlert(session,"exampleAlert")
        }
        else if(samp[[2]]==1) { #get treatment/condition groups with only one sample
          if(length(samp[[1]])>1) #if there are more than one treatment/condition group with only one sample
          {
            x<-strsplit(samp[[1]],' ') #Get names of all treatment/condition groups with only one sample
            # print('k')
            # print(x)
            # print(length(x))
            if(!is.null(x)) closeAlert(session,"exampleAlert")
            sam<-NULL
            for (i in 1:(length(x)-1))
            {
              sam<-paste0(sam,x[[i]],',')
            }
            sam<-paste0(substr(sam,1,nchar(sam)-1),' and ',x[[length(x)]])
            #print(sam)
            #Alert user with the  names of all treatment/condition groups with only one sample
            createAlert(session,"alert", "exampleAlert", title = "Oops!",
                        content = paste0("The conditions ",sam," have only 1 sample.Either add additional samples for these conditions or remove the condition."), append = T)
          }
          else{
            #Alert user with the  name of  treatment/condition group with only one sample
            createAlert(session,"alert", "exampleAlert", title = "Oops!",
                        content = paste0("The condition ",samp[[1]]," has only 1 sample.Either add additional samples for this condition or remove the condition."), append = T)
          }
          
          
        }
        else if(samp[[2]]==2) #get treatment/condition groups with only two samples
        {
          if(length(samp[[1]])>1)#if there are more than one treatment/condition groups with only two samples
          {
            x<-strsplit(samp[[1]],' ')#Get names of all treatment/condition groups with only two samples
            # print('k')
            # print(x)
            # print(length(x))
            print("module 2 line 192")
            if(!is.null(v_click1())) print(x)
            if(!is.null(x)) closeAlert(session,"exampleAlert")
            sam<-NULL
            for (i in 1:(length(x)-1))
            {
              sam<-paste0(sam,x[[i]],',')
            }
            print(sam)
            sam<-paste0(substr(sam,1,nchar(sam)-1),' and ',x[[length(x)]])
            print(sam)
            
            #Alert users with names of all treatment/condition groups with only two samples
            createAlert(session,"alert", "exampleAlert", title = "Warning!",
                        content = paste0("The conditions ",sam," have only 2 samples"), append = T)
            
          }
          else
          {
            #Alert user with name of treatment/condition group with only two samples
            createAlert(session,"alert", "exampleAlert", title = "Warning!",
                        content = paste0("The condition ",samp[[1]]," has 2 samples"), append = T)
            
          }
          
        }
        
  
        #Step 2: Trimming: Remove genes with low counts
        #part 1: Compute for each gene the maximum of average mean in each treatment/condition group 
        #part 2: Filter genes with low counts (max average value across all conditions<threshold(user input) )
        #part two is computed in the reactive dds.fc
        
        #This reactive only computes part 1
        
        #Apply row and column statistics
        dt1<-as.data.frame(sapply( levels(dds.norm$condition), function(lvl)
        {
          #Get Average eaxpression of a gene across each treatment/condition group
          rowMeans( counts(dds.norm,normalized=TRUE)[,dds.norm$condition == lvl] )
          
        }))
        
        #The dataframe dt1 contains average expression of a gene in each treatment/condition group
        #row->gene and column->all treatment/condition groups
        #print(head(dt1))
        #Add new column labelled as "max" to dt1
        #Pick maximum of average mean among all conditions for each gene and assign it to "max"
        dt1[,"max"]<-as.numeric(apply(dt1,1,max))#min))#max))
        #print(head(dt1))
        
        #return deseq2 dataset(contains all input for statistical testing(hypothesis testing and ANOVA))
        #       and dataframe conintaing the average expression of a gene in each treatment/condition group + maximum of this average value
        list(dds,dt1)
        
      }
      #It is possible that two variable in the annotation table can coontain the same information
      #meaning the variables explain the same variance in the dataset in which case DESEq2 throws an error
      #"the model matrix is not full rank". In such a case the user may have to remove one of these
      #variables from the annotation table or combine it as one. For more information rgearding 
      #the problem and what to do please check for the error in DESeq2 vignette.
      , error = function(c) {
        print(typeof(c$message))
        print(substr(c$message,1,33))
        if("the model matrix is not full rank" == substr(c$message,1,33))  c$message <- paste0("Error:The model matrix is not full rank. To know more about this please click the help button on this page.")
        else c$message
        stop(c)
      })
    
  })
#Prior to trimming( remove genes with low counts), we need to get the threshold from the user
#The user decides the the threshold based on the cut off plot.
# The cut off plot is an interactive plot which indicates the number of genes remaining in the dataset after
# trimming. The user can hover over the plot to see the threshold value and the corresponding genes that would remain in the 
# dataset if the value is chosen
  
#This reactive computes the input for the cut off plot  
  cutof<-reactive({
    #if(!is.null(input$conchoice)){
      #call reactive trim
      # get dataframe conintaing the (average expression of a gene 
      # in each treatment/condition group + maximum of this average value)
      matrix<-trim()[[2]]
      print("Inside module 2 line 269 cutof reactive")
      dat<-as.data.frame(matrix(NA, nrow = 200, ncol = 2))
      colnames(dat)<-c("x","y")
      
      #Loop to compute the number of genes that would remain in the dataset
      # for thresholds one to 200.
      for(i in 0:200)
      {
        #print(i)
        dat[i,1]<-i #Threshold value (value for x-axis)
        dat[i,2]<-length(which(matrix$max>i)) #Get number of remaining genes after trimming (value for y-axis)
      }
      #print(head(x))
      #print(head(y))
      #dat<-as.data.frame(x,y)
      #print(head(dat))
      dat
    #}
  }) 
# Display cut off plot
  output$cut_off_plot <-renderPlotly({
    
    #req(conchoice)
    # if(!is.null(conchoice))
    #  {
      dat<-cutof()
      #print(head(dat))
      # Plotly plot
      plot_ly(x=dat$x,y=dat$y,type = "scattergl",mode = 'markers',
              hoverinfo = 'text',
              text = ~paste('Cut-off:', dat$x, 
                            '</br> Total gene counts:', dat$y)) %>%
        #add_trace(x=dat$x[1],y =dat$y[1], mode = "lines")%>%
        #add_markers(type = "scattergl")%>%
        layout(
          xaxis = list(title = "Cut-off", gridcolor = "#bfbfbf", domain = c(0, 0.98)),
          yaxis = list(title = "Total gene counts", gridcolor = "#bfbfbf"))
     # }
     # else plotly_empty()
    
  })
  
#Trim expression/count data using threshold if provided by user and update DESeq2 dataset 
dds.fc<-reactive({
    
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Preparing normalized table", value = 0)
    
    # Number of times we'll go through the loop to update progress bar
    n <- 3
    #Step 1: Retieve DESeq2 dataset object 
    dds<-trim()[[1]]#isolate(trim()[[1]])
    #if(!is.null(v_click1)) dds<-trim()[[1]]#incase samples have been selected to be deleted
    print('module 2 line 325 dds')
    #print(dds)
    
    #In order to visualize PCA we need to transform count data. Hence we perform rlogtransformation.
    #define variable to contain transformed data
    rld<-NULL
    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    #Step 2: check if user has provided threshold to trim data, if not then skip trimming and proceed to normalization
    #print("cutoff")
    #print(input$cutoff)
    if(input$cutoff=="")
    {
      print("pass")
      dds.fc<-dds
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      rld<-rlogTransformation(dds.fc)
    }
    else{
      #Prior to trimming we need to get threshold from the user 
      #req(input$cutoff)
      
      #Step 3: Extract count table from dds to prepare for trimming
      data<-assay(dds)
      #set order of columns in expression data as same as order of sample  ID in annotation data
      #Extract annotation table
      pheno<-colData(dds)
      data<-data[,as.vector(pheno[,1])]
      #define variable to contain filtered data
      filtered<-NULL
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      #Obtain indexes of genes remaining after cut off
      filtered_genes<-which(trim()[[2]]$max>=as.numeric(input$cutoff))
      #filter count data of low expressed genes using indices in the pervious step
      filtered<-data[filtered_genes,]
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      #create a new DESeq2 dataset with filtered genes
      #Get inputs
      #Input 1: expression/count table (containined in filtered variable computed in previous step)
      #Input 2: annotation table provided by user (already extracted in line 387)
      #Input 3: design formula constructed using treatment variable specified by user (check reactive labelled trim)
      design<-design(dds)
      #create new DESEq2 dataset 
      dds.fc<-DESeqDataSetFromMatrix(filtered,colData=pheno,design)
      # print('dds.fc')
      # print(colData(dds.fc))
      # print(dds.fc)
      #Also perform rlog transformation so that the transformed data can be used later
      rld<-rlogTransformation(dds.fc)
     }

    # Increment the progress bar, and update the detail text.
    progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
    
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    list(dds.fc,rld)
  })

#Once dataset has been trimmed by removing genes with low counts we proceed to the next step of pre-processing
#Perform normalization
#Get normalized data
normal <-reactive({
   print("inside module 2 line 405 normal")

   # if(!is.null(input$conchoice))
   # {
    dds<-dds.fc()[[1]]
    
    #Estimate size factors
    dds.norm=estimateSizeFactors(dds)
    
    #Get normalized data
    norm<-counts(dds.norm, normalized=TRUE)
    #print("normalized table")
    #print(head(norm))
    if(!is.null(norm))
    {
      createAlert(session,"message","exampleAlert1", title="Message: To proceed to next step",
                  content = "Click on on outlier and batch effect check button.", append=FALSE)
      
    }
    else closeAlert(session,"exampleAlert1")
    norm
 # }
})

# Display normalized table
output$normal <-
  DT::renderDataTable({
    #print("inside normalized table")
    #print(input$conchoice)
    #req(input$conchoice)
     # if(!is.null(input$conchoice))# & input$cutoff !="")
     # {
      DT::datatable(normal(),class = 'cell-border stripe',
                    selection = list(mode='single',target = 'row'),
                    extensions = list('Scroller'=NULL,'Buttons'=NULL),
                    options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 250,scroller = TRUE,dom = 'Bfrtip',
                                   buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download only genes on clipboard'))#I('colvis')
                    )
      )
    #}
  })
#display gene expression across conditon
observeEvent(input$normal_rows_selected,{
  #get which gene was clicked
  #get the row  number of the gene clicked
  print('hey')
  print(input$normal_rows_selected)
  selected_row <- input$normal_rows_selected
  print(selected_row)
  
  #get the normalized data
  full_data<-normal()
  print(head(full_data))
  print('hey')
  
  #get the gene
  an_gene<-rownames(full_data)[selected_row]
  print(an_gene)
  counts<-as.vector(full_data[selected_row,])#[which(an_gene %in% rownames(full_data)),])
  print(counts)
  print(as.factor(colData(dds.fc()[[1]])[,as.numeric(conchoice)]))
  cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
  #rownames(cond)<-colData(dds.fc()[[1]])[,1]
  library('data.table')
  print(length(counts))
  print(length(cond))
  df<-data.frame(counts,cond)
  colnames(df)<-c('count','condition')
  print(head(df))
  callModule(gene_count_module,"module",NULL,reactive({dds.fc()[[1]]}),reactive({an_gene}))
})

#download normalized data
output$download <- downloadHandler(
  filename = function() { 'normalized_table.xlsx' },
  content = function(file) {
    #req(input$conchoice)
    # if(!is.null(input$conchoice))# & input$cutoff !="")
    # {
    #write.csv(normal(), file)
    write.xlsx2(normal(), file, sheetName = "Normalized table",
                col.names = TRUE, row.names = TRUE, append = FALSE)
   # }
  }
)

  dds<-list(dds.fc=reactive({dds.fc()}),normal=reactive({normal()}))
  return(dds)


}










