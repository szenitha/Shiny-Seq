#install.packages("shiny")
# library(shiny)
# #library(shinysky)
# #install.packages("shinyjs")
# library(shinyjs)
# #install.packages("shinyBS")
# library(shinyBS)

#raw data tab
#tabPanel("Raw Data",

Module_Raw_data_UI<- function(id)
{
  ns<-NS(id)
  tagList(
    
    #Select organism
    fluidRow(column(3,
                    selectInput(ns("organism"), label = h5("Select organism"), 
                                choices = list("","Human" = 1, "Mouse" = 2),
                                selected = NULL))),
    uiOutput(ns("inp_ch")),
    
    fluidRow(column(4,uiOutput(ns("file3")))),


    #Read expression table
    fluidRow(column(4,uiOutput(ns("file1")))),

    #Read annotation table
    # #conditionalPanel("input.file1",
    fluidRow(column(4,uiOutput(ns("file2"))))#,
    #
  )
}

Module_Raw_data_Input<-function(input,output,session)
{
  
  #Step 1: Raw data tab
  #Aim: Display count/expression table and annotation table
  #Procedure:
  #Option 1: Generate count table from kallisto file and upload annotation table
  #Take kallisto file path as input
  #Upload annotation table
  #Option 2: Upload Count/expression table and annotation table
  #Upload count table/expression table and annotation table from file
  #observeEvent(input$reset_button, {js$reset()})  
  
  #Function that reads data from a file
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  print("module 1 line 56 ok")
  #print("hey")
  input_file <- function(in_file)
  {
    if(is.null(in_file))
    {
      return (NULL)
    }
    else
      #print(in_file$datapath)
      return(read.csv(in_file$datapath, header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")) 
  }
  
  
   output$inp_ch<-renderUI({
    # req(input$organism)
    #validate(need(try(input$organism!=""),"Please choose organism under study"))
    if(input$organism!="")
    {
      #print(input$organism)
      radioButtons(session$ns("filechoice"), "Choose starting point", choices = list("kallisto"=1, "count table"=2),selected = 2 )
    }
      })
   
  #Option 1:  start with kallisto
  #Input kallisto file path to compute count table/expression table to be displayed
               # {
   output$file3<-renderUI({
                    
                    # validate(need(try(input$organism!=""),"Please choose organism under study"))
                    #validate(need(try(as.numeric(input$filechoice)!=0),"Please choose input source"))
                    # print(as.numeric(input$filechoice))
                    req(input$filechoice)
                    if(as.numeric(input$filechoice)==1 && input$organism!="")
                    {
                      #shinyDirButton("dir", "Chose directory","Upload")
                      textInput(session$ns("kallisto_path"), "Kallisto path", value = "")#E:/Zenitha/kallisto")
                      
                    }
                  })



  #compute expression table/count table
  #t2g--->txi
  #This reactive prepeares input for the reactive txi which generates count table from kallisto files
  t2g<-reactive({
    #validate(need(try(input$organism!="","Please choose organism under study")))
    if (!exists("mart")) {
      org<-NULL
      if(as.numeric(input$organism==1)) org<-"hsapiens_gene_ensembl"
      else if(as.numeric(input$organism==2)) org<-"mmusculus_gene_ensembl"
      mart <-
        biomaRt::useMart(host = "jul2015.archive.ensembl.org",#"uswest.ensembl.org",#"grch37.ensembl.org",#,"ensembl.org",#
                         biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = org)
      print("mart")
      print(mart)
      t2g <- biomaRt::getBM(
        attributes = c("ensembl_transcript_id",#"ensembl_gene_id",#
                       "external_gene_name"),
        mart = mart
      )
      t2g
    }
  })
  # create the count matrix out of the kallisto abundance table
  txi <-
    reactive( {
      #validate(need(try(input$kallisto_path!=""),"Please enter path of kallisto file"))
      if (as.numeric(input$filechoice)==1 && input$kallisto_path!=""){
        print("module 1 line 126 ok")
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())

        progress$set(message = "Processing kallisto files", value = 0)

        print(dir(input[["kallisto_path"]]))

        #print(list.files(input[["kallisto_path"]], pattern="*.tsv"))#, full.names=TRUE)
        print(sapply(dir(input[["kallisto_path"]]),
                     function(id)
                       #file.path(input[["kallisto_path"]], id,"*.tsv")
                       #grep("tsv", file.path(input[["kallisto_path"]], id))
                       list.files(file.path(input[["kallisto_path"]],id), pattern="*.tsv", full.names=TRUE)
        ))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/2, detail = paste("Doing part", 1,"/",2))
        print(t2g())
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        d<-as.data.frame(tximport(sapply(dir(input[["kallisto_path"]]),
                                         function(id)
                                         {
                                           #file.path(input[["kallisto_path"]], id,"*.tsv")
                                           list.files(file.path(input[["kallisto_path"]],id), pattern="*.tsv", full.names=TRUE)
                                         }
                                           ),#kallisto/abundance.tsv
        "kallisto",
        tx2gene = t2g(),ignoreTxVersion = TRUE)$counts)#
        # Increment the progress bar, and update the detail text.
        progress$inc(2/2, detail = paste("Doing part", 2,"/",2))

        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        d
      }
    })

  #Option2: Upload count/expression table and annotation table from user

  #Input file containing count table/expression table as input from user
             output$file1<-renderUI({
                   
                   #validate(need(try(input$organism)!="","Please choose organism under study"))
               req(input$filechoice)
                   if(as.numeric(input$filechoice)==2)
                   {
                     fileInput(session$ns("file1"), label = h5('Choose Expression file to upload'),
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv'
                               )
                     )
                   }
                   
                 })
                 #Input file containing annotation file as input from the user
                 output$file2<-renderUI({
                   # validate(need(try(input$organism!=""),"Please choose organism under study"),
                   #          need(try(!is.null(input$filechoice)),"Please choose startng point of input"))
                   # if( input$organism!="" && !is.null(input$filechoice))
                   # {
                   req(input$filechoice)
                   print("line 194 ok")
                   #print(as.numeric(input$filechoice))
                   
                   
                   if(as.numeric(input$filechoice)==2)#!is.null(input$file1) && 
                   {
                     req(input$file1)
                     fileInput(session$ns("file2"), label = h5('Choose Annotation file to upload'),
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv'
                               )
                     )
                   }
                   else if( as.numeric(input$filechoice)==1)#input$kallisto_path!="" && 
                   {
                     req(input$kallisto_path)
                     fileInput(session$ns("file2"), label = h5('Choose Annotation file to upload'),
                               accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv',
                                 '.tsv'
                               )
                     )
                   }
                   # }
                   
                 })
   
              # })
 
  

  #compose count /expression table obtained as input to be displayed
  #edata reactive contains the count/expression table
  edata<-reactive({
    #setting up expression data
    #Get expression table from file
    if(as.numeric(input$filechoice)==2)
    {
      input_exp<-input_file(input$file1)
      if(!is.null(input_exp))
      {
        data<-input_exp[,-1]
        #convert entries to integers
        data=as.matrix(data)
        storage.mode(data)="integer"
        data= data.frame(data)
        rownames(data)=input_exp[,1]
        # print("before")
        # print(nrow(data))
        # idx.nz <- apply(data, 1, function(x) { all(x > 0)})
        # print("after")
        # print(head(idx.nz))
        # print(length(which(idx.nz)==T))
        # print(length(idx.nz))
        # sum(idx.nz)
        # data <- subset(data, idx.nz)
        #data <- data[ rowSums(data) > 1, ]#change
        data<-data
      }
    }
    else if(as.numeric(input$filechoice)==1)
    {
      dat<-txi()
      #convert entries to integers
      data=as.matrix(dat)
      storage.mode(data)="integer"
      data= data.frame(data)
      rownames(data)=rownames(dat)
      print("module 1 line 271 ok")
      print(head(data))
      data <- data[ rowSums(data) > 1, ]

    }

  })

  #Compose annotation data from file
  #pdata contains annotation file
  pData<-reactive({
    #setting up annotation data
    #Get annotation table from file
    pheno<-input_file(input$file2)
    print("module 1 line 285 ok")
    #print(pheno)

    #Get expression table
    data<-edata()
    #Expression table and annotation table should not be null
    if(!is.null(data) && !is.null(pheno))
    {

      #Assuming that the first column of annotation table is sample id ,extract all sample IDs
      sample_id = pheno[,1]
      #print(sample_id)
      #Get all sample IDs from expression table(sample ID refer to column names of expression table)
      exp_sample_id = colnames(data)
      #print('all')
     # print(exp_sample_id)
      #Check if all sample ID in expression table are present in the annotation table
      if (all(exp_sample_id %in% sample_id))
      {
        #print("Yay")

        #set all variables of annotation table as factors
        col<-1:ncol(pheno)
        for (i in col)
        {
          pheno[,i]<-as.factor(pheno[,i])
        }
        #Remove those sample IDs that are present in the expression table but absent in the annotation table from the annotation table
        idx<-NULL
        if(!all(sample_id %in% exp_sample_id))
        {
          idx<-which(!(sample_id %in% exp_sample_id))
          pheno <- pheno[-idx, ]
        }

        pheno<-pheno

      }
      else 
        {
          pheno<-NULL
          #print('nay')

        }
      pheno
    }
  })
out<-list(file2=reactive({input$file2}),
          edata=reactive({edata()}),
          pData=reactive({pData()}),
          organism=reactive({input$organism}))#c(1,2,3,4)}))
  
  #input$file2}))#,edata=reactive({edata()}),pData=reactive({pData()}))
return(out)
  }
 
  # })
  
