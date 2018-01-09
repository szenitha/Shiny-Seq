source("./functions.r")
source("./PCA.r")
source("./gene_count_module.r")
Module_Batch_effect_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    #Step 1: creates a checkbox widget to select the batch effect
    #is batch effect present? if so then is it due to a known variable/unknown variable?
    #This widget gets the input from the user
    fluidRow(column(4,
                    radioButtons(ns("batch_choice"), "Choose the kind of batch effect",
                                 choices = list("no batch effect"=1,"known batch"=2,"hidden batch/unknown batch"=3),
                                 selected = 1))),
    conditionalPanel(sprintf("input['%s']==2", ns("batch_choice")),
                      fluidRow(column (4,uiOutput(ns("batch")))),
                      actionButton(ns("helperr2"), "Help")),
    conditionalPanel(sprintf("input['%s']==3", ns("batch_choice")),
                     fluidRow(column (4,uiOutput(ns("text")))),
                     actionButton(ns("go"), "Go")
                     
                     ),
    actionButton(ns("PCA"),"PCA"),
    fluidRow(column(10,
                    # batch corrected data table
                    downloadButton(ns('download_batch_corrected_Data'), 'Download full Data'),
                    DT::dataTableOutput(ns("batch_corrected")),
                    gene_count_module_UI(ns("module1"))
                    #uiOutput("plot_option_gene_bc"),
                    # downloadButton('download_gene_bc', 'Download plot'),
                    # plotOutput("gene_plot_bc")
             )),
    bsModal(ns("modalExample"), "Batch corrected PCA", ns("go"),size = "large",#uiOutput("plots")),
            selectInput(ns("top_pca_sva"), label = h5("Should PCA be computed for:"), 
                        choices = list("Top 500 genes" = 1, "All" = 2),
                        selected = 2),
            fluidRow(textOutput(ns("plottitle1")),column(10,plotlyOutput(ns("plot1")))),#,height = "400px", width = "100%")),
            fluidRow(textOutput(ns("plottitle2")),column(10,plotlyOutput(ns("plot2")))),
            fluidRow(textOutput(ns("plottitle3")),column(10,plotlyOutput(ns("plot3")))),
            fluidRow(textOutput(ns("plottitle4")),column(10,plotlyOutput(ns("plot4")))),
            fluidRow(textOutput(ns("plottitle5")),column(10,plotlyOutput(ns("plot5")))),
            fluidRow(textOutput(ns("plottitle6")),column(10,plotlyOutput(ns("plot6"))))),
 
  bsModal(ns("modalpca"),"batch corrected PCA",ns("PCA"),size="large",
             PCA_UI(ns("pca_module_batch"))),
  bsModal(ns("modalhelperr2"), "Help for error model matrix not full rank", ns("helperr2"),size = "large",
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
          )
  
          )
  
}

Module_Batch_effect<-function(input,output,session,conchoice,designchoice,dds.fc,normal)
{
  print("inside module 3 line 18")
 # dds.fc<-dds()
  print(colData(dds.fc()[[1]]))
  #Step 2: If user selects option 2 in the above widget. It means that one/two variables in the
  #annotation variable is causing the batch effect
  #Hence the below checkbox widget is created to enable user to select the variable from annotation 
  #table known to cause batch effect.This variable is identified by visualizing the PCA of normalized data in the previous step
  output$batch <-
    renderUI({
      req(input$batch_choice)
      if(as.numeric(input$batch_choice)==2) {
        print("create Checkbox")
        checklist = list()
        for (i in seq_along(colnames(colData(dds.fc()[[1]])))[-1]) { 
          if(as.numeric(i)!=as.numeric(conchoice)){
            checklist[[colnames(colData(dds.fc()[[1]]))[[i]]]] = i}
        }
        selectizeInput(session$ns("batch"), "Choose the batch variable", checklist,selected = NULL, options = list(maxItems = 2))
      }
    })
  #Step 3: If user selects option 3 in step 1, then it is an indication that an unknown variable is
  #causing the batch effect. The package used in this application to identify unknown/hidden batch effects is SVA
  #There is another package called RUVSeq which is also designed for the same purpose. Only SVA is supported here.
  #However you can check out the package RUVSeq in bioconductor. Since SVA does not need housekeeping genes
  #to look for unknown batch effects and the fact that the results of SVA and RUVSeq is same, we use SVA here.
  #The application identifies the unknown vaiables by constructing from the dataset
  #The constructed variables are called as surrogate variables. The total number of surrogate variables are
  #computed in the below step.
  #Calculates the total number of surrogate variables
  num_sv<-reactive({
    #Check the vignette of package SVA to get a better understanding of computing the surrogate variables
    
    #Load package SVA
    library('sva')
    #Inputs neede to compute the number of surrogate variables are
    #Input 1: Normalized data 
    #Input 2:Full model
    
    #Input 1:
    #Get rlog data 
    rld <- assay(dds.fc()[[2]])#normal()
    data_normalized<-normal()
    #print(ncol(data_normalized))
    
    #Input 2:
    #compute full model
    mod =model.matrix(design(dds.fc()[[1]]), data=colData(dds.fc()[[1]]))
    #compute null modal
    mod0 = model.matrix(~1,data=colData(dds.fc()[[1]]))
    
    #Get the total number of surrogate variables
    n.sv = num.sv(data_normalized,mod,method="leek")# use method="be" if u want less no. of surrogate variables 
    #get surrogate variables
    svaobj <-sva(rld, mod, mod0, n.sv=n.sv-2)#n.sv)
    print("inside module3 line 86")
    #print(head(svaobj$sv))
    n.sv<-dim(svaobj$sv)[2]
    list(n.sv,rld,mod,mod0,svaobj)
  })
  
  #Displays the total number of surrogate variables to the user as a text.
  output$text<-renderUI({
    if(as.numeric(input$batch_choice)==3) {
      print(paste('total number of sv',num_sv()[[1]]))
      n.sv<-num_sv()[[1]]
      if(n.sv==0) h3('No hidden batch') #If there are no surrogate variables (Means no batch effec)
      else
      {
        print('hi')
        print(n.sv)
        x<-paste("Total number of surrogate variables is ",n.sv)
        print(x)
        if(n.sv==1) textInput(session$ns("num"),label = list(h4(x),h4("Enter number of surrogate values to be removed")), 
                              value = 1) #If total number of surrogate variables is one
        
        else if(n.sv>1) #If the total number of surrogate variables is more than one 
        {               
          if( n.sv<6) #If total number of surrogate variables is between 1-6. the we set the default number of
          {           #Surrogate variables to be removed as 1. It is adviced not to remove too many surrogate variables as it could
                      #result in removing variation due to biological variable of interest. In this case 1 is good.
            selectInput(session$ns("num"),label = list(h5(x),h5("Enter number of surrogate values to be removed") ),
                        choices = 0:n.sv,selected = 1)#,min=0,max=n.sv)  
          }#If total number of surrogate variables is more than 6 then a max of 3 surrogate variables is adviced to be removed
          #henece we set the default as 3
          else selectInput(session$ns("num"),label = list(h5(x),h5("Enter number of surrogate values to be removed") ),
                           choices = 0:n.sv,selected = 3)#,min=0,max=n.sv)  
        }
      }
    }
    
    
  })

  # #This reactive adjusts the DESeq2 object considering batch variable
  # #Ideally the design considered id ~condition but if batch effect is present then ~batch+condition
  # #check DESeq2 manual for more details
  batch_design<-reactive({
    #batch_design<-NULL
    req(input$batch_choice)
    if(as.numeric(input$batch_choice)==2){ #If variables in annotation table have been identified to cause batch effect
      req(input$batch)
      library(formula.tools)
      #get
      full_design<-design(dds.fc()[[1]]) #get current design (default: ~condition ) unless design choice provided by user

      print(as.numeric(input$designchoice) == as.numeric(input$batch))
      
      if(length(input$batch)==1) #If only one known variable is causing batch effect
      {
        if(!(as.numeric(input$batch) %in% as.numeric(input$designchoice))) #Check if the variable was already specified by user by selcting design choice
        { #If true then add batch variable to design. New design formula will be ~batch variable + condition
          print('heya')
          var1<-substring(as.character(design(dds.fc()[[1]])),1,1)
          #print(var1)
          var2<-substring(as.character(design(dds.fc()[[1]])),2)
          #print(var2)
          full_design<-paste(var1,colnames(colData(dds.fc()[[1]])[as.numeric(input$batch)]),' + ',var2,sep = '')

          full_design<- formula(full_design)
  
        }
        #print(full_design)
        #Here we construct design formula only containing design (~batch variable)
        #This value is used during ANOVA
        only_batch_design<-paste('~',colnames(colData(dds.fc()[[1]])[as.numeric(input$batch)]),sep = '')
        data<-assay(dds.fc()[[1]])
        pheno<-colData(dds.fc()[[1]])
        
        print(full_design)
        #set order of columns in expression data as same as order of sample  ID in pheno data
        #data=data[,as.vector(pheno[,1])]
        dds.fc_local<-DESeqDataSetFromMatrix(data,colData=pheno,design = full_design)
        rld<-rlogTransformation(dds.fc_local) #Perform rlog tranformation. This format is an input to PCA
        #print(dds.fc)
        #Return (DESeq2 datatset,batch design formula containing only batch variable,rlog transformed data)
        #DESeq2 dataset contains count data, annotation table and design formula(~condition+batch(if present))
        list(dds.fc_local,formula(only_batch_design),rld)

      }

      else if(length(input$batch)==2) #If two known variables cause batch effect
      {
        print('heya')
        var1<-substring(as.character(design(dds.fc()[[1]])),1,1)
        #print(var1)
        var2<-substring(as.character(design(dds.fc()[[1]])),2)
        #print(var2)
        
        #Check which of these variables were already selected by user as design choice
        idx<-which(!(as.numeric(input$batch) %in% as.numeric(designchoice)))
        print(idx)
        #if one of these variables were selected as design choice. In this case we only add variable not present to design formula.
        if(length(idx)==1)
        {
          full_design<-paste(var1,colnames(colData(dds.fc()[[1]])[as.numeric(input$batch[idx])]),' + ',var2,sep = '')

          full_design<- formula(full_design)
          print("line 253")
          print(full_design)
        }
        #If none of these variables were specified in design choice. In this case we add both variables to design formula.
        else if(length(idx)==2)
        {
          full_design<-paste(var1,colnames(colData(dds.fc()[[1]])[as.numeric(input$batch[1])]),' + ',
                             colnames(colData(dds.fc()[[1]])[as.numeric(input$batch[2])]),'+',
                             var2,sep = '')

          full_design<- formula(full_design)
        }

        temp<-paste('~',colnames(colData(dds.fc()[[1]])[as.numeric(input$batch[1])]),sep = '')
        only_batch_design<-paste(temp,colnames(colData(dds.fc()[[1]])[as.numeric(input$batch[2])]))
        data<-assay(dds.fc()[[1]])
        pheno<-colData(dds.fc()[[1]])
        print('hey2')
        print(full_design)
        #set order of columns in expression data as same as order of sample  ID in pheno data
        #data=data[,as.vector(pheno[,1])]
        #Check if the batch variables are confounded with biological variables of interest.
        #If so the throw error.
        tryCatch(dds.fc<-DESeqDataSetFromMatrix(data,colData=pheno,design = full_design)
                 ,error = function(c) {
                   print(typeof(c$message))
                   print(substr(c$message,1,33))
                   if("the model matrix is not full rank" == substr(c$message,1,33))  c$message <- paste0("The model matrix is not full rank. To know more about this please click the help button on this page.")
                   else c$message
                   stop(c)
                 })
        rld<-rlogTransformation(dds.fc)
        #print(dds.fc)
        #return (DESeq2 dataset, design formula containing batch variable,rlog tranformed data)
        list(dds.fc,formula(only_batch_design),rld)
      }

    }
    #If hidden/unknown batch effect present. SVA package is used. For more information please check vignette of package
    else if(as.numeric(input$batch_choice)==3)# && (!is.null(as.numeric(input$num))))
    {
      #First we load SVA package
      library('sva')
     
      #get the number of surrogate variables (Indicated by user. Once the total number of surrogate variables
      #are displayed, user selected the number of surrogate variables to be removed) we take this variable. 
      
      #Input 1: Assign number of surrogate variables provided by user
      n.sv=as.numeric(input$num) 
      print(n.sv)
      #Is number of surrogate variables > 0?
      if(n.sv>0){
        
        #To compute the surrogate variables, we need the following inputs
        
        
        # Input 1: rlog transformed data 
        rld<-num_sv()[[2]]
        # Input 2: compute full model
        mod<-num_sv()[[3]]
        # Input 3: Compute null model
        mod0<-num_sv()[[4]]
        
        #We now use the sva function twice to get the surrogate variables. Each one is computed for a specific purpose
        #purpose 1: We want to visualize the PCAs of expression data with 1-6 surrogate variables removed, provided
        # the total number of surrogate variables compted is >6. If the total number of surrogate variables is <6 then
        # we take the PCAs of expression data with all surrogate variables removed
                svall<-NULL
        #if total number of surrogate variables is more than 6 then pick only first 6 surrogate variables 
        if(num_sv()[[1]] >6) svall<-sva(rld, mod, mod0, n.sv= 6)#4)#change 6)
        else svall<-sva(rld, mod, mod0, n.sv= num_sv()[[1]]) #else pick total
        print(svall$sv)
        
        #purpose 2: We compute the expression data with the user specified number of surrogate variables removed
        #Get surrogate variables already computed in reactive n.sv
        svaobj <-num_sv()[[5]]
        print("inside module 3 line 262")
        print(head(svaobj$sv[,1:n.sv]))
        #svaobj<-svaseq(data_normalized,mod,mod0,n.sv = n.sv)
        sv<-as.data.frame(svaobj$sv[,1:n.sv])
         
        #prepare DESeq2 dataset
        data<-assay(dds.fc()[[1]])#edata
        # print('data')
        # print(head(data))
        pData<-colData(dds.fc()[[1]])#annotation
        # print('pdata')
        # print(pData)
        
        #Here we need to add surrogate variables to design formula of DESeq2 dataset
        #This is done to consider the surrogate variables as covariates in the differential expression
        #analysis (new design formula (~Surrogate variable (1-values specified by user)+condition))
        var1<-substring(as.character(design(dds.fc()[[1]])),1,1)[1]#as.character(design(dds.fc()[[1]]))[1]
        print("line298")
        print(as.character(design(dds.fc()[[1]])))
        print(var1)
        var2<-substring(as.character(design(dds.fc()[[1]])),2) #as.character(design(dds.fc()[[1]])) [2]
        rownames(sv)<-pData[,1]
        temp<-''
        for (i in 1:n.sv){
          temp<-paste(temp,'SV',i,' + ',sep = '')
          colnames(sv)[i]<-paste('SV',i,sep = '')
        }
        print(temp)
        print(substr(temp, 1, nchar(temp)-2))
        print(nchar(temp))
        full_design<-''
        batch_design<-''
        if(n.sv==1)
        {
          full_design<-paste('~SV1 + ',var2,sep = '')
          batch_design<-'~SV1'
        }
        else
        {
          print("line 297")
          print(var1)
          print(var2)
          print(temp)
          full_design<-paste(var1,temp,var2,sep = '')#,temp)
          print("full")
          print(full_design)
          batch_design<-paste('~',substr(temp, 1, nchar(temp)-2),sep=' ')
        }
        print(full_design)
        print(batch_design)
        print(rownames(sv))
        # print('merge')
        # print(pData)
        # print(sv)
        #add surrogate variables to current annotation
        pheno<-merge(pData,sv,by=0,all=T,sort=FALSE)[,-1]
        print('sva')
        #print(head(data))
        #print(pheno)
        #set order of columns in expression data as same as order of sample  ID in pheno data
        #data=data[,as.vector(pheno[,1])]
        dds.fc<-DESeqDataSetFromMatrix(data,colData=pheno,design = formula(full_design))
        rld<-rlogTransformation(dds.fc)
        #print(colData(rld))
        svobj<-NULL
        if(num_sv()[[1]]!=1) svobj<-svaobj$sv[,1:n.sv]
        else svobj<-svaobj$sv
        #return(DESeq2 dataset,design formula with batch design, rlog data,full modal,
        #Object containing all user selected surrogate varaibles,Only first 1-6 surrogate variables) 
        #
        list(dds.fc,formula(batch_design),rld,mod,svobj,svall$sv)#~SV1+SV2...+condition
      }

    }
    else if(as.numeric(input$batch_choice)==1){#if no batch effect option is selected
      list(dds.fc()[[1]],formula('~1'),dds.fc()[[2]])
      #return (DESeq2 dataset,design formula ~1(for ANOVA),rlog data)
      #DESeq2 and rlog data that was already computed in normalization step
    }
  })
  observeEvent(input$help, {
    toggleModal(session, "modalqc", toggle = "close")
  })

  #reactive calculates and returns the batch corrected table
  batch_corrected<-reactive({
    req(input$batch_choice)
    #req(input$batch)
    withProgress(message = 'Processing:',
                 detail = 'Removing batch effect .This may take a while...',value = 0, {
                   # Number of times we'll go through the loop
                   n <- 3
                   #if known variable is identified to cause batch effect
                   if(as.numeric(input$batch_choice)==2)#
                     {    
                     req(input$batch)#&&(!is.null(input$batch))){
                     dds<-dds.fc()[[1]]#batch_design()[[1]]
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 1))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     #Input 1: get rlog transformed data
                     rld<-dds.fc()[[2]]
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 2))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     #We remove  batch effect caused by known variable using LIMMA package
                     #Load LIMMA package
                     library('limma')
                     
                     #Step 2: To remove batch: Get which variable is batch
                     batch<-colData(rld)[,as.numeric(input$batch)]
                     batch2<-NULL 
                     #If there are two batch variables
                     if(length(input$batch)==2)
                     {
                       batch<-colData(rld)[,as.numeric(input$batch[1])] #Get batch variable 1
                       batch2<-colData(rld)[,as.numeric(input$batch[2])] #get batch variable 2
                     }
                     #Prepare input for batch removal
                     #Get condition/variable of interest
                     condition<-colData(rld)[,as.numeric(conchoice)]
                     print(condition)
                     print(model.matrix(~condition))
                     #remove batch effect
                     #Inputs provided: rlog data,batch variable 1, batch variable 2(if present),condition model matrix
                     batch_coorected<-removeBatchEffect(as.matrix(assay(rld)),batch = batch,batch2=batch2,design = model.matrix(~condition))

                     #Get most variable genes. This is done as an input preparation for PCA of batch corrected data
                     #500 most variable genes
                     rv <- rowVars(assay(rld))
                     genes<-order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]#top 500 genes or all
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 3))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     #return (batch corrected table(expression table-batch effect), top most variable genes)
                     list(batch_coorected,genes)
                     #print(head(r))
                     #r

                   }
                   else if(as.numeric(input$batch_choice)==3)# && as.numeric(input$num)>0)#(!is.null(input$num)))
                   {
                     req(input$num)
                     #Get DESeq2 dataset computed in normalization step
                     dds<-dds.fc()[[1]]
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 1))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     #Step 1: get rlog transformed data
                     rld<-dds.fc()[[2]]#rlogTransformation(dds)#,blind = FALSE)
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 2))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     #print(head(assay(rld)))
                     
                     #get the number of surrogate variables to be removed
                     mod<-batch_design()[[4]] #Get full design
                     svobj<-batch_design()[[5]] #Get object containing all user selected surrogate variables
                     svall<-batch_design()[[6]] #Get 1 to 6 surrogate variables
                     #print(svall)
                     print("hey")
                     #print(svall[,1:3])
                     
                     #Compute top 500 most variable genes for PCA of batch corrected data
                     rv <- rowVars(assay(rld))
                     genes<-order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]
                     clean_list<-list()
                     
                     #perform batch correction on rlog data removing 1-6 surrogates variables
                     #if total number of surrogate variables is <6
                     #If total number of surrogate variables is >6 then perform batch correction
                     #on rlog data by removingfirst 6 surrogate variables
                     
                     temp<-6
                     if(dim(as.matrix(svall))[2]<6) temp<-dim(as.matrix(svall))[2]
                     if(temp == 1) clean_list[[length(clean_list)+1]]<-cleanY(assay(rld),mod,svall)
                     else
                     {
                       for(i in 1:temp)#change
                       {
                         clean_list[[length(clean_list)+1]]<-cleanY(assay(rld),mod,svall[,1:i])
                       }
                     }
                     print('clean_list')
                    # print(head(clean_list))
                     # Increment the progress bar, and update the detail text.
                     incProgress(1/n, detail = paste("Doing part", 3))

                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     #Return (batch corrected data with user specified number of surrogate variables,
                     #top 500 most variable gene list, batch corrected data of fitst 6 or 'n' surrogate
                     #variables provided n<6)
                     list(cleanY(assay(rld),mod,svobj),genes,clean_list)
                     #genes<-select(500,rld)

                     #cleany =  cleanY(assay(rld),mod,svobj$sv[,1:2])
                   }
                 })
  })
  # Display the batch corrected table
  output$batch_corrected <- DT::renderDataTable({

    DT::datatable(batch_corrected()[[1]],class = 'cell-border stripe',
                  selection = list(mode='single',target = 'row'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 250,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download only genes on clipboard'))#I('colvis')
                  )
    )
  })

  #download batch corrected table
  output$download_batch_corrected_Data <- downloadHandler(
    filename = function() { 'batch_corrected.xlsx' },
    content = function(file) {
      #write.csv(batch_corrected()[[1]], file)
      write.xlsx2(batch_corrected()[[1]], file, sheetName = "Batch corrected data",
                  col.names = TRUE, row.names = TRUE, append = FALSE)
    }
  )
  #####gene plot count for batch corrected table (boxplot for batch corrected table)
  #display gene expression across conditon
  observeEvent(input$batch_corrected_rows_selected,{
    #get which gene was clicked
    #get the row  number of the gene clicked
    print('hey')
    print(input$batch_corrected_rows_selected)
    #print(input$ANOVA_rows_clicked)
    selected_row <- input$batch_corrected_rows_selected
    print(selected_row)
    #get the normalized data
    full_data<-normal()
    print(head(full_data))
    print('hey')
    
    #get the table where we clicked the gene
    result<-batch_corrected()[[1]]
    print(head(result))
    #get the gene
    an_gene<-rownames(result)[selected_row]
    print(an_gene)
    counts<-as.vector(full_data[which(an_gene %in% rownames(full_data)),])
    #print(counts)
    #print(as.factor(colData(dds.fc()[[1]])[,as.numeric(input$conchoice)]))
    cond<-as.vector(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
    #rownames(cond)<-colData(dds.fc()[[1]])[,1]
    library('data.table')
    print(length(counts))
    print(length(cond))
    df<-data.frame(counts,cond)
    colnames(df)<-c('count','condition')
    print(head(df))
    
    if(as.numeric(input$batch_choice)==1) callModule(gene_count_module,"module1",NULL,reactive({dds.fc()[[1]]}),reactive({an_gene}))
    else callModule(gene_count_module,"module1",result[selected_row,],reactive({dds.fc()[[1]]}),reactive({an_gene}))
  
    })
  # 
  # 
  # ############################
  # # Insert the right number of plot output objects into the web page
  # #this section is activated when batch effect tab--->click go button
  output$plots <- renderUI({

    plot_output_list <- lapply(1:6, function(i) {
      plotname <- paste("plot", i, sep="")
      plottitle <- paste("plottitle", i, sep="")
      tags$div(class = "group-output",
               textOutput(session$ns(plottitle), container = h6),
               plotlyOutput(session$ns(plotname),height = 50, width = "25%")
      )

    })
      #Convert the list to a tagList - this is necessary for the list of items
      #to display properly.
      do.call(tagList, plot_output_list)
    })
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  #Display the PCA after removeing  1 to 6 surrogate variables
  for (i in 1:6) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")
      plottitle <- paste("plottitle", my_i, sep="")
      output[[plottitle]] <- renderText({
        if(num_sv()[[1]]>6)
        {
          paste("No. of surrogates removed ", my_i, sep = "")
        }
        else if(num_sv()[[1]]<=6)
        {
          print('hey')
          # if(num_sv()[[1]]>=2)
          # {
            if(my_i<num_sv()[[1]])
            {
              paste("No. of surrogates removed ", my_i, sep = "")
            }
          #}
        }
      })
      output[[plotname]] <-renderPlotly({
        p<-plotly_empty()
        if(num_sv()[[1]]>6)
        {
          # Create a Progress object
          progress <- shiny::Progress$new()
          # Make sure it closes when we exit this reactive, even if there's an error
          on.exit(progress$close())

          progress$set(message = "Making plot", value = 0)

          # Number of times we'll go through the loop
          prog <- 3
         
          #Get batch corrrected and annotation table
          data<-batch_corrected()[[3]][[my_i]]
          genes<-batch_corrected()[[2]]
          pheno<-colData(dds.fc()[[1]])#batch_design()[[1]])
          # Increment the progress bar, and update the detail text.
          progress$inc(1/prog, detail = paste("Doing part", 1,"/",prog))

          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          #print(head(data,1))
          #Get which column u want to plot by
          #var<- as.numeric(conchoice)#as.numeric(input[[paste("batch_corrected_pca",i)]])
          dat<-list(assay(dds.fc()[[1]]),pheno,data)
          p<-pcaplot(dat,conchoice,input$top_pca_sva,"3D",resize_factor=1)[[1]]

          # Increment the progress bar, and update the detail text.
          progress$inc(1/prog, detail = paste("Doing part", 2,"/",prog))

          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)

          # Increment the progress bar, and update the detail text.
          progress$inc(1/prog, detail = paste("Doing part", 3,"/",prog))

          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          p
        }
        else if(num_sv()[[1]]<=6)
        {
          if(num_sv()[[1]]>=2)
          {
            if(my_i<num_sv()[[1]])
            {
              print('my_i')
              print(my_i)
              # Create a Progress object
              progress <- shiny::Progress$new()
              # Make sure it closes when we exit this reactive, even if there's an error
              on.exit(progress$close())

              progress$set(message = "Making plot", value = 0)

              # Number of times we'll go through the loop
              prog <- 3
             
              #Get batch corrrected and annotation table
              data<-batch_corrected()[[3]][[my_i]]
              genes<-batch_corrected()[[2]]
              pheno<-colData(dds.fc()[[1]])#batch_design()[[1]])
              # Increment the progress bar, and update the detail text.
              progress$inc(1/prog, detail = paste("Doing part", 1,"/",prog))

              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
              #print(head(data,1))
              #Get which column u want to plot by
              #var<-as.numeric(conchoice)#as.numeric(input[[paste("batch_corrected_pca",i)]])
              dat<-list(assay(dds.fc()[[1]]),pheno,data)
              p<-pcaplot(dat,conchoice,input$top_pca_sva,"3D",resize_factor=1)[[1]]
              
              # Increment the progress bar, and update the detail text.
              progress$inc(1/prog, detail = paste("Doing part", 2,"/",prog))

              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
             
              
              # Increment the progress bar, and update the detail text.
              progress$inc(1/prog, detail = paste("Doing part", 3,"/",prog))

              # Pause for 0.1 seconds to simulate a long computation.
              Sys.sleep(0.1)
            }
          }
          p
        }
        #p

      })
    })
  }
 b<-reactiveValues(top=NULL)
observeEvent(input$PCA,
             {
               batch_pca<-normal()
               if(as.numeric(input$batch_choice)>1) batch_pca<-batch_corrected()[[1]]
               batch<-list(assay(dds.fc()[[1]]),colData(dds.fc()[[1]]),batch_pca)
               b$top<-callModule(PCA,"pca_module_batch",batch,resize_factor=1)
             })  
#req(input$batch_choice)

return(list(batch_data_for_DESeq=reactive({batch_design()}),
            batch_corrected_data=reactive({batch_corrected()[[1]]}),
            batch_choice=reactive({input$batch_choice}),
            top_batch_pca=reactive({b$top})
            ))
}









# 
# 

