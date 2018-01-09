
heatmap_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    fluidRow(column(3,
                    selectInput(ns("Distance"), label = h5("Select distance method"), 
                                choices = list("Eucledian" = 1, "Manhattan" = 2),
                                selected = 1)
                    )),
    fluidRow(column(3,
                    selectInput(ns("Linkage"), label = h5("Select clustering method"), 
                                choices = list("Average" = 1,"complete"=2,"Ward.D2"=3,"Ward.D" = 4,"Singular"=5),
                                selected = 1)
                    )),
    uiOutput(ns("heat_comb")),
    downloadButton(ns('download_heatmap'), 'Download Plot'),
    fluidRow(column(8,plotOutput(ns("heatmap"))
                    )
             )
  )
}



heatmap_module<-function(input,output,session,heatmap_call,dds,batch,de_genes,combination,
                         wgcna_output)#batch_choice,batch_corrected)
{
  print("line 29 heatmap")
  output$heat_comb <- renderUI({
    if(heatmap_call!="ANOVA")
    {
      print("heatmap")
      result<-de_genes()
      print(result())
      print(is.null(result()))
      combo<-combination()
      if(length(combo())>0)
      {
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Making plot", value = 0)
        
        # Number of times we'll go through the loop
        n <- 2
        
        rows<-length(combo())
        modules<-NULL
        res<-data.frame(matrix(NA, nrow = length(combo()), ncol = 3))
        
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
        if(!is.null(wgcna_output()))
        {
          if(length(wgcna_output()$modules())>0)
          {
            modules<-as.data.frame(table(wgcna_output()$modules()))
            colnames(modules)<-c("Var1","number")
            entry<-c(combo(), levels(modules$Var1))
            print(entry)
            checklist<-list()
            for (i in seq_along(entry)) {
              checklist[[entry[[i]]]] = i
            }
            
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            
            selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                        choices = checklist,selected = 1)
            
          }
        }
        else{
          
          comb<-lapply(1:length(combo()), function(i) {
            combo()[[i]]
            #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
            
          })
          checklist<-list()
          for (i in seq_along(comb)) {
            checklist[[comb[[i]]]] = i
          }
          
          # Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
          
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          
          selectInput(session$ns("heat_choice"),label = h5("Choose comparison") ,
                      choices = checklist,selected = 1)
          
        }
        
      }
    }
    
  })
    print("mod heatmap")
    
        #heatmap DE
        #combinations

        heatmap<-reactive({
          
          dds.fc<-batch()$batch_data_for_DESeq()[[1]]
          rld<- assay(dds.fc)
          
          if(as.numeric(batch()$batch_choice())!=1)
          {
            rld<-batch()$batch_corrected()
          }
          print("rld")
          print(head(rld))
          if(heatmap_call=="ANOVA")
          {
            print("line 128 heatmap")
            combo<-combination()
            print(combo())
            num<-length(combo())
            heatmap_name<-"Heatmap of top 1000 most variable genes"
            heatmap_genes(heatmap_call,dds(),dds.fc,rld,NULL,
                          result,NULL,NULL,
                          num,heatmap_name,
                          input$Distance,input$Linkage)
          }
          else #if (heatmap_call=="DE")
          {
            result<-de_genes()
            print(head(result()))
            combo<-combination()
            print(combo())
            num<-length(combo())
            
            heatmap_name<-" "
            if(!is.null(wgcna_output()))
            {
              if(length(wgcna_output()$modules())>0)
              {
                mod<-wgcna_output()$modules()
                modules<-as.data.frame(table(mod))
                colnames(modules)<-c("Var1","numbers")
                entry<-c(as.vector(combo()), as.vector(modules$Var1))
                
                heatmap_name<-paste("Heatmap of ",heatmap_call,entry[as.numeric(input$heat_choice)])
                h<-heatmap_genes(heatmap_call,NULL,dds.fc,rld,input$heat_choice,
                              result,
                              mod,wgcna_output()$WGCNA_matrix(),
                              num,heatmap_name,
                              input$Distance,input$Linkage)
                #insert validate
                validate(need(!is.null(h),'Heat map cannot be dislayed for one gene/trascription factor'))
                h
              }
            }
            else
            {
              heatmap_name<-paste("Heatmap of ",combo()[[as.numeric(input$heat_choice)]])
              
              h<-heatmap_genes(heatmap_call,NULL,dds.fc,rld,input$heat_choice,
                            result,NULL,NULL,
                            num,heatmap_name,
                            input$Distance,input$Linkage)
              #insert validate
              validate(need(!is.null(h),'Heat map cannot be dislayed for one gene/trascription factor'))
              h
            }
          }
          
        })

        #Display heatmap as .png with dimension 980 X 680
        output$heatmap <- renderImage({
          # A temp file to save the output. It will be deleted after renderImage
          # sends it, because deleteFile=TRUE.
          outfile <- tempfile(fileext='.png')
          # Generate a png
          if(!is.null(input$heat_choice)){
            # Create a Progress object
            progress <- shiny::Progress$new()
            # Make sure it closes when we exit this reactive, even if there's an error
            on.exit(progress$close())
            
            progress$set(message = "Making plot", value = 0)
            
            # Number of times we'll go through the loop
            n <- 2
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
            png(outfile, width=980, height=680)
            heatmap()
            dev.off()
            # Increment the progress bar, and update the detail text.
            progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
            
            # Pause for 0.1 seconds to simulate a long computation.
            Sys.sleep(0.1)
          }
          else if(heatmap_call=='ANOVA')
          {
            png(outfile, width=980, height=680)
            heatmap()
            dev.off()
          }
          else
          {
            png(outfile, width=980, height=680)
            dev.off()
          }
          # Return a list
          list(src = outfile,
               alt = "This is alternate text")
          
        }, deleteFile = TRUE)
        
        
        #download heatmap
        output$download_heatmap <- downloadHandler(
          
          filename = function(){
            combo<-combination()
            choice<-as.numeric(input$heat_choice)
            heatmap_names<-""
            
            if(heatmap_call=="ANOVA")  heatmap_name<-"Heatmap of top 1000 most variable genes"
            else heatmap_name<-paste("Heatmap of ",heatmap_call,combo()[[choice]])
            paste(heatmap_name,'.pdf')
            
          },
          content = function(file) {
            
            png(file,width = 980, height = 680)
            # pdf(file)
             heatmap()
            #ggsave(file,heatmap())
            dev.off()
            
          })   
      
         
         # result<-NULL
         # 
         # result<-DE$de_genes()
         
         
         
return(list(
  input_Distance=reactive({input$Distance}),
  input_Linkage=reactive({input$Linkage})
))         
        
}
  
  


  
  