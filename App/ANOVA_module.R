
source("./gene_count_module.r")
ANOVA_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    radioButtons(ns("lfc"), "Show log2 foldchange?:", choices = c("Yes", "No"),selected = "No"),
    downloadButton(ns('download_ANOVA_Table'), 'Download full Data'),
    DT::dataTableOutput(ns("ANOVA")),
    gene_count_module_UI(ns("module5"))
    # uiOutput("plot_option_gene"),
    # downloadButton('download_gene_anova', 'Download full Data'),
    # plotOutput("gene_plot")
  )
}
ANOVA_module<-function(input,output,session,
                       batch_choice,dds,design,combination,
                       conchoice,DE_genes,normal,batch_corrected)
{
  if(!is.null(DE_genes()))
  {
    print("inside anova module")
    #print(batch)
    #Perform ANOVA 
    anova<-reactive({
      res<-NULL
      # batch_choice<-batch$batch_choice()
      # dds<-batch$batch_data_for_DESeq()[[1]]
      # design<-batch$batch_data_for_DESeq()[[2]]
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      
      progress$set(message = "Preparing ANOVA table", value = 0)
      
      # Number of times we'll go through the loop
      n <- 2
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      if(as.numeric(batch_choice())==2)
      {
        
        res<-DESeq(dds(),test = "LRT",reduced = design())#,parallel = TRUE)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      else if(as.numeric(batch_choice())==3)
      {
        #print(dds)
        print(design())
        #print(dds())
        res<-DESeq(dds(),test = "LRT",reduced = design())#,parallel = TRUE)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      else #if(as.numeric(batch_choice())==1)
      {
        res<-DESeq(dds(),test = "LRT",reduced = design())#,parallel = TRUE)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        
      }
      res
    })
    #prepare ANOVA table for display(This reactive merges the anova table with DE analysis results)
    anova_table<-reactive({
      
      dds_anova<-anova()
      res<-results(dds_anova)
      
      #Get base mean of each condition
      base_mean<-sapply( levels(dds_anova$condition), function(lvl) rowMeans( counts(dds_anova,normalized=TRUE)[,dds_anova$condition == lvl] ) )
      
      #Merge the basemeans with the anova table
      anova<-merge(res,base_mean,by=0,all=TRUE)
      for(i in 8:length(colnames(anova)))
      {
        colnames(anova)[i]<-paste0(colnames(anova)[i],'_mean')
      }
      print(head(anova))
      genes<-anova[,1]
      anova<-anova[-1]
      rownames(anova)<-genes
      colnames(anova)[1]<-"Overall mean"
      result<-DE_genes()
      combo<-combination()
      num <-length(combo())
      l<-list()
      l[[1]]<-anova[,-3]
      for(i in 1:num)
      {
        print(nrow(as.data.frame(result()[[i]][4])))
        res<-as.data.frame(result()[[i]][4])
        print(head(res))
        ######change####
        a_tab<-res[,-4]
        print(colnames(a_tab))
        c<-colnames(a_tab)
        print(length(c))
        print('howdy')
        print(c)
        res <- a_tab[,c(c[2],c[3],c[5],c[6],c[1],c[4])]
        ####################
        
        colnames(res)[5]<-"Overall mean"
        for (j in 1:length(colnames(res)))
        {
          colnames(res)[j]<-paste0(combo()[[i]],' ',colnames(res)[j])
        }
        
        l[[length(l)+1]]<-res
        
        #print(head(anova,1))
      }
      print(head(l[[1]]))
      #Condense results of ANOVA and DE analysis into one data
      anova<-transform(Reduce(merge, lapply(l, function(x) data.frame(rn = rownames(x),x))), row.names=rn, rn=NULL)
      
    })
    
    #library(DT)
    #Display ANOVA table
    output$ANOVA <- DT::renderDataTable({
      a_tab<-anova_table()[,-c(2,3)]
      ########change#######
      cond<-unique(colData(dds())[,as.numeric(conchoice())])
      print(cond)
      c<-colnames(a_tab)
      print(length(c))
      temp<-as.vector(c[4:(3+length(cond))])
      temp2<-as.vector(c[(length(cond)+4):length(c)])
      print(temp2)
      # #temp<-as.vector(c[7:length(c)])
      print(c(temp,c[2],c[3],temp2,c[1]))
      print('howdy')
      # print(c)
      # anova <- 
      print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
      
      ######################
      a_tab=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
      anova<-NULL
      if(input$lfc =="No")
      {
        names<-colnames(a_tab)
        idx<-grep("log2FoldChange",names)
        anova=a_tab[,-idx]
      }
      else anova=a_tab
      
      DT::datatable(anova,class = 'cell-border stripe',
                    selection = list(mode='single',target = 'row'),
                    extensions = list('Scroller'=NULL,'Buttons'=NULL),
                    options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 300,scroller = TRUE,dom = 'Bfrtip',
                                   buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),text = 'Download only genes on clipboard'))#I('colvis')
                    )
      )
      
    })
    #display boxplot when row(gene) in ANOVA table is clicked
    observeEvent(input$ANOVA_rows_selected,{
      print('hey')
      print(input$ANOVA_rows_selected)
      #print(input$ANOVA_rows_clicked)
      selected <- input$ANOVA_rows_selected
      print(selected)
      full_data<-normal()
      print(head(full_data))
      #get the gene
      print("howdy")
      print(head(anova_table()))
      #print(rownames(anova_table()))
      print(head(rownames(anova_table())),1)
      print(typeof(rownames(anova_table())))
      temp<-strsplit(rownames(anova_table())," ")
      print(head(temp))
      print(temp[[1]])
      an_gene<-temp[[selected]]
      print("hey")
      print(an_gene)
      library('data.table')
      print(an_gene)
      counts<-as.vector(full_data[which(an_gene %in% rownames(full_data)),])
      #print(counts)
      #print(as.factor(colData(dds.fc()[[1]])[,as.numeric(input$conchoice)]))
      cond<-as.vector(colData(dds())[,as.numeric(conchoice())])
      print(length(counts))
      print(length(cond))
      df<-data.frame(counts,cond)
      colnames(df)<-c('count','condition')
      print(head(df))
      print(batch_choice())
      print(as.numeric(batch_choice()))
      batch<-batch_choice()
      if(as.numeric(batch)==1) callModule(gene_count_module,"module5",NULL,reactive({dds()}),reactive({an_gene}))
      else
        {
          idx<-which(rownames(batch_corrected()) %in% an_gene)
          print(idx)
          callModule(gene_count_module,"module5",batch_corrected()[idx,],reactive({dds()}),reactive({an_gene}))
        }
    })
    #Download ANOVA table
    output$download_ANOVA_Table <- downloadHandler(
      
      filename = function() { 'ANOVA_table.xlsx' },
      content = function(file) {
        
        anova<-anova_table()
        a_tab<-anova_table()[,-c(2,3)]
        ########reordering columns in anova table for display#######
        cond<-unique(colData(dds())[,as.numeric(conchoice())])
        print(cond)
        c<-colnames(a_tab)
        print(length(c))
        temp<-as.vector(c[4:(3+length(cond))])
        temp2<-as.vector(c[(length(cond)+4):length(c)])
        print(temp2)
        # #temp<-as.vector(c[7:length(c)])
        print(c(temp,c[2],c[3],temp2,c[1]))
        print('howdy')
        # print(c)
        # anova <- 
        print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
        
        ######################
        anova=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
        #write.csv(anova, file)
        write.xlsx2(anova, file, sheetName = "ANOVA",
                    col.names = TRUE, row.names = TRUE, append = FALSE)
      }
    )
  }
 
return(list(
  anova_table=reactive({anova_table()}),
  dds=reactive({anova()})
)
  )  
}
