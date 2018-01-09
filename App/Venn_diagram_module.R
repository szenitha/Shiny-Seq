
Venn_diagram_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    DT::dataTableOutput(ns("de_venn")),
    downloadButton(ns('download_venny'), 'Download plot'),
    plotOutput(ns("venn_diag"))
   
  )
}

Venn_diagram_module<-function(input,output,session,DE_genes,
                                      combination,wgcna_output)
{
  #combination()
  combo<-combination()
  print(combo())
  num<-length(combo())
  
  result<-DE_genes()
 
#######venn diagram#####
#display interactive table that summarizes DE genes identified for comparisons
output$de_venn <- DT::renderDataTable({
  
    rows<-num
    modules<-NULL
    WGCNA_matrix<-NULL
    res<-data.frame(matrix(NA, nrow = num, ncol = 3))
    
    if(!is.null(wgcna_output()))
    {
      if((length(wgcna_output()$modules())>0))
      {
        mod<-wgcna_output()$modules()
        modules<-as.data.frame(table(mod))
        colnames(modules)<-c("Var1","numbers")
      # print(modules)
      # WGCNA_matrix<-wgcna()[[2]]
        entry<-c(combo(), levels(modules$Var1))
      # print(modules$Var1)
        print(entry)
        rows<-length(entry)
        res<-data.frame(matrix(NA, nrow = rows, ncol = 3))
        colnames(res)<-c('Up regulated','Down regulated','Both')
      
      # entry<-c(as.vector(input$combination), as.vector(modules$Var1))
      # print(modules$Var1)
      # print(entry)
        rownames(res)<-lapply(1:rows, function(i) {
         entry[[i]]
        
        })
       for(i in 1:num)
       {
        #print(nrow(as.data.frame(result()[[i]][1])))
        res[i,1]<-nrow(as.data.frame(result()[[i]][1]))
        res[i,2]<-nrow(as.data.frame(result()[[i]][2]))
        res[i,3]<-nrow(as.data.frame(result()[[i]][3]))
       }
       for(i in num+1:nrow(modules))
       {
        #print(nrow(as.data.frame(result()[[i]][1])))
        print('res')
        print(result()[[i]][3])
        res[i,1:2]<-0
        res[i,3]<-result()[[i]][3]
       }
      }
    }
    else{
      rownames(res)<-lapply(1:num, function(i) {
        combo()[[i]]
        #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
        
      })
      colnames(res)<-c('Up regulated','Down regulated','Both')
      for(i in 1:num)
      {
        #print(nrow(as.data.frame(result()[[i]][1])))
        res[i,1]<-nrow(as.data.frame(result()[[i]][1]))
        res[i,2]<-nrow(as.data.frame(result()[[i]][2]))
        res[i,3]<-nrow(as.data.frame(result()[[i]][3]))
      }
    }
    print(res)
    DT::datatable(res,class = 'cell-border stripe',
                  selection = list(target = 'cell'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),
                                                             text = 'Download table'))),#I('colvis')
                  
                  escape = FALSE
                  
    )
})
#$(this).toggleClass('selected');
#Shiny.onInputChange('cells',
#table.cell('.selected').data().toArray());

#Display the venn diagram when user selects the sets to get the overlap of the sets
observeEvent(input$de_venn_cell_clicked,{
  print('hey')
  print(input$de_venn_cells_selected)
  print(input$de_venn_cell_clicked)
  selected <- input$de_venn_cells_selected
  print(nrow(selected))
  venn_list<-list()
  #result()<-DE_genes()
  
  #num <- length(input$combination)
  comp<-lapply(1:num, function(i) {
    combo()[[i]]
    
  })
  
  if(!is.null(wgcna_output()))
  {
    if((length(wgcna_output()$modules())>0))
    {
      mod<-wgcna_output()$modules()
      modules<-as.data.frame(table(mod))
      colnames(modules)<-c("Var1","numbers")
      #WGCNA_matrix<-wgcna_output()$WGCNA_matrix()
      entry<-c(combo(), levels(modules$Var1))
      print(entry)
      rows<-length(entry)
      entry<-c(as.vector(combo()), as.vector(modules$Var1))
      comp<-lapply(1:rows, function(i) {
        entry[i]
    })
    }
  }
  reg<-c('up regulated genes','down regulated genes','All DE genes')
  pal<-c('pink','light green','light blue','light yellow')
  col<-c()
  circ<-c()
  venn_names<-list()
  if(nrow(selected) %in% 2:4)
  {
    for(i in 1:nrow(selected))
    {
      print(selected[i])
      row<-selected[i,1]
      print('row')
      print(row)
      col<-selected[i,2]
      print('col')
      print(col)
      print('hey')
      genes<-NULL
      if(row>num)
      {
        # mod<-wgcna_output()$modules()
        # modules<-as.data.frame(table(mod))
        # colnames(modules)<-c("Var1","numbers")
        mod<-wgcna_output()$modules()
        modules<-as.data.frame(table(mod))
        colnames(modules)<-c("Var1","numbers")
        WGCNA_matrix<-wgcna_output()$WGCNA_matrix()
        print(modules$Var1[row-num])
        #print(head(colnames(WGCNA_matrix)))
        idx_w<-which(mod==modules$Var1[row-num])
        print(head(idx_w))
        genes<-colnames(WGCNA_matrix)[idx_w]
        print(head(as.data.frame(genes)))
      }
      else
      {
        df<-as.data.frame(result()[[row]][col])
        genes<-df[,1]
        rownames(df)<-genes
      }
      venn_list[[length(venn_list)+1]]<-genes
      venn_names[[length(venn_names)+1]]<-paste(comp[row]," ",reg[col])
      
    }
    alpha<-c()
    i<-nrow(selected)
    if(i==2)
    {
      col<-c(pal[1],pal[2])
      circ<-c(pal[1],pal[2])
      alpha<-c(0.5,0.5)
    }
    else if(i==3)
    {
      col<-c(pal[1],pal[2],pal[3])#c('lightpink1','lightgreen','lightgoldenrod1')
      circ<-c(pal[1],pal[2],pal[3])
      alpha<-c(0.5,0.5,0.5)
    }
    else if(i==4)
    {
      col<-c(pal[1],pal[2],pal[3],pal[4])#c('lightpink1','lightgreen','lightgoldenrod1','lightskyblue1')
      circ<-c(pal[1],pal[2],pal[3],pal[4])
      alpha<-c(0.5,0.5,0.5,0.5)
    }
    print(col)
    
    output$venn_diag <- renderPlot({
      venn.plot <- venn.diagram(venn_list, NULL, fill=col,alpha=alpha,col=circ,
                                cex = 2, cat.fontface=4, category.names=venn_names)
      grid.draw(venn.plot)                     
      
    })
    #download venn diagram
    output$download_venny <- downloadHandler(
      filename = paste(" Venn diagram of DE genes: ",nrow(selected),"comparisons",'.pdf'),
      content = function(file) {
        #png(file,width=800, height=500)
        pdf(file)
        venn.plot <-venn.diagram(venn_list, NULL, fill=col,alpha=alpha,col=circ,
                                 cex = 2, cat.fontface=4, category.names=venn_names)
        grid.draw(venn.plot) 
        dev.off()
        
      }
      
    )
  }
})
}