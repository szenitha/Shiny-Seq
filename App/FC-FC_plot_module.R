

FC_FC_plot_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    DT::dataTableOutput(ns("de_fc")),
    downloadButton(ns('download_fc'), 'Download plot'),
    #fluidRow(column(8,
    div( style = "position:relative",
         plotOutput(ns("fc_plot"), click = ns("plot_click_1"),
                    brush = ns("plot_brush"),
                    hover = hoverOpts(ns("plot_hover_1"), delay = 100, delayType = "debounce")),
         uiOutput(ns("hover_info"))),
    verbatimTextOutput(ns("info_0")),
    p(h5("Point clicked")),
    verbatimTextOutput(ns("info_1")),
    p(h5("Points selected")),
    verbatimTextOutput(ns("info_2"))
  )
}


FC_FC_plot_module<-function(input,output,session,DE_genes,
                            combination)
{
  result<-DE_genes()
  combo<-combination()
  num <- length(combo())
  ########fc fc plot#########
  #display an interactve table that summarizes the DE genes identified for each comparisons
  output$de_fc <- DT::renderDataTable({
    
    res<-data.frame(matrix(NA, nrow = length(combo()), ncol = 3))
    rownames(res)<-lapply(1:length(combo()), function(i) {
      combo()[[i]]
      #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
      
    })
    colnames(res)<-c('up-regulated_genes','down-regulated_genes','All-DE')
    print(res)
    print(res[1,1])
    for(i in 1:length(combo()))
    {
      print(nrow(as.data.frame(result()[[i]][1])))
      res[i,1]<-nrow(as.data.frame(result()[[i]][1]))
      res[i,2]<-nrow(as.data.frame(result()[[i]][2]))
      res[i,3]<-nrow(as.data.frame(result()[[i]][3]))
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
  
  #Displays the FC-FC plot for the comparisons selected
  observeEvent(input$de_fc_cell_clicked,{
    print('hey')
    print(input$de_fc_cells_selected)
    print(input$de_fc_cell_clicked)
    selected <- input$de_fc_cells_selected
    print(nrow(selected))
    venn_list<-list()
    
    #set order of columns in expression data as same as order of sample  ID in pheno data
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Processing Data", value = 0)
    
    comp<-lapply(1:length(combo()), function(i) {
      combo()[[i]]
      #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
      
    })
    if(nrow(selected) == 2)
    {
      row1<-selected[1,1]
      print(row1)
      col1<-selected[1,2]
      row2<-selected[2,1]
      col2<-selected[2,2]
      gene<-list('up reg genes','down reg genes','All DE')
      xlab<-paste(comp[row1],gene[col1])
      print(xlab)
      ylab<-paste(comp[row2],gene[col2])
      l<-list()
      A<-as.data.frame(result()[[row1]][4])
      #print(head)
      #genes<-A[,1]
      #A<-A[-1]
      #rownames(A)<-genes
      print(head(A))
      B<-as.data.frame(result()[[row2]][4])#input[[paste0("p-value",i)]]
      #genes<-B[,1]
      #B<-B[-1]
      #rownames(B)<-genes
      print(head(B))
      AB_merge = merge(A,B,by=0,all=T)
      print(head(AB_merge))
      # AB_merge<-transform(Reduce(merge, lapply(l, function(x) data.frame(rn = rownames(x),x))), row.names=rn, rn=NULL)
      # print(head(AB_merge))
      rownames(AB_merge)<-AB_merge$Row.names
      AB_merge<-AB_merge[,-1]
      # #both
      A.idx<-NULL
      B.idx<-NULL
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/5, detail = paste("Doing part", 1,"/",5))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      if(col1==1)
      {
        A.idx<-which(AB_merge$padj.x<input[[paste0("p-value",row1)]] & AB_merge$log2FoldChange.x>0)
      }
      else if(col2==2)A.idx<-which(AB_merge$padj.x<input[[paste0("p-value",row1)]] & AB_merge$log2FoldChange.x<0)
      else A.idx<-which(AB_merge$padj.x<input[[paste0("p-value",row1)]])
      print(length(A.idx))
      print(input[[paste0("p-value",row2)]])
      if(col2==1)
      {
        B.idx<-which(AB_merge$padj.y<input[[paste0("p-value",row2)]] & AB_merge$log2FoldChange.y>0)#which(rownames(AB_merge) %in% rownames(B))
      }
      else if(col2==2) B.idx<-which(AB_merge$padj.y<input[[paste0("p-value",row2)]] & AB_merge$log2FoldChange.y<0)
      else B.idx<-which(AB_merge$padj.y<input[[paste0("p-value",row2)]])
      
      print(length(B.idx))
      DE_both<-intersect(A.idx,B.idx)
      print(length(DE_both))
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/5, detail = paste("Doing part", 2,"/",5))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      AB_merge[       ,"DE_genes"] = "not_DE"
      AB_merge[A.idx,"DE_genes"] = paste('only',xlab)#partek
      #print(head(AB_merge))
      
      AB_merge[B.idx,"DE_genes"] = paste('only',ylab)#DE_B"
      
      AB_merge[DE_both  ,"DE_genes"] = "Both"
      
      AB_merge[       ,"DE_genes"] = as.factor(AB_merge$DE_genes)
      
      AB_merge$DE_genes = factor(AB_merge$DE_genes ,levels = c("not_DE",paste('only',xlab),paste('only',ylab),"Both"))
      print(head(AB_merge))
      print('hey')
      #print(head(AB_merge))
      AB_merge_df<-DF <- data.frame(matrix(unlist(AB_merge), nrow=nrow(AB_merge)),stringsAsFactors=FALSE)
      Merged_FCs  = AB_merge[complete.cases(AB_merge_df),]
      Merged_FCs$DE_genes = factor(Merged_FCs$DE_genes ,levels = c("not_DE",paste('only',xlab),paste('only',ylab),"Both"))
      # print(subset(Merged_FCs,Merged_FCs$DE=="DE_A"))
      # print(subset(Merged_FCs,Merged_FCs$DE=="DE_B"))
      # print(subset(Merged_FCs,Merged_FCs$DE=="DE_both"))
      p<-subset(Merged_FCs,Merged_FCs$DE_genes!="not_DE")
      print(head(p))
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/5, detail = paste("Doing part", 3,"/",5))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      lim<-max(max(abs(p$log2FoldChange.x)),max(abs(p$log2FoldChange.y)))
      #print(abs(log2(input[[paste0("FC_cutoff",row1)]])))
      #print(log2(input[[paste0("FC_cutoff",row2)]]))
      #p<-subset(AB_merge,AB_merge$DE=="DE_B")
      fc_row<-input[[paste0("FC_cutoff",row1)]]
      fc_col<-input[[paste0("FC_cutoff",row2)]]
      if(input[[paste0("hypothesis_choice", row1)]]==1) fc_row <-1
      # fc_row<-input[[paste0("FC_cutoff",row1)]]
      #else fc_row <-1
      if(input[[paste0("hypothesis_choice", row2)]]==1) fc_col <-1#fc_row<-input[[paste0("FC_cutoff",row2)]]
      #else fc_col <-1
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/5, detail = paste("Doing part", 4,"/",5))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      output$fc_plot <- renderPlot({
        #plot(p$log2FoldChange.x,p$log2FoldChange.y)
        ggplot(p, aes(x = log2FoldChange.x, y = log2FoldChange.y,color = DE_genes,shape=DE_genes)) +
          geom_point(alpha=0.5) +
          #   # scale_color_manual(values = c("Green", "Red", "Black")) +
          #   # scale_shape_manual(values = c(20,20,20)) +
          geom_vline(xintercept = -log2(fc_row))+#input[[paste0("FC_cutoff",row1)]]))+
          geom_vline(xintercept =log2(fc_row))+#input[[paste0("FC_cutoff",row1)]])) +
          geom_hline(yintercept = -log2(fc_col))+#input[[paste0("FC_cutoff",row2)]]))+
          geom_hline(yintercept = log2(fc_col))+#input[[paste0("FC_cutoff",row2)]])) +
          geom_abline(slope = 1, intercept = 0) +
          geom_smooth(method = lm,aes(color = NULL, shape = NULL), se = F)+
          scale_y_continuous(limits = c(-lim,lim))+
          scale_x_continuous(limits = c(-lim,lim))+
          #   # # scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000), limits = c(0.0001,1000),labels = comma) + 
          #   #scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000), limits = c(0.0001,1000),labels = comma) +
          theme_bw() +
          xlab(xlab) + ylab(ylab)+
          ggtitle("FC-FC-plot")
        
        
        # venn.plot <- venn.diagram(venn_list, NULL, fill=col,alpha=alpha,col=circ,
        #                           cex = 2, cat.fontface=4, category.names=venn_names)
        # grid.draw(venn.plot)                     
        
      })
      
      # Increment the progress bar, and update the detail text.
      progress$inc(1/5, detail = paste("Doing part", 5,"/",5))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      #print(brushedPoints(p, input$plot_brush))
      output$info_0 <- renderText({
        if(!is.null(p))
        {
          paste0('x corresponds to ',xlab,"\n",
                 'y corresponds to ',ylab)
        }
      })
      #dispyas the gene selected in the plot
      output$info_1 <- renderPrint({
        if(!is.null(p))
        {
          nearPoints(p[,-c(4,5,11,12,15)], input$plot_click_1, threshold = 1, maxpoints = 1,addDist = FALSE)
        }
      })
      #dispyas the genes selected in the plot
      output$info_2 <- renderPrint({
        if(!is.null(p))
        {
          brushedPoints(p[,-c(4,5,11,12,15)], input$plot_brush)
        }
      })
      #dispyas the gene on hover in the plot
      output$hover_info <- renderUI({
        if(!is.null(p))
        {
          hover <- input$plot_hover_1
          point <- nearPoints(p, hover, threshold = 5, maxpoints = 1, addDist = FALSE)
          if (nrow(point) == 0) return(NULL) 
          # calculate point position INSIDE the image as percent of total dimensions
          # from left (horizontal) and from top (vertical)
          left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
          top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
          # calculate distance from left and bottom side of the picture in pixels
          left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
          top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
          # create style property fot tooltip 
          # background color is set so tooltip is a bit transparent
          # z-index is set so we are sure are tooltip will be on top
          style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                          "left:", left_px + 2, "px; top:", top_px + 2, "px;") 
          # actual tooltip created as wellPanel 
          wellPanel(       style = style,       p(HTML(paste0("<b> Gene name: </b>", rownames(point), "<br/>",
                                                              "<b> foldchange x: </b>",point$FoldChange.x,"<br/>",
                                                              "<b> log2foldchange x: </b>", point$log2FoldChange.x, "<br/>",
                                                              "<b> foldchange y: </b>",point$FoldChange.y,"<br/>",
                                                              "<b> log2foldchange y:</b>", point$log2FoldChange.y, "<br/>"       
                                                              #"<b> Distance from left: </b>", 
                                                              #left_px, "<b>, from top: </b>", top_px
          )))     )
          
        }
      })
      #Compute FC-FC plot
      fcplot<-reactive(
        {
          ggplot(p, aes(x = log2FoldChange.x, y = log2FoldChange.y,color = DE_genes,shape=DE_genes)) +
            geom_point(alpha=0.5) +
            #   # scale_color_manual(values = c("Green", "Red", "Black")) +
            #   # scale_shape_manual(values = c(20,20,20)) +
            geom_vline(xintercept = -log2(fc_row))+#input[[paste0("FC_cutoff",row1)]]))+
            geom_vline(xintercept =log2(fc_row))+#input[[paste0("FC_cutoff",row1)]])) +
            geom_hline(yintercept = -log2(fc_col))+#input[[paste0("FC_cutoff",row2)]]))+
            geom_hline(yintercept = log2(fc_col))+#input[[paste0("FC_cutoff",row2)]])) +
            geom_abline(slope = 1, intercept = 0) +
            geom_smooth(method = lm,aes(color = NULL, shape = NULL), se = F)+
            scale_y_continuous(limits = c(-lim,lim))+
            scale_x_continuous(limits = c(-lim,lim))+
            #   # # scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000), limits = c(0.0001,1000),labels = comma) + 
            #   #scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000), limits = c(0.0001,1000),labels = comma) +
            theme_bw() +
            xlab(xlab) + ylab(ylab)+
            ggtitle("FC-FC-plot")
        })
      #download FC-FC plot
      output$download_fc <- downloadHandler(
        filename = paste(" fc plot of DE genes: ",'comparisons.svg'),#xlab," ",ylab,
        content = function(file) {
          #png(file,height=500, width=800)
          ggsave(file,fcplot())#, height=5, width=5, units='in', dpi=600)
          #dev.off()
          
        }
        
      )
    }
  })
}
