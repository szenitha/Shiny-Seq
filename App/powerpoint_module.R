########################################################################################
########################################################################################
#####download ppt
source("./functions.r")
powerpoint_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    downloadButton(ns('downloadppt'), 'Download ppt'),
    downloadButton(ns('download_all_tables'), 'Download All tables')
  )
}
#combo<-combination()

powerpoint_module<-function(input,output,session,
                            normal,dds.fc,conchoice,top_unorm,top_norm,pca_input_norm,pca_input_raw,
                            top_b,batch_corrected,batch_choice,combination,ok3,DE_genes,DE_TF,
                            input_scale,p_values,input_scale_volx,input_scale_voly,hypothesis_choice,
                            input_Distance_de,input_Linkage_de,input_Distance_tf,input_Linkage_tf,
                            dds,input_Distance_anova,input_Linkage_anova,
                            Enriched_Kegg_table,Enriched_Kegg_obj,
                            Enriched_BP_table,Enriched_BP_obj,
                            Enriched_hall_table,Enriched_hall_obj)
  
{
  output$downloadppt <- downloadHandler(
    filename = "result.pptx",
    content = function(file) {
      doc = pptx( )
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Processing Data", value = 0)
      combo<-combination()
      n<-7+length(combo())
      
      # Slide 1 : Title slide
      #+++++++++++++++++++++++
      doc <- addSlide(doc, "Title Slide")
      doc <- addTitle(doc,"Data Analysis results of project")
      doc <- addSubtitle(doc, "Presented by:")
      # 
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 1,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      # Slide 2 : Add boxplot
      #+++++++++++++++++++++++
      doc <- addSlide(doc, "Comparison")
      doc <- addTitle(doc, "Box Plot")
      doc <- addParagraph(doc, "Raw data")
      
      doc <- add_slide(doc,layout = "Two Content", master = "Office Theme") %>%
        ph_with_text(type = "body", str = "Raw Data", index = 1) %>%
        ph_with_text(type = "body", str = "Normalized Data", index = 2) %>%
        ph_with_text(type = "title", str = "Box Plot")
      
      doc<-addParagraph(doc,"")
      # #boxplot_raw<-
      doc <- addPlot(doc,fun=function() boxplot_output(assay(dds.fc()),colData(dds.fc()),as.numeric(conchoice())),
                      offx = 0.5, offy = 2, width = 6, height = 5 )
      doc <- addParagraph(doc, "Normalized data")
      # #boxplot_norm<-
      doc <- addPlot(doc,fun=function() boxplot_output(normal(),colData(dds.fc()),as.numeric(conchoice())),
                     offx = 7, offy = 2, width = 6, height = 5 )
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 2,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      
      # Slide 3 : Add 2d pca raw and normalized
      #+++++++++++++++++++++++
      # Sys.setenv("plotly_username" = "szenitha")
      # Sys.setenv("plotly_api_key" = "BgkfmFBgWVxaU6rtWkFc")
      # m <- plot_ly(x = 1:10)
      # saveWidget(as.widget(raw_pca()[[2]]), "E:/Zenitha/temp.html")
      # webshot("E:/Zenitha/temp.html", file = "E:/Zenitha/test.png",
      #         cliprect = "viewport")
      ## save to html
      # htmlwidgets::saveWidget(widget = raw_pca()[[2]], file = "E:/Zenitha/temp.html",selfcontained = FALSE)
      
      ## save html to png 
      #  webshot::webshot(url = "file:///E:/Zenitha/temp.html",
      #                   file = "E:/Zenitha/test.png", delay = 10 )#
      #@param x a plotly object
      # screenshot_plotly = function(x, file = 'webshot.png', ...) {
      #   file2 = normalizePath(file, mustWork = FALSE)
      #   #d = tempfile()
      #   #dir.create(d)
      #   #owd = setwd(d)
      #   owd = getwd()
      #   on.exit({
      #     setwd(owd)
      #     #unlink(d, recursive = TRUE)
      #   }, add = TRUE)
      #   htmlwidgets::saveWidget(x, 'index.html', FALSE)
      #   file.copy(webshot::webshot('index.html', ...), file2)
      #   file
      # }
      # screenshot_plotly(raw_pca()[[2]])
      #2d raw
      raw_pca<-pcaplot(pca_input_raw(),conchoice(),top_unorm(),"2D",NULL)[[2]]
      #ggsave(filename="./un_normalized_2D_pca.png", plot=raw_pca,dpi = 600)
      #2d normalized
      norm_pca<-pcaplot(pca_input_norm(),conchoice(),top_norm(),"2D",NULL)[[2]]
      #ggsave(filename="./normalized_2D_pca.png", plot=norm_pca,dpi = 600)
      title<-""
      subtitle1<-"Raw data"
      subtitle2<-"Normalized data"
      if(top_norm()==2 && top_unorm()==2) title<-"2D PCA of All genes"
      else if(top_norm()==1 && top_unorm()==1) title<-"2D PCA of top 500 most variable genes"
      else
      {
        title<-"2D PCA"
        if(top_unorm()==2)
        {
          subtitle1<-"Raw data of all genes"
        }
        else
        {
          subtitle1<-"Raw data of top 500 most variable genes"
        }
        if(top_norm()==2)
        {
          subtitle2<-"Raw data of all genes"
        }
        else
        {
          subtitle2<-"Raw data of top 500 most variable genes"
        }
      }
      doc <- addSlide(doc, "Comparison")
      doc <- addTitle(doc, title)
      doc <- addParagraph(doc, subtitle1)
      doc <- addParagraph(doc, "")
      if( capabilities(what = "png") )
        doc <- ph_with_gg(doc, value = raw_pca )
      #doc <- addImage(doc, "./un_normalized_2D_pca.png",offx = 0.6, offy = 2.7, width = 6, height = 4)
      doc <- addParagraph(doc, subtitle2)
      if( capabilities(what = "png") )
        doc <- ph_with_gg(doc, value = norm_pca )
      #doc <- addImage(doc, "./normalized_2D_pca.png",offx = 7, offy = 2.7, width = 6, height = 4 )
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 3,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      # Slide 4 : Add 3d pca raw and normalized
      #+++++++++++++++++++++++
      #3d pca raw and normalized
      # plotly_IMAGE(raw_pca()[[2]], format = "png", out_file = "un normalized 3D pca.png")
      #3d normalized
      # plotly_IMAGE(norm_pca()[[2]], format = "png", out_file = "normalized 3D pca.png")
      title<-""
      subtitle1<-"Raw data"
      subtitle2<-"Normalized data"
      if(input$top_norm==2 && input$top_unorm==2) title<-"3D PCA of All genes"
      else if(input$top_norm==1 && input$top_unorm==1) title<-"3D PCA of top 500 most variable genes"
      else
      {
        title<-"3D PCA"
        if(input$top_unorm==2)
        {
          subtitle1<-"Raw data of all genes"
        }
        else
        {
          subtitle1<-"Raw data of top 500 most variable genes"
        }
        if(input$top_norm==2)
        {
          subtitle2<-"Raw data of all genes"
        }
        else
        {
          subtitle2<-"Raw data of top 500 most variable genes"
        }
      }

      doc <- addSlide(doc, "Comparison")
      doc <- addTitle(doc, title)
      doc <- addParagraph(doc, subtitle1)
      doc <- addParagraph(doc, "")
      #doc <- addImage(doc, "un normalized 3D pca.png",offx = 0.6, offy = 2.7, width = 6, height = 4)
      doc <- addParagraph(doc, subtitle2)
      #doc <- addImage(doc, "normalized 3D pca.png",offx = 7, offy = 2.7, width = 6, height = 4 )
      # Increment the progress bar, and update the detail text.
      progress$inc(1/(n), detail = paste("Doing part", 4,"/",(n)))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)

      # #  # Slide 5 & 6: Add batch and normalized(2D and 3D)
      #  #+++++++++++++++++++++++
      
      if(!is.null(batch_choice()) && as.numeric(batch_choice()!=1))
      {
        pca_input_batch<-list(assay(dds.fc()),colData(dds.fc()),batch_corrected())
        #print(head(pca_input_batch))
        print(top_b())
        top_pca_batch<-1
        top_batch<-top_b()

        if(!is.null(top_batch))
          {
          top_pca_batch<-top_batch()
        }

        batch_corrected_pca<-pcaplot(pca_input_batch,conchoice(),top_pca_batch,"2D",NULL)[[2]]
        ggsave(filename="./batch_corrected_2D_pca.png", plot=batch_corrected_pca,dpi = 600)
        title<-""
        subtitle1<-"Normalized data"
        subtitle2<-"Batch corrected data"
        if(top_norm()==2 && top_pca_batch==2) title<-"2D PCA of All genes"
        else if(top_norm()==1 && top_pca_batch==1) title<-"2D PCA of top 500 most variable genes"
        else
        {
          title<-"2D PCA"
          if(top_norm()==2)
          {
            subtitle1<-"Normalized data of all genes"
          }
          else
          {
            subtitle1<-"Normalized data of top 500 most variable genes"
          }
          if(top_pca_batch==2)
          {
            subtitle2<-"Batch corrected data of all genes"
          }
          else
          {
            subtitle2<-"Batch corrected data of top 500 most variable genes"
          }
        }
        #    #slide 5
        doc <- addSlide(doc, "Comparison")
        doc <- addTitle(doc, title)
        doc <- addParagraph(doc, subtitle1 )
        doc <- addParagraph(doc, "")
        doc <- addImage(doc, "./normalized_2D_pca.png",offx = 0.6, offy = 2.7, width = 6, height = 4 )
        doc <- addParagraph(doc, subtitle2)
        doc <- addImage(doc, "./batch_corrected_2D_pca.png",offx = 7, offy = 2.7, width = 6, height = 4)

        #    #slide 6
      #   title<-""
      #   subtitle1<-"Normalized data"
      #   subtitle2<-"Batch corrected data"
      #   if(top_norm()==2 && top_pca_batch==2) title<-"3D PCA of All genes"
      #   else if(top_norm()==1 && top_pca_batch==1) title<-"3D PCA of top 500 most variable genes"
      #   else
      #   {
      #     title<-"3D PCA"
      #     if(top_norm()==2)
      #     {
      #       subtitle1<-"Normalized data of all genes"
      #     }
      #     else
      #     {
      #       subtitle1<-"Normalized data of top 500 most variable genes"
      #     }
      #     if(top_pca_batch==2)
      #     {
      #       subtitle2<-"Batch corrected data of all genes"
      #     }
      #     else
      #     {
      #       subtitle2<-"Batch corrected data of top 500 most variable genes"
      #     }
      #   }
      #   #plotly_IMAGE(batch_corrected_pca()[[2]], format = "png", out_file = "batch corrected 3D pca.png")
      #   doc <- addSlide(doc, "Comparison")
      #   doc <- addTitle(doc, title)
      #   doc <- addParagraph(doc, subtitle1 )
      #   doc <- addParagraph(doc, "")
      #   #doc <- addImage(doc, "normalized 3D pca.png",offx = 0.6, offy = 2.7, width = 6, height = 4 )
      #   doc <- addParagraph(doc, subtitle2)
      #   #doc <- addImage(doc, "batch corrected 3D pca.png",offx = 7, offy = 2.7, width = 6, height = 4)
       }
      # # Increment the progress bar, and update the detail text.
      # progress$inc(1/(n), detail = paste("Doing part", 5,"/",(n)))
      # # Pause for 0.1 seconds to simulate a long computation.
      # Sys.sleep(0.1)

      if(!is.null(ok3()) && ok3()>0)
      {
        # Slide 7: Add heatmap of 1000 most variable genes
        #+++++++++++++++++++++++
        rld<-NULL
        if(!is.null(batch_choice()) && as.numeric(batch_choice()==1))
        {
          rld<-pca_input_norm()[[3]]
        }
        else rld<-batch_corrected()
        png("./Heatmap_of_1000_most_variable_genes.png",width = 980, height = 680)

        heatmap_genes("ANOVA",dds(),dds.fc(),rld,NULL,
                      DE_genes,NULL,NULL,
                      num,"Heatmap_of_1000_most_variable_genes",
                      input_Distance_anova(),input_Linkage_anova())
        dev.off()
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc, "Heatmap of 1000 most variable genes")
        doc <- addImage(doc,"./Heatmap_of_1000_most_variable_genes.png",vector.grahics=T,offx = 2, offy = 2, width = 8, height = 5 )#vector graphic to make image editable in ppt
        # Increment the progress bar, and update the detail text.
        progress$inc(1/(n), detail = paste("Doing part", 6,"/",(n)))
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        # Slide 8: summary slide
        #+++++++++++++++++++++++
        # Slide 9: comparison
        #+++++++++++++++++++++++
        combo<-combination()
        num<- length(combo())
        #enrichment kegg data
        res<-Enriched_Kegg_obj()
        result<-Enriched_Kegg_table()
        res_bp<-Enriched_BP_obj()
        result_bp<-Enriched_BP_table()
        res_hall<-Enriched_hall_obj()
        result_hall<-Enriched_hall_table()
        result_de<-DE_genes()
        result_tf<-DE_TF()

        #DE genes
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc,"summary: DE genes")
        dat<-summary_analysis(result_de,result_tf,result,result_bp,result_hall,combo())
        MyFTable <- vanilla.table(dat[[1]],add.rownames=TRUE)
        MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
        doc <- addFlexTable( doc, MyFTable)
        #TF
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc,"summary: Transcription factors")
        MyFTable <- vanilla.table(dat[[2]],add.rownames=TRUE)
        MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
        doc <- addFlexTable( doc, MyFTable)
        #Kegg
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc,"summary: Enriched Kegg pathways")
        MyFTable <- vanilla.table(dat[[3]],add.rownames=TRUE)
        MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
        doc <- addFlexTable( doc, MyFTable)
        #GO Terms
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc,"summary: Enriched GOTerms(Biological processes)")
        MyFTable <- vanilla.table(dat[[4]],add.rownames=TRUE)
        MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
        doc <- addFlexTable( doc, MyFTable)
        #Hallmark
        doc <- addSlide(doc, "Title and Content")
        doc <- addTitle(doc,"summary: Enriched Hallmark")
        MyFTable <- vanilla.table(dat[[5]],add.rownames=TRUE)
        MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
        doc <- addFlexTable( doc, MyFTable)

        # Increment the progress bar, and update the detail text.
        progress$inc(1/(n), detail = paste("Doing part", 7,"/",(n)))
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      #
        for (i in 1:num)
        {
      #     #p-value plot
           name<-str_replace_all(combo()[[i]],"[^[:alnum:]]"," ")
          doc <- addSlide(doc, "Title Only")
          doc <- addTitle(doc,paste0("comparison ",i,": ",name))
          doc <- addSlide(doc, "Title and Content")
          doc <- addTitle(doc, paste("p value plot",name))
          p<-function(){p_value_all(result_de,i,combo())}
          doc <- addPlot(doc,p)
      #     #ma plot
          doc <- addSlide(doc, "Title and Content")
          doc <- addTitle(doc, paste("MA plot",name))
          ma<-function(){ma_plot(i,combo(),result_de,input_scale(),p_values())}
          doc <- addPlot(doc,ma)

      #     #volcano plot
      #     #nam<-paste("condition",str_replace_all(combination()[[i]][2],"[^[:alnum:]]","."),sep = "")
          doc <- addSlide(doc, "Title and Content")
          doc <- addTitle(doc, paste("Volcano plot",name))
          nam<- str_replace_all(paste("./volcano plot ",name,".png",sep = "")," ","_")
          print(input_scale_volx())
          print(input_scale_voly())
          print("line 364")
          lim_y<-input_scale_voly()
          if(is.null(lim_y))
          {
            num<- length(combo())
            #Get the deseq2 dataset
            res_vol<-result_de[[as.numeric(i)]][5] [[1]]
            num<-(-log10(na.omit(res_vol$padj)))
            idx<-which(is.finite(num))
            max<-max(round(num[idx],1))
            #max<-(-log10(min(na.omit(res$padj))))
            #max<-round(max,1)
            lim_y<-max+5
          }
          lim_x<-input_scale_volx()
          if(is.null(lim_x))
          {
            num<- length(combo())
            #Get the deseq2 dataset
            res_vol<-result_de[[as.numeric(i)]][5] [[1]]
            max<-max(na.omit(abs(res_vol$log2FoldChange)))
            max<-round(max,1)
            print(max)
            lim_x<-max+5
          }
          print(lim_x)
          print(lim_y)
          p<-volcano_plot(i,combo(),result_de,
                          lim_x,lim_y,
                          p_values(),hypothesis_choice())[[1]]
          ggsave(filename=nam, plot=p,
                 device="png"

                   )
          doc <- addImage(doc,nam)
          #comparison of heatmap of de genes and heatmap of tf
          nam1<- str_replace_all(paste("./Heatmap of DE genes for ",name,".png",sep = "")," ","_")
          rld<-NULL
          print(batch_choice())
          if(!is.null(batch_choice()) && as.numeric(batch_choice()==1))
          {
            print("ppt 416")
            rld<-pca_input_norm()[[3]]
          }
          else rld<-batch_corrected()
          print("ppt 419")
          print(head(rld))
          h<-heatmap_genes("DE",NULL,dds.fc(),rld,i,
                           DE_genes,NULL,NULL,
                           num,paste("Heatmap of DE genes for ",name),
                           input_Distance_de(),input_Linkage_de())
          #print(h)
          if(is.null(h)) {
            png(nam1,width = 980, height = 680)
            dev.off()
          }
          else
          {
            png(nam1,width = 980, height = 680)
            heatmap_genes("DE",NULL,dds.fc(),rld,i,
                          DE_genes,NULL,NULL,
                          num,paste("Heatmap of DE genes for ",name),
                          input_Distance_de(),input_Linkage_de())
            dev.off()
          }

      #     #tf heatmap

          if(!is.null(result_tf))
          {
            h<-heatmap_genes("TF",NULL,dds.fc(),rld,i,
                             DE_TF,NULL,NULL,
                             num,paste("Heatmap of TF genes for ",name),
                             input_Distance_tf(),input_Linkage_tf())
            print("ppt line 434")
            #print(h)
            if(is.null(h))
            {
              # png(nam1,width = 980, height = 680)
              # dev.off()
              doc <- addSlide(doc, "Title and Content")
              doc <- addTitle(doc, "Heatmap of DE genes")
              #doc <- addParagraph(doc, "DE genes")
              #doc <- addParagraph(doc, "")
              doc <- addImage(doc,nam1,offx = 0.6, offy = 2.7, width = 6, height = 4)
            }
            else
            {
              # heatmap_genes("TF",NULL,dds.fc(),rld,i,
              #                  DE_TF,NULL,NULL,
              #                  num,paste("Heatmap of TF genes for ",name),
              #                  input_Distance_tf(),input_Linkage_tf())
              nam<- str_replace_all(paste("./Heatmap of TF for ",name,".png",sep = "")," ","_")
              png(nam,width = 980, height = 680)
              heatmap_genes("TF",NULL,dds.fc(),rld,i,
                            DE_TF,NULL,NULL,
                            num,paste("Heatmap of TF genes for ",name),
                            input_Distance_tf(),input_Linkage_tf())
              dev.off()
              doc <- addSlide(doc, "Comparison")
              doc <- addTitle(doc, "Heatmap")
              doc <- addParagraph(doc, "DE genes")
              doc <- addParagraph(doc, "")
              doc <- addImage(doc,nam1,offx = 0.6, offy = 2.7, width = 6, height = 4)
              doc <- addParagraph(doc, "Transcription factors")
              doc <- addImage(doc,nam,offx = 7, offy = 2.7, width = 6, height = 4)
            }


          }
          else
          {
            doc <- addSlide(doc, "Title and Content")
            doc <- addTitle(doc, "Heatmap of DE genes")
            doc <- addImage(doc,nam1,offx = 0.6, offy = 2.7, width = 6, height = 4)

          }
      #
          #up regulated genes

          #barplot of kegg pathways
          #barplot of go terms
          #barplot of hallmark
          df<-as.data.frame(result_de[[i]][1])
          print(head(df))
          print(head(df$FoldChange))
          genes<-df[,1]
          df<-df[-1]
          rownames(df)<-genes
          print("temp")
          print(head(order(unlist(df$FoldChange),decreasing=TRUE)))
          df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
          colnames(df_final)[1]<-"Overall mean"
          print(colnames(df_final))
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          a_tab<-as.data.frame(df_final[,-4])
          print(colnames(a_tab))
          c<-colnames(a_tab)
          print(length(c))
          temp<-as.vector(c[7:length(c)])
          print(length(temp))
          print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
          print('howdy')
          print(c)
          df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
          data<-df_de[,1:(3+length(temp))]
          print(head(data))
          ######
          if(nrow(data)>0)
          {
            doc <- addSlide(doc, "Title and Content")
            doc <- addTitle(doc,"Top 10 Up regulated genes")
            MyFTable <- vanilla.table(head(data,10),add.rownames=TRUE)
            MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
            doc <- addFlexTable( doc, MyFTable)
          }
          #up regulated tf

          df<-as.data.frame(result_tf[[i]][[1]])
          # genes<-df[,1]
          # df<-df[-1]
          # rownames(df)<-genes

          ######
          if(nrow(df)>0)
          {
            df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
            colnames(df_final)[1]<-"Overall mean"
            for (j in 8:length(colnames(df_final)))
            {
              colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
            }
            data<-as.data.frame(df_final[,-c(1,3,4,5)])
            print(head(data))

            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(data)>=10) title<-paste("Top 10 Up regulated Transcription Factors",name)
            else title<-paste("All Up regulated Transcription Factors",name)
            print(title)
            doc <- addTitle(doc,title)
            MyFTable <- vanilla.table(head(data,10),add.rownames=TRUE)
            MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
            doc <- addFlexTable( doc, MyFTable)
          }

          #up regulated kegg
          r<-result[[i]][[1]]#as.data.frame(result[[i]][[1]])
          print('kegg')
          print(typeof(r))
          print(head(r))
          if(nrow(r)>0)
          {
            #kegg<-function(){keggplot_all(i,1,res)}
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 up regulated kegg pathways",name)
            else title<-paste("All up regulated kegg pathways",name)
            print(title)
            #png(paste(title,'.png'))
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),
                   enrichment_plot("kegg",res,result,i,1,10,"Barplot",""),
                   dpi = 1000,height = 6,width = 12)
            #dev.off()
            doc<-addTitle(doc,title)
            #doc <- addPlot(doc,kegg)#
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }

          #biological process
          #up regulated
          r<-as.data.frame(result_bp[[i]][[1]])
          if(nrow(r)>0)
          {
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 up regulated biological processes",name)
            else title<-paste("./All up regulated biological processes",name)
            print(title)
            p<-function(){enrichment_plot("biological process",res_bp,result_bp,i,1,10,"Barplot","")}
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),p(),dpi = 1000,height = 6,width = 12)
            doc<-addTitle(doc,title)
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }
          #up regulated hallmark
          r<-result_hall[[i]][[1]]#as.data.frame(result[[i]][[1]])
          print('hallmark')
          print(typeof(r))
          print(head(r))
          if(nrow(r)>0)
          {
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 up regulated Hallmark gene sets",name)
            else title<-paste("All up regulated Hallmark gene sets",name)
            print(title)
            #png(paste(title,'.png'))
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),
                   enrichment_plot("hallmark",res_hall,result_hall,i,1,10,"Barplot",""),
                   dpi = 1000,height = 6,width = 12)
            #dev.off()
            doc<-addTitle(doc,title)
            #doc <- addPlot(doc,kegg)#
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }
          #Down regulated genes
          #pepare table
          df<-as.data.frame(result_de[[i]][2])
          genes<-df[,1]
          df<-df[-1]
          rownames(df)<-genes
          df_final<- df[order(unlist(df$FoldChange)),]
          colnames(df_final)[1]<-"Overall mean"
          print(colnames(df_final))
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          a_tab<-as.data.frame(df_final[,-4])
          print(colnames(a_tab))
          c<-colnames(a_tab)
          print(length(c))
          temp<-as.vector(c[7:length(c)])
          print(length(temp))
          print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
          print('howdy')
          print(c)
          df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
          data<-df_de[,1:(3+length(temp))]
          print(head(data))
          ######
          #down regulated genes
          if(nrow(data)>0)
          {
            doc <- addSlide(doc, "Title and Content")
            doc <- addTitle(doc,"Top 10 Down regulated genes")
            MyFTable <- vanilla.table(head(data,10),add.rownames=TRUE)
            MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
            doc <- addFlexTable( doc, MyFTable)
          }
          #down regulated tf

          df<-as.data.frame(result_tf[[i]][[2]])
          # genes<-df[,1]
          # df<-df[-1]
          # rownames(df)<-genes

          ######
          if(nrow(df)>0)
          {
            df_final<- df[order(unlist(df$FoldChange)),]
            colnames(df_final)[1]<-"Overall mean"
            for (j in 8:length(colnames(df_final)))
            {
              colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
            }
            data<-as.data.frame(df_final[,-c(1,3,4,5)])
            print("ppt line 642")
            print(head(data))
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(data)>=10) title<-paste("Top 10 Down regulated Transcription Factors",name)
            else title<-paste("All Down regulated Transcription Factors",name)
            print(title)
            doc <- addTitle(doc,title)
            MyFTable <- vanilla.table(head(data,10),add.rownames=TRUE)
            MyFTable <- setZebraStyle(MyFTable, odd = 'lightblue', even = 'white')
            doc <- addFlexTable( doc, MyFTable)
          }
          #down regulated kegg
          r<-as.data.frame(result[[i]][[2]])
          #up reg col
          r2<-as.data.frame(result[[i]][[1]])
          if(nrow(r)>0)
          {
            #kegg<-function(){keggplot_all(i,2,res)}
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 down regulated kegg pathways",name)
            else title<-paste("All down regulated kegg pathways",name)
            print(title)
            #png(paste(title,'.png'))
            num<-2
            if(nrow(r2)==0) num<-1
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),
                   enrichment_plot("kegg",res,result,i,num,10,"Barplot",""),
                   dpi = 1000,height = 6,width = 12)
            #dev.off()
            doc<-addTitle(doc,title)
            #doc <- addPlot(doc,kegg)#
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }
          #down regulated bp
          r<-as.data.frame(result_bp[[i]][[2]])
          #up reg col
          r2<-as.data.frame(result_bp[[i]][[1]])
          if(nrow(r)>0)
          {
            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 down regulated biological processes",name)
            else title<-paste("All down regulated biological processes",name)
            print(title)
            num<-2
            if(nrow(r2)==0) num<-1
            p<-function(){enrichment_plot("biological process",res_bp,result_bp,i,num,10,"Barplot","")}
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),p(),dpi = 1000,height = 6,width = 12)
            doc<-addTitle(doc,title)
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }
          #down regulated hall
          r<-as.data.frame(result_hall[[i]][[2]])
          #up reg col
          r2<-as.data.frame(result_hall[[i]][[1]])
          if(nrow(r)>0)
          {
            #kegg<-function(){keggplot_all(i,2,res)}

            doc <- addSlide(doc, "Title and Content")
            title<-NULL
            if(nrow(r)>=10) title<-paste("Top 10 down regulated Hallmark gene sets",name)
            else title<-paste("All down regulated Hallmark gene sets",name)
            print(title)
            #png(paste(title,'.png'))
            num<-2
            if(nrow(r2)==0) num<-1
            nam<-str_replace_all(paste('./',title,sep = "")," ","_")
            ggsave(paste(nam,'.png',sep = ""),
                   enrichment_plot("hallmark",res_hall,result_hall,i,num,10,"Barplot",NULL),
                   dpi = 1000,height = 6,width = 12)
            #dev.off()
            doc<-addTitle(doc,title)
            #doc <- addPlot(doc,kegg)#
            doc<-addImage(doc,paste(nam,'.png',sep = ""))
          }
          # Increment the progress bar, and update the detail text.
          progress$inc(1/(n), detail = paste("Doing part", 7+i,"/",(n)))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
        }
      }

      writeDoc(doc,file)
      
    }
  )
###########
 
########  
  #download all tables
  #incase this throws error
  #download and install rtools and set system variable 'PATH' in control panel
  #http://stackoverflow.com/questions/29129681/create-zip-file-error-running-command-had-status-127
  output$download_all_tables <- downloadHandler(
    filename = 'all_tables.zip',
    content = function(fname) {
      if(ok3()>0)
      {
        #tmpdir <- tempdir()
        setwd(getwd())#tmpdir)
        #print(tempdir())
        fs<-NULL
        ########################processing data for download#######################
        fs<-c(fs,'normalized_table.xlsx')
        #normalized table
        #write.csv(normal(),'normalized_table.csv')
        write.xlsx2(normal(), file='normalized_table.xlsx', sheetName = "Normalized table",
                    col.names = TRUE, row.names = TRUE, append = FALSE)
        #batch corrected table
        if(!is.null(batch_choice()) && as.numeric(batch_choice()!=1))
        {
          fs<-c(fs,'batch_corrected.xlsx')
          #write.csv(batch_corrected()[[1]], 'batch_corrected.csv')
          write.xlsx2(batch_corrected()[[1]], file='batch_corrected.xlsx', sheetName = "Batch corrected table",
                      col.names = TRUE, row.names = TRUE, append = FALSE)
        }
        
        #total number of comparisons
        combo<-combination()
        num<- length(combo())
        #de genes table
        result_de<-DE_genes()
        #tf table
        result_tf<-DE_TF()
        #kegg pathways
        res<-Enriched_Kegg_obj()
        result<-Enriched_Kegg_table()
        #go terms(biological process)
        res_bp<-Enriched_BP_obj()
        result_bp<-Enriched_BP_table()
        #hallmark
        res_hall<-Enriched_hall_obj()
        result_hall<-Enriched_hall_table()
        
        
        for (i in 1:num)
        {
          #up regulated 
          #de genes
          df<-as.data.frame(result_de[[i]][1])
          genes<-df[,1]
          df<-df[-1]
          rownames(df)<-genes
          df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
          colnames(df_final)[1]<-"Overall mean"
          print(colnames(df_final))
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          a_tab<-as.data.frame(df_final[,-4])
          print(colnames(a_tab))
          c<-colnames(a_tab)
          print(length(c))
          temp<-as.vector(c[7:length(c)])
          print(length(temp))
          print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
          print('howdy')
          print(c)
          df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
          data<-df_de[,1:(3+length(temp))]
          print(head(data))
          ######
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            file<-paste(condition,".xlsx")
            sheet_name<-str_replace_all(paste('Up regulated genes for ',condition)," ","_")
            fs<-c(fs,file)
            #write.csv(data, file)
            # Write the first data set in a new workbook
            write.xlsx(data, file = paste(condition,".xlsx"),
                       sheetName = sheet_name, append = FALSE)
          }
          #up regulated tf
          
          df<-as.data.frame(result_tf[[i]][[1]])
          # genes<-df[,1]
          # df<-df[-1]
          # rownames(df)<-genes
          df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
          colnames(df_final)[1]<-"Overall mean"
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          data<-as.data.frame(df_final[,-c(1,3,4,5)])
          print(head(data))
          ######
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Up regulated TF for ',condition)," ","_")
            #fs<-c(fs,file)
            #write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          
          #up regulated kegg
          data<-result[[i]][[1]]#as.data.frame(result[[i]][[1]])
          print('kegg')
          print(typeof(data))
          print(head(data))
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Up regulated kegg pathways for ',condition)," ","_")
            #fs<-c(fs,file)
            #write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a third data in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          
          #biological process
          #up regulated
          data<-as.data.frame(result_bp[[i]][[1]])
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Up regulated BP for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          #up regulated hallmark
          data<-result_hall[[i]][[1]]#as.data.frame(result[[i]][[1]])
          print('hallmark')
          print(typeof(data))
          print(head(data))
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Up regulated Hallmark gene sets for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a fourth data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          #Down regulated
          #de genes
          df<-as.data.frame(result_de[[i]][2])
          genes<-df[,1]
          df<-df[-1]
          rownames(df)<-genes
          df_final<- df[order(unlist(df$FoldChange),decreasing=TRUE),]
          colnames(df_final)[1]<-"Overall mean"
          print(colnames(df_final))
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          a_tab<-as.data.frame(df_final[,-4])
          print(colnames(a_tab))
          c<-colnames(a_tab)
          print(length(c))
          temp<-as.vector(c[7:length(c)])
          print(length(temp))
          print(c(temp,c[2],c[3],c[5],c[6],c[1],c[4]))
          print('howdy')
          print(c)
          df_de <- a_tab[,c(c[2],c[5],c[6],temp,c[3],c[1],c[4])]
          data<-df_de[,1:(3+length(temp))]
          print(head(data))
          ######
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Down regulated genes for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          #Down regulated tf
          
          df<-as.data.frame(result_tf[[i]][[2]])
          # genes<-df[,1]
          # df<-df[-1]
          # rownames(df)<-genes
          df_final<- df[order(unlist(df$FoldChange),decreasing = TRUE),]
          colnames(df_final)[1]<-"Overall mean"
          for (j in 8:length(colnames(df_final)))
          {
            colnames(df_final)[j]<-paste(colnames(df_final)[j],"mean")
          }
          data<-as.data.frame(df_final[,-c(1,3,4,5)])
          print(head(data))
          ######
          if(nrow(data)>0)
          {
            condition<-str_replace_all(str_replace_all(combo()[[i]],"[^[:alnum:]]",".")," ","_")
            sheet_name<-str_replace_all(paste('Down regulated TF for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          
          #Down regulated kegg
          data<-result[[i]][[2]]#as.data.frame(result[[i]][[1]])
          print('kegg')
          print(typeof(data))
          print(head(data))
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Down regulated kegg pathways for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          
          #biological process
          #Down regulated
          data<-as.data.frame(result_bp[[i]][[2]])
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Down regulated BP for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          #up regulated hallmark
          data<-result_hall[[i]][[2]]#as.data.frame(result[[i]][[1]])
          print('hallmark')
          print(typeof(data))
          print(head(data))
          if(nrow(data)>0)
          {
            condition<-str_replace_all(combo()[[i]],"[^[:alnum:]]",".")
            sheet_name<-str_replace_all(paste('Down regulated Hallmark gene sets for ',condition)," ","_")
            # fs<-c(fs,file)
            # write.csv(data, file)
            file<-paste(condition,".xlsx")
            # Add a second data set in a new worksheet
            write.xlsx(data, file = file, 
                       sheetName=sheet_name, append=TRUE)
          }
          
        }
        #anova table
        #preparing anova table for download
        anova<-anova_table()
        a_tab<-anova_table()[,-c(2,3)]
        ########change#######
        cond<-unique(colData(dds.fc()[[1]])[,2])
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
        
        fs<-c(fs,'ANOVA_table.xlsx')
        
        # Add a second data set in a new worksheet
        write.xlsx(anova, file = 'ANOVA_table.xlsx', 
                   sheetName="ANOVA table", append=FALSE)
        ###########################################################################
        print (fs)
        
        zip(zipfile=fname, files=fs)
      }
    },
    contentType = "application/zip"
  )    
}
