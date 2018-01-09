##When the remove samples button is clicked in the normalized tab
#####Display interactive boxplot
boxplot_output<-function(edata,pData,col)
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  data<-edata
  pheno<-pData
  print(head(data))
  # #order pheno data by condition
  # pheno<-pheno[order(pheno[,2]),]
  
  #Rename sample id with row just for display
  i<-1:nrow(pheno)
  print("pheno")
  print(pheno)
  
  #get which column u want to plot by
  var<-as.numeric(col)
  print(var)
  print(as.vector(pheno[,1]))
  # #set order of columns in expression data as same as order of sample  ID in pheno data
  data=data[,as.vector(pheno[,1])]#1
  colnames(data)<- as.vector(pheno[,var])
  #print(var)
  # boxplot.matrix(as.matrix(data),outline=FALSE,xlab='Rows',ylab='Value',col=colors[pheno[,var]],boxwex=0.25,names=lapply(i,function(x) paste('row',x[1])))
  # legend("topright", legend = unique(pheno[,var]),xpd=TRUE, pch = 16, col = colors[unique(pheno[,var])], cex=0.85,inset=0.0005)
  print(head(melt(data)))
  ggplot(data=melt(data), aes(as.factor(Var2), value,fill=as.factor(Var2))) + geom_boxplot(outlier.shape = NA)+
    labs(x="Samples",y="Gene Counts") 
}
   
pca_components<-function(dataset,column_no,annotation,resize_factor=NULL){
  #obtain PC (principle components)
  pca = prcomp(t(dataset),scale=F)
  pcaVars2=signif(((pca$sdev))/(sum((pca$sdev))),3)*100
  signed = ifelse(max(pca$x[,2] > 70), 1, -1) # make same sign
  signed1 = ifelse(max(pca$x[,3] > 70), 1, -1) # make same sign
  total_variance = sum(pcaVars2[1:3])
  c<-annotation[,column_no]
  #print(c)
  #print(typeof(c))
  return(list(pca$x[,1],signed*pca$x[,2],signed1*pca$x[,3],pcaVars2[1],pcaVars2[2],pcaVars2[3],c))
  #o/p returned is a list of (PC1,PC2,PC3,Percentage of variance contributed by PC1,
  #                                       Percentage of variance contributed by PC2,
  #                                       Percentage of variance contributed by PC3,
  #                                       annotation)
}

pcaplot<-function(data,annotation,top,plotType,resize_factor=NULL)
{
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Making plot", value = 0)

  # Number of times we'll go through the loop to update progres bar
  n1 <- 3
  #========user input=======#
  pheno<-data[[2]] #annotation table
  dat<-data[[3]] #rlog transformed data
  var<-as.numeric(annotation)#Get which column u want to plot by
  #========================#
  # print("input")
  # print(head(pheno))
  # print(head(dat))
  #compute variance of all genes. this is done by taking the row variance (rowvars function) for each gene.
  rv <- rowVars(dat)
  #pick the top 500 most variable genes. This list will be used if user wants PCA of top 500 genes to be displayed
  genes<-order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]
 
  # Increment the progress bar once input recieved, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 1,"/",n1))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Get the components of PCA
  parameters<-list()
  if(top==2) parameters<-pca_components(dat,var,pheno)# pca components computed for all genes #input
  else parameters<-pca_components(dat[genes,],var,pheno)#pca compnents of top 500 most variable genes 
  
  # print("line 84")
  # print(parameters)
  # Increment the progress bar after obtaining the components, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 2,"/",n1))
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  
  #Define color pallete. This is needed to color the sample points on the pca plot.
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  #Obtain the percentage of variance explained by each component
  xlab = paste("PC1:",parameters[[4]],"%")
  ylab = paste("PC2:",parameters[[5]],"%")
  zlab = paste("PC3:",parameters[[6]],"%")
  #condense the principle component vectors(PCA1,PCA2,PCA3) and colors assigned to each sample into a dataframe
  #df is matrix with rows correspond to samples and columns (PCA1,PCA2,PCA3,Colors assigned to each sample)
  df<-data.frame(PCA1=parameters[[1]],PCA2=parameters[[2]],PCA3=parameters[[3]],an=parameters[[7]])
  print(colnames(pheno)[var])
  print(dim(df))
  print(length(colors[unique(df$an)]))
  p<-plotly_empty() %>% layout(autosize = F,width=1000,height=800)
  if(!is.null(resize_factor)) p<-plotly_empty() %>% layout(autosize = F,width=800,height=400)
  #check if the color column of df is not empty
  print("line 105 -functions")
  p0<-NULL
  if(length(colors[unique(df$an)])!=0)
  {
    print("hey")
    if (identical(plotType, "2D")) { #variable plotType is input obtained from user: 2d/3d
      print(colors[unique(df$an)])
      p0 <- ggplot(df, aes(x = PCA1, y = PCA2,color=df$an)) + 
        geom_point()+
        labs(x=xlab,y=ylab)+
        scale_color_manual(name=colnames(pheno)[var],values=colors[unique(df$an)])+ theme_bw() #2d PCA plot
      #req(p0)
      if(!is.null(resize_factor)) p<-ggplotly(p0) %>% layout(autosize = F,dragmode ="select",width=800,height=400)
      else p<-ggplotly(p0) %>% layout(autosize = F,width=1000,height=800,dragmode = "select")
    } else {
      #print()
      if(!is.null(resize_factor))
      {
        p<- plot_ly(x=df$PCA1,y=df$PCA2,z=df$PCA3) %>%    #3D PCA plot
          add_markers(type = "scatter3d",color=df$an,colors = unique(colors[df$an])) %>%
          layout(autosize = F,width=800,height=400,scene = list(xaxis = list(title = xlab),
                                                                 yaxis = list(title = ylab),
                                                                 zaxis = list(title = zlab)))
      }
      else
      {
        p<- plot_ly(x=df$PCA1,y=df$PCA2,z=df$PCA3) %>%    #3D PCA plot
          add_markers(type = "scatter3d",color=df$an,colors = unique(colors[df$an])) %>%
          layout(autosize = F,width=1000,height=800,scene = list(xaxis = list(title = xlab),
                                                                 yaxis = list(title = ylab),
                                                                 zaxis = list(title = zlab)))
      }
      
    }
  }
  else
    {
      if(!is.null(resize_factor)) p<-plotly_empty() %>% layout(autosize = F,width=800,height=400)
      else p<-plotly_empty() %>% layout(autosize = F,width=1000,height=800)
    }

  # Increment the progress bar when PCA plot has been successfully generated, and update the detail text.
  progress$inc(1/n1, detail = paste("Doing part", 3,"/",n1))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  return(list(p,p0)) #return PCA plot
}

# The following function ensures each treatment/condition has a min of 3 samples 
# else returns treatment/condition groups containing either one or two samples
min_samples_three <- function(dds.norm) {
  print('Inside function min_samples_three')
  dt1<-as.vector(sapply( levels(dds.norm$condition), function(lvl)
  {
    #Obtain number of samples present in each treatmen/conditiont group
    num<-ncol(counts(dds.norm,normalized=TRUE)[,dds.norm$condition == lvl] )
    if(is.null(num)) num<-1
    num
    
  }))
  #dt1 is a dataframe with two columns namely treatment/condition group and number of samples in a treatment group
  print(dt1)
  #Which treatment/condition groups in dt1 have 2 samples (Obtain index)
  two_samples<-which(dt1 %in% 2)
  print('two')
  #Which treatment/condition groups in dt1 has 1 sample (Obtain index)
  one_sample<-which(dt1 %in% 1)
  print('one sample')
  
  #Get a list of all treatment/condition groups
  conditions<-levels(dds.norm$condition)
  
  #return only those treatment/condition groups having only one or two samples
  if(length(one_sample)!=0) return (list(conditions[one_sample],1))#return treatment/condition groups containing only one sample
  else if(length(two_samples)!=0) return (list(conditions[two_samples],2))#return treatment/condition groups containing only two samples
  else if(length(one_sample)!=0 && length(two_samples)!=0) NULL #return null if there are no treatment/condition groups with neither one or two samples
  
  
}

#The following plots (density and ECDF plots) is used as a quality control for normalized data
#It is triggered by clicking the QC button in normalized table tab(represented as a sub tab under normalization tab)
#Plot 1: Density plot function: returns density plot
density_plot<-function(normal,x,y,zoom) #input(normalized table,x and y limits of plots if parameter zoom =True
                                        # Zoom parameter is true when user wants to zoom in,this results in resizing the plot )
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #print('ya')
  library(geneplotter)
  #To assess whether the normalization has worked,
  #we plot the densities of counts for the different samples.
  #Since most of the genes are (heavily) affected by the experimental conditions,
  #a succesful normalization will lead to overlapping densities.
  if(zoom==F) multidensity( normal,xlab="mean counts",xlim=c(0,1000),col=colors)
  else multidensity( normal,xlab="mean counts",xlim=x,ylim=y,col=colors)
  
}
#plot 2: ECDF plot function
ECDF_plot<-function(normal)
{
  library(RColorBrewer)
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  #print('ya')
  library(geneplotter)
  #In an ECDF plot, the estimated probility is plotted on the y-axis and the count values
  #on the x-axis. I.e. you can read of the median and other quantiles from this plot.
  #As already mentioned, if the normalization has worked, the ECDFs of the different samples
  #should be overlapping.
  multiecdf(normal, normalized = T,xlab="Normalized gene counts", xlim=c(0, 1000),col=colors)
  
  }
#batch correction function
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
heatmap_sample<-function(sampleDists,vsd)
{
  library("RColorBrewer")
  library("pheatmap")
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Making plot", value = 0)
  
  # Number of times we'll go through the loop
  n <- 3
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, sep="-")
  colnames(sampleDistMatrix) <- paste(paste(vsd$condition, sep="-"),"(",colnames(vsd),")",sep=" ")#NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  # Increment the progress bar, and update the detail text.
  progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
  
  # Pause for 0.1 seconds to simulate a long computation.
  Sys.sleep(0.1)
 p<- pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
 # Increment the progress bar, and update the detail text.
 progress$inc(1/n, detail = paste("Doing part", 3,"/",n))
 
 # Pause for 0.1 seconds to simulate a long computation.
 Sys.sleep(0.1)
 p
}

#ma plot
ma_plot<-function(ma_choice,combination,DE_genes,scale,p_values)
{
  if(!is.null(ma_choice))
  {
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Processing Data", value = 0)
    
    num<- length(combination)
    
    #Get the deseq2 dataset
    #dds.fc<-batch_design()[[1]]
    res<-DE_genes[[as.numeric(ma_choice)]][5] [[1]]
    print(head(res))
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 1,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    #filter most significant genes with fdr cut off 0.05
    resSig <- subset(res, padj < p_values[[as.numeric(ma_choice)]])# & abs(res$log2FoldChange)>1)
    print(head(resSig))
    print("hey")
    #get down regulated genes
    head(resSig[ order(resSig$log2FoldChange), ])
    d<-resSig[resSig$log2FoldChange<1,]
    topGene_down<-rownames(head(d[order(d$log2FoldChange),],2))
    print(topGene_down)
    #get up regulated genes
    head(resSig[ order(resSig$log2FoldChange,decreasing = TRUE), ])
    u<-resSig[resSig$log2FoldChange>1,]
    topGene_up<-rownames(head(u[order(u$log2FoldChange,decreasing = TRUE),],2))
    print(topGene_up)
    #top variable genes
    #plotMA(res, ylim=c(-10,10),xlim=c(0.01,100))
    topGene <- rownames(resSig)[which.min(resSig$padj)]# only most significant gene
    #topGene<-rownames(res)[which(res$padj<0.05)]#all top DEG between conditions in contrast
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 2,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    print('hey')
    print(scale)
    print('ok')
    plotMA(res,ylim=c(-1*scale,scale))
    # with(res[topGene, ], {
    #     points(baseMean, log2FoldChange, col="green", cex=2, lwd=2)
    #     text(baseMean, log2FoldChange, topGene, pos=2, col="green")
    #   })
    #   #display both up and down regulated genes
    #   with(res[topGene_down, ], {
    #     points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
    #     text(baseMean, log2FoldChange, topGene_down, pos=2, col="dodgerblue")
    #   })
    #   with(res[topGene_up, ], {
    #     points(baseMean, log2FoldChange, col="hot pink", cex=2, lwd=2)
    #     text(baseMean, log2FoldChange, topGene_up, pos=2, col="hot pink")
    #   })
      # Increment the progress bar, and update the detail text.
      progress$inc(1/3, detail = paste("Doing part", 3,"/",3))
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
  }
}

#volcano plot
volcano_plot<-function(vol_choice,combination,DE_genes,scale_volx,scale_voly,p_values,hypothesis_choice)
{
  # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    num<- length(combination)
    #Get the deseq2 dataset
    res<-DE_genes[[as.numeric(vol_choice)]][5] [[1]]
    print(head(res))
    
    progress$set(message = "Processing Data", value = 0)
    
    padj_cutoff<-p_values[[as.numeric(vol_choice)]]
    
    logfc_cutoff<-0
    if(hypothesis_choice[[as.numeric(vol_choice)]]!=1)
    {
      logfc_cutoff<-log2(input[[paste0("FC_cutoff",as.numeric(vol_choice))]])
    }
    
    gene_list<-as.data.frame(res[c("log2FoldChange", "padj")])
    ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    gene_list$group = as.factor(abs(gene_list$log2FoldChange) > logfc_cutoff & gene_list$padj < padj_cutoff)
    #gene_list<-complete.cases(gene_list)
    print(head(gene_list))
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 1,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    # # Find and label the top peaks..
    gene_list_complete<-na.omit(gene_list)
    top_peaks <- gene_list_complete[with(gene_list_complete, order(log2FoldChange,padj)),][1:10,]
    top_peaks <- rbind(top_peaks, gene_list_complete[with(gene_list_complete, order(-log2FoldChange,padj)),][1:10,])
    # Add gene labels for all of the top genes we found
    # here we are creating an empty list, and filling it with entries for each row in the dataframe
    # each list entry is another list with named items that will be used by Plot.ly
    print(top_peaks)
    a <- list()
    for (i in seq_len(nrow(top_peaks))) {
      m <- top_peaks[i, ]
      a[[i]] <- list(
        x = m[["log2FoldChange"]],
        y = -log10(m[["padj"]]),
        text = rownames(m),
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40
      )
    }
    print("a")
    print(a)
    print(nrow(subset(na.omit(gene_list),group=="FALSE")))
    print(nrow(subset(na.omit(gene_list),group=="TRUE")))
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 2,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    ##Construct the plot object
    g = ggplot(data=na.omit(gene_list), aes(x=log2FoldChange, y=-log10(padj)) )+
      geom_point(aes(col=group,text=sprintf("Gene name: %s<br>p-value: %s", rownames(na.omit(gene_list)),padj)))+#alpha=0.4, size=1.75) +
      #geom_point(aes(col = factor(Significant))) +  
      #scale_color_manual(values=c("black", "red")) + 
      xlim(c(-1*scale_volx,scale_volx)) + 
      ylim(c(0,scale_voly)) +
      #scale_y_continuous(trans = "log1p")+
      xlab("log2 fold change") + ylab("-log10 p-value")+ labs(color = "Significant")+ theme_bw()#+
    # geom_text_repel(
    #   data = top_peaks,#subset(na.omit(gene_list),abs(log2FoldChange) > logfc_cutoff & padj < padj_cutoff),
    #   mapping=aes(label =rownames(top_peaks)),
    #   box.padding = unit(0.5, "lines"),
    #   point.padding = unit(0.5, "lines"),
    #   arrow = arrow(length = unit(0.01, 'npc')))
    # #g
    g1<-g+geom_text_repel(
      data = top_peaks,#subset(na.omit(gene_list),abs(log2FoldChange) > logfc_cutoff & padj < padj_cutoff),
      aes(label =rownames(top_peaks)),
      box.padding = unit(0.25, "lines"),
      point.padding = unit(0.1, "lines"),
      arrow = arrow(length = unit(0.01, 'npc')))
    
    # Increment the progress bar, and update the detail text.
    progress$inc(1/3, detail = paste("Doing part", 3,"/",3))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
   g0<-ggplotly(g)%>%layout(annotations = a)                        
    return(list(g1,g0))
    #g 
}
get_pheno<-function(data,pheno)
{
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
    return(pheno)
  }
}
#Enrichment function(returns object of enrichment function called)
#if Kegg anallysis is performed the enrichKegg function is called and its o/p is returned
enrichment_function<-function(enrichment_type,enrichment_input)
{
  obj<-NULL
    if(enrichment_type=="kegg")
    {
      obj<-enrichKEGG(gene         = enrichment_input[[1]],
                      organism     = enrichment_input[[2]],
                       pvalueCutoff = 1)
    }
    else if(enrichment_type=="biological process")
    {
      obj<-enrichGO(gene     = enrichment_input[[1]],
               universe      = enrichment_input[[2]],
               OrgDb         = enrichment_input[[3]],
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 1,
               qvalueCutoff  = 1, 
               readable      = TRUE)
      
    }
    else if(enrichment_type=="hallmark") 
    {
      
      obj<-enricher(enrichment_input[[1]],
                    TERM2GENE = enrichment_input[[2]],#c1_hallmark,
                    universe  = enrichment_input[[3]],
                    pvalueCutoff = 1,
                    #pAdjustMethod = "none",
                    qvalueCutoff = 1)
    }
    
}


#enrichment (result=DE_genes(),
             #organism=input$organism, 
             #dds.fc=batch_design()[[1]]
            #num<- length(input$combination),
            #wgcna_modules<-wgcna_output()$modules)
# WGCNA_matrix<-wgcna()[[2]]
# mod<- heat_wgcna()[[1]]

enrichment_main<-function(enrichment_type,result,input_organism,dds.fc,num,mod,WGCNA_matrix,c1_hallmark)
{
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  progress$set(message = "Processing data", value = 0)
  
  #get all de genes
  #result
  org<-NULL
  organism<-NULL
  universe<-NULL
  All_genes<-NULL
  #get the organism
  if(as.numeric(input_organism)==1)
  {
    org<-org.Hs.eg.db
    if(enrichment_type=="kegg") organism<-'hsa'
    else
    {
      All_genes=AnnotationDbi::select(org,rownames(assay(dds.fc)),"ENTREZID","SYMBOL",multiVals='first')
      
    }
    
  }
  else if (as.numeric(input_organism)==2)
  {
    org<-org.Mm.eg.db
    
    if(enrichment_type=="kegg") organism<-'mmu'
    
    else if (enrichment_type=="hallmark")
    {
      human = useMart(host = "jul2015.archive.ensembl.org",#"ensembl",
                      biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl")
      mouse = useMart(host = "jul2015.archive.ensembl.org",#"ensembl",
                      biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "mmusculus_gene_ensembl")
      
      All_genes=rownames(assay(dds.fc))
      genes_backgroud = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.character(All_genes), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
      print(head(genes_backgroud))
      universe_genes <- genes_backgroud[,2]
      
      entrez_background = bitr(universe_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      idx_all <- match(universe_genes, entrez_background$SYMBOL)
      universe<- entrez_background[idx_all,]$ENTREZID
    }
    else if (enrichment_type=="biological process")
    {
      All_genes=AnnotationDbi::select(org,rownames(assay(dds.fc)),"ENTREZID","SYMBOL",multiVals='first')
      
    }
    
  }
  
  #All genes
  #dds.fc
  #All_genes=rownames(assay(dds.fc))
  #all de comparisons
  #num <- length(combination())
  
  n<-num*2
  Enriched_list<-list()
  Enriched_obj<-list()
  Enriched_Kegg_gene<-list()#needed to display user selecetd kegg batch way.(contains a list of genes
                            #pertaining to the selected pathway and the corresponding foldchange values)
  k<-0
  
#if wgcna is performed then we perform enrichment analysis on the 
if(length(mod)>0)
    {
      modules<-as.data.frame(table(mod))
      n<-n+nrow(modules)
    }
#print(head(result)) 
for(i in 1:num) # looping through all DE comparisons
  {
    e_list<-list()
    e_obj<-list()
    kegg_genelist<-list() #used only for kegg
    for(j in 1:2) #two loops. First loop for up regulated de genes in a comparison 
    {                         #Second loop for down regulated de genes in a comparison
      k<-k+1
      print("k")
      print(num)
      print(i)
      res<-as.data.frame(result[[i]][j]) #get up/down regulated de genes for a comparison
      print(nrow(res))
      if(nrow(res)!=0)
      {
        genes<-res[,1]
        df<-res[,-1]
        rownames(df)<-genes
        print(head(df))
        
        entrez_id<-NULL
        gene_symbol<-NULL
        obj<-NULL
        geneList<-NULL #only popuated when kegg is computed
        if((enrichment_type=="hallmark") && (as.numeric(input_organism)==2))
        {
          genes_de = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.character(rownames(df)), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
          reg = bitr(genes_de[,2], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
          idx <- match(genes_de[,2], reg$SYMBOL)
          entrez<-reg[idx,]$ENTREZID
          #gene_symbol<-reg[idx,]$SYMBOL
          gene_symbol<-reg$SYMBOL[idx]
          entrez_id<-entrez[!is.na(entrez)]
        }
        else
        {
          #convert de gene names to entrez id(input format required to perform enrichment analysis) 
          reg=AnnotationDbi::select(org,rownames(df),"ENTREZID","SYMBOL",multiVals="first")
          idx <- match(rownames(df), reg$SYMBOL)
          df$entrez<-reg[idx,]$ENTREZID
          #gene_symbol<-reg[idx,]$SYMBOL
          gene_symbol<-reg$SYMBOL[idx]
          #print(gene_symbol)
          entrez_id<-df$entrez[!is.na(df$entrez)]
        }
        if(enrichment_type=="kegg")
          {
           input<-list(entrez_id,organism)
           obj<-enrichment_function("kegg",input)
           #compute gene list for kegg
           temp<-data.frame(entrez=c(reg[idx,]$ENTREZID),foldchange=c(df[idx,]$log2FoldChange))
           #temp<-temp[!is.na(temp)]
           #temp<-complete.cases(temp)
           temp<-temp[complete.cases(temp), ]
           #print(temp)
           geneList<-temp$foldchange
           #names(geneList)<-temp[,1]
           names(geneList)<-temp$entrez
          }
        else if(enrichment_type=="hallmark")
          {
            input<-list(entrez_id,c1_hallmark,universe)
            obj<-enrichment_function("hallmark",input)
          }
        else if(enrichment_type=="biological process") 
        {
          input<-list(entrez_id,All_genes$ENTREZID,org)
          obj<-enrichment_function("biological process",input)
        }
        
       #print(head(summary(kk)))
        print('obj')
        #print(obj)
        print('summary')
        #print(summary(obj))
        print(nrow(as.data.frame(summary(obj))))
        
        
        df_obj<-NULL
        if(!is.null(obj))
        {
          if(nrow(as.data.frame(summary(obj)))>0)
          {
            df_obj<-as.data.frame(summary(obj))[1:8]
            if(enrichment_type!="biological process")
            {
              df_obj<-as.data.frame(summary(obj))[1:8]
              print('df_obj')
              #print(df_obj)
              for(x in 1:length(df_obj[,8]))
              {
                print('kegg_sum')
                #print(typeof(x))
                #print(df_obj[x,8])
                temp<-strsplit(df_obj[x,8],"/")
                #print(temp)
                print('temp1')
                #print(temp[[1]])
                id<-which(entrez_id %in% temp[[1]] )
                #print(id)
                #print(gene_symbol)
                #print(paste(reg$SYMBOL[id], collapse = '/'))
                df_obj[x,8]<-paste(gene_symbol[id], collapse = '/')
                #print(sapply(reg$SYMBOL[id],function(y) paste(y,'/')))
            }
            
            }
            e_list[[length(e_list)+1]]<-df_obj
            e_obj[[length(e_obj)+1]]<-obj
            kegg_genelist[[length(kegg_genelist)+1]]<-geneList
          }
          else
          {
            e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
            e_obj[[length(e_obj)+1]]<-NULL
            kegg_genelist[[length(kegg_genelist)+1]]<-NULL
          }
          
        }
        else
        {
          e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          e_obj[[length(e_obj)+1]]<-NULL
          kegg_genelist[[length(kegg_genelist)+1]]<-NULL
        }
      }
      else
      {
        e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
        e_obj[[length(e_obj)+1]]<-NULL
        kegg_genelist[[length(kegg_genelist)+1]]<-NULL
      }
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", k,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
    }
    e_list[[length(e_list)+1]]<-data.frame(matrix(NA, nrow = 0, ncol = 8))
    e_obj[[length(e_obj)+1]]<-NULL
    kegg_genelist[[length(kegg_genelist)+1]]<-NULL
    #print(head(kegg))
    Enriched_list[[i]]<-e_list
    Enriched_obj[[i]]<-e_obj
    Enriched_Kegg_gene[[i]]<-kegg_genelist
  }
if(length(mod)>0) 
    {
      modules<-as.data.frame(table(mod))
      colnames(modules)<-c("Var1","number")
      #mod
      #WGCNA_matrix
      for(i in 1:nrow(modules))
      {
        #print(nrow(as.data.frame(result[[i]][1])))
        print('freq')
        print(modules$Freq[i])
        
        print(modules$Var1[i])
        idx_w<-which(mod==modules$Var1[i])
        print(head(idx_w))
        DEG<-colnames(WGCNA_matrix)[idx_w]
        print(head(as.data.frame(DEG)))
        
        e_list<-list()
        e_obj<-list()
        
        k<-k+1
        print("k")
        
        if(!identical(character(0),DEG))
        {
          entrez_id<-NULL
          gene_symbol<-NULL
          obj<-NULL
          
          if((enrichment_type=="hallmark") && (as.numeric(input_organism)==2))
          {
            genes_de = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = as.character(DEG), mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
            reg = bitr(genes_de[,2], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
            idx <- match(genes_de[,2],reg$SYMBOL)
            entrez<-reg[idx,]$ENTREZID
            gene_symbol<-reg$SYMBOL[idx]
            entrez_id<-entrez[!is.na(entrez)]
          }
          else
          {
            #convert de gene names to entrez id(input format required to perform enrichment analysis) 
            #get entrez id of genes 
            reg=AnnotationDbi::select(org,DEG,"ENTREZID","SYMBOL",multiVals="first")
            idx <- match(DEG,reg$SYMBOL) #match gene name to entrez id
            entrez<-reg[idx,]$ENTREZID
            gene_symbol<-reg$SYMBOL[idx]
            print(head(gene_symbol))
            print('entrez')
            print(entrez)
            entrez_id<-entrez[!is.na(entrez)]
          }
          if(enrichment_type=="kegg")
          {
            input<-list(entrez_id,organism)
            obj<-enrichment_function("kegg",input)
          }
          else if(enrichment_type=="hallmark")
          {
            input<-list(entrez_id,c1_hallmark,universe)
            obj<-enrichment_function("hallmark",input)
          }
          else if(enrichment_type=="biological process") 
          {
            input<-list(entrez_id,All_genes$ENTREZID,org)
            obj<-enrichment_function("biological process",input)
          }
          
          
          print('obj')
          print(obj)
          print('summary')
          print(summary(obj))
          print(nrow(summary(obj)))
          print(head(as.data.frame(summary(obj))[,8]))
          #print(nrow(as.data.frame(summary(kk))[,8]))
          
          if(!is.null(obj))
          {
            if(nrow(as.data.frame(summary(obj)))>0)
            {
              df_obj<-as.data.frame(summary(obj))[1:8]
              print('df_obj')
              print(df_obj)
              if(enrichment_type!="biological process")
              {
                for(x in 1:length(df_obj[,8]))
                {
                  print('kegg_sum')
                  #print(typeof(x))
                  print(df_obj[x,8])
                  temp<-strsplit(df_obj[x,8],"/")
                  #print(temp)
                  print('temp1')
                  #print(temp[[1]])
                  print(gene_symbol)
                  id<-which(entrez_id %in% temp[[1]] ) #get entrez id of genes in a pathway
                  #print(id)
                  print(paste(gene_symbol[id], collapse = '/'))
                  df_obj[x,8]<-paste(gene_symbol[id], collapse = '/')#convert entrez id to genes
                  #print(sapply(reg$SYMBOL[id],function(y) paste(y,'/')))
                }
              }
              
              d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
              Enriched_list[[length(Enriched_list)+1]]<-list(d,d,df_obj)
              Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,obj)
              Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
            }
            else
            {
              d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
              Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
              Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
              Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
            }
            
          }
          else
          {
            d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
            Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
            Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
            Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
          }
        }
        else
        {
          d<-data.frame(matrix(NA, nrow = 0, ncol = 8))
          Enriched_list[[length(Enriched_list)+1]]<-list(d,d,d)
          Enriched_obj[[length(Enriched_obj)+1]]<-list(NULL,NULL,NULL)
          Enriched_Kegg_gene[[length(Enriched_Kegg_gene)+1]]<-list(NULL,NULL,NULL)
        }
        # Increment the progress bar, and update the detail text.
        progress$inc(1/n, detail = paste("Doing part", k,"/",n))
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      
    }
  
  return(list(Enriched_list,Enriched_obj,Enriched_Kegg_gene))
}

#Display barplot of top 10 kegg pathway identified for selected comparison
#result<-Enriched_Kegg()[[2]] #obj
# res<-Enriched_Kegg()[[1]]
#input$category
#input$plot_k
#enrichment input
#category_go
#input$plotType_bp (plot_type)
#input$category_bp category
enrichment_plot<-function(enrichment_type,result,res,row,col,category,plot_type,category_go)
  {
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
  #p<-NULL
  
  #print(input$plot_k)
  #print(as.numeric(category))
  if(category!=""){
    print('heyho_ho')
    #print(result)
   
    #up regulated column
    up<-unique(nrow(res[[row]][[1]]))
    print('up')
    print(up)
    
    #down regulated column
    down<-unique(nrow(res[[row]][[2]]))
    
    #which column is null
    col_idx<-NULL
    if(length(up)==1 && up==0) col_idx<- 1
    else if (length(down)==1 && down==0) col_idx<-2
    #print(result)
    #kk<-result[[row]][[col]]
    print('colidx')
    print(col_idx)
    kk<-NULL
    
    if(enrichment_type!="biological process")
    {
      if(nrow(res[[row]][[col]])!=0)
      {
        if(plot_type=="Barplot") 
        {
          kk<-result[[row]][[col]]
          p<-barplot(kk, showCategory=as.numeric(category))
          
          # Increment the progress bar, and update the detail text.
          progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
          
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
          #print(head(x))
          return(p)
        }
      }
    }
    
    else
    {
      if(nrow(res[[row]][[col]])!=0)
      {
        ego<-result[[row]][[col]]
      }
      else ego<-NULL 
      
      print(ego)
      x<-NULL
      if(category_go!=""){
        x<-gofilter(ego,level=as.numeric(category_go))#dropGO(ego, level = as.numeric(input$category_go), term = NULL)
      }
      else x<-ego
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 1,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      #print(head(x))
      print(x)
      if(plot_type=="Barplot")
      {
        if(!is.null(x))
        {
          p<-barplot(x, showCategory=as.numeric(category))
        }
      }
      #else if(input$plotType_bp=="Dotplot")p<- dotplot(x,showCategory=as.numeric(input$category_bp))
      else if(plot_type=="Network plot" )
      {
        if(!is.null(x))
        {
          p<- enrichMap(x,n=as.numeric(category), vertex.label.cex=1)#, layout=igraph::layout.kamada.kawai)
        }
        
      }
      else if(plot_type=="GO term plot")
      {
        if(!is.null(x))
        {
          #x<-dropGO(ego, level = as.numeric(input$category_go), term = NULL)
          p<- plotGOgraph(x)
        }
        
      }
      # Increment the progress bar, and update the detail text.
      progress$inc(1/n, detail = paste("Doing part", 2,"/",n))
      
      # Pause for 0.1 seconds to simulate a long computation.
      Sys.sleep(0.1)
      return(p)
    }
    
    # if(is.null(col_idx))
    # {
    #   kk<-result[[row]][[col]]
    #   print(col)
    #   print(summary(kk))
    #   if(!is.null(kk))
    #   {
    #     if(input$plot_k=="Barplot") 
    #     {
    #       p<-barplot(kk, showCategory=as.numeric(input$category))
    #       p
    #     }
    #     #else p<-dotplot(ck)#(kk, showCategory=as.numeric(input$category))
    #   }
    # }
    # else
    # {
    #   if(col!=col_idx)
    #   {
    #     #print(result)
    #     #print(result[[row]][[1]])
    #     kk<-result[[row]][[1]]
    #   }
    #   print(col)
    #   #print(result[[row]][[col]])
    #   kk<-result[[row]][[col]] 
    #   if(!is.null(kk))
    #   {
    #     print('hey_kk')
    #     print(kk)
    #     if(input$plot_k=="Barplot" && !is.null(kk)) 
    #     {
    #       p<-barplot(kk, showCategory=as.numeric(input$category))
    #       p
    #     }
    #   }
    #   
    # }
    
  }
    
}
####Compute transcription factors from a list
Enriched_trasncription_factors<-function(TF_list,result,num,input_organism,mod,
                                         WGCNA_matrix,anova_table,conchoice,dds.fc)
{
  if (!is.null(TF_list)){  
    #get all de genes
    
    DE_TF<-list()
    for(i in 1:num)
    {
      TF<-list()
      for(j in 1:3)
      {
        #print(head(DE_TF[[i]][j]))
        print("k")
        res<-as.data.frame(result[[i]][j])
        genes<-res[,1]
        df<-res[,-1]
        rownames(df)<-genes
        print(head(df))
        
        if(as.numeric(input_organism)==1)
        {
          print(which(rownames(df) %in% TF_list$Human))
          TF[[length(TF)+1]]<-df[which(rownames(df) %in% TF_list$Human),]
        }
        else if (as.numeric(input_organism)==2)
        {
          print(which(rownames(df) %in% TF_list$Mouse))
          print(head(df[which(rownames(df) %in% TF_list$Mouse),]))
          TF[[length(TF)+1]]<-df[which(rownames(df) %in% TF_list$Mouse),]
          #print(head(DE_TF[[i]][j]))
        }
        print(head(TF))
      }
      DE_TF[[i]]<-TF
    }
    if((length(mod)>0))
    {
      modules<-as.data.frame(table(mod))
      colnames(modules)<-c("Var1","Freq")
      print(modules)
      #########preparing the anova table in the order as output#########
      a_tab<-anova_table[,-c(2,3)]
      cond<-unique(colData(dds.fc)[,as.numeric(conchoice)])
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
      all_genes=a_tab[,c(temp,c[2],c[3],temp2,c[1])]
      ##################################################################
      for(i in 1:nrow(modules))
      {
        print('freq')
        print(modules$Freq[i])
        
        print(modules$Var1[i])
        idx_w<-which(mod==modules$Var1[i])
        print(head(idx_w))
        DEG<-colnames(WGCNA_matrix)[idx_w]
        print(head(as.data.frame(DEG)))
        if(as.numeric(input_organism)==1)
        {
          print(which(DEG %in% TF_list$Human))
          genelist<-DEG[which(DEG %in% TF_list$Human)]
          ########
          print(head(all_genes))
          anova_genes<-rownames(all_genes)
          #anova_genes<-rownames(all_genes[which(rownames(all_genes) %in% TFs),])
          print(head(which(anova_genes %in% genelist)))
          print(head(all_genes[which(anova_genes %in% genelist),]))
          df_final<-all_genes[which(anova_genes %in% genelist),]
          #######
          DE_TF[[num+i]]<-list(character(0),character(0),genelist,df_final)
        }
        else if (as.numeric(input_organism)==2)
        {
          print(which(DEG %in% TF_list$Mouse))
          print("op")
          print(head(DEG[which(DEG %in% TF_list$Mouse)]))
          genelist<-DEG[which(DEG %in% TF_list$Mouse)]
          # if(op!=character(0))
          # {
          ########
          print(head(all_genes))
          anova_genes<-rownames(all_genes)
          #anova_genes<-rownames(all_genes[which(rownames(all_genes) %in% TFs),])
          print(head(which(anova_genes %in% genelist)))
          print(head(all_genes[which(anova_genes %in% genelist),]))
          df_final<-all_genes[which(anova_genes %in% genelist),]
          #######
          DE_TF[[num+i]]<-list(character(0),character(0),genelist,df_final)
          
        }
        
      }
      
      
    }
    return(DE_TF)
  }
}
#Heatmap function for either DE genes/DE transcription factors and ANOVA
heatmap_genes<-function(heatmap_call,dds,dds.fc,rld,heat_choice,
                        DE_genes,mod,WGCNA_matrix,num,heatmap_name,
                        Distance,Linkage)
{
  
  topVarGenes=NULL
  print("function line 1186")
  print(nrow(rld))
  
  if(heatmap_call=="ANOVA")
  {
    resSig<-results(dds)
    
    topVar<-head(rownames(resSig[order(resSig$padj),]),1000)
    topVarGenes<-which(rownames(rld) %in% topVar)
    
  }
  else if(heatmap_call=="DE")
  {
    if(!is.null(heat_choice))
    {
      result<-DE_genes()
      DEG<-NULL
      
      if(as.numeric(heat_choice)>num)
      {
        modules<-as.data.frame(table(mod))
        colnames(modules)<-c("Var1","numbers")
        print(modules$Var1[as.numeric(heat_choice)-num])
        idx_w<-which(mod==modules$Var1[as.numeric(heat_choice)-num])
        print("inside function line 1161")
        print(head(idx_w))
        DEG<-colnames(WGCNA_matrix)[idx_w]#all module genes
        
        print(head(as.data.frame(DEG)))
      }
      else
      {
        print("heatmap function line 1169")
        print(head(result[[as.numeric(heat_choice)]][[3]]))
        print(head(rownames(as.data.frame(result[[as.numeric(heat_choice)]][[3]]))))
        print(head(result[[as.numeric(heat_choice)]][[3]][,1]))
        #DEG<-rownames(result[[as.numeric(heat_choice)]][[3]])#as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]
        DEG<-as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]
        }
      print("DEG")
      print(length(DEG))
      print(typeof(DEG))
      print(length(DEG))
      topVarGenes<-which(rownames(rld) %in% DEG)
      print(head(topVarGenes))
      print(length(topVarGenes))
      print(length(topVarGenes)==0)
      
      # if(length(topVarGenes)==0)
      # {
      #   print(typeof(result[[as.numeric(heat_choice)]][[3]]))
      #   DEG<-as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]#as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]
      #   
      # }
        
    }
    
  }
  else if(heatmap_call=="TF")
  {
    if(!is.null(heat_choice))
    {
      print("line 1207")
      result<-DE_genes()
      #print(head(result))
      DEG<-NULL
      
      if(as.numeric(heat_choice)>num)
      {
        
        DEG<-result[[as.numeric(heat_choice)]][[3]]#colnames(WGCNA_matrix)[idx_w]#all module genes
        
        print(head(as.data.frame(DEG)))
      }
      else
      {
        print("heatmap function line 1211")
        print(head(result[[as.numeric(heat_choice)]][[3]]))
        print(head(rownames(as.data.frame(result[[as.numeric(heat_choice)]][[3]]))))
        print(head(result[[as.numeric(heat_choice)]][[3]][,1]))
        DEG<-rownames(result[[as.numeric(heat_choice)]][[3]])#as.data.frame(result[[as.numeric(heat_choice)]][[3]])[,1]
      }
      print("DEG")
      print(length(DEG))
      print(typeof(DEG))
      print(length(DEG))
      topVarGenes<-which(rownames(rld) %in% DEG)
      print(head(topVarGenes))
     
    }
    
  }
  validate(need(!is.null(topVarGenes),"The number of differentially expressed genes for this condition is Zero"))
  library(RColorBrewer)
  my_palette <- colorRampPalette(c("blue","white", "red"))(n = 15)
  print(my_palette)
  #Color key
  library(RColorBrewer)
  colors <-brewer.pal(8, "Dark2")
  n <- 30
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  mat<-rld[topVarGenes,]
  print(typeof(mat))
  print(head(mat))
  print(nrow(mat))
  p<-NULL
  if(nrow(mat)>1 && !is.null(nrow(mat)))
  {
  #mat<-mat-rowMeans(mat)
  mat <- t(scale(t(mat))) # scale and center rows
  #mat <- scale(mat) # scale and center columns
  
  dist_method<-c("euclidean", "manhattan")
  print(Distance)
  print(Linkage)
  print(dist_method[as.numeric(Distance)])
  
  distance = dist(mat,method = dist_method[as.numeric(Distance)])
  link_method<-c("average", "complete","Ward.D2","Ward.D","single")
  print(link_method[as.numeric(Linkage)])
  print("functions line 1305")
  print(head(distance))
  cluster = hclust(distance, method = link_method[as.numeric(Linkage)])
  
  #library('gplots')
  names<-FALSE
  if(nrow(mat)<20) names<-TRUE#rownames(mat)
  
  library(pheatmap)
  
  annotation <- data.frame(condition = colData(dds.fc)[,"condition"])
  rownames(annotation) <-  colnames(mat)# check out the row names of annotation
  
  Var1        <- unique(colors[colData(dds.fc)[,"condition"]])
  print(Var1)
  print(colData(dds.fc)[,"condition"])
  names(Var1) <- levels(colData(dds.fc)[,"condition"])
  print(Var1)
  anno_colors <- list(condition = Var1)
  
  p<-pheatmap(mat, annotation = annotation,color = my_palette,
              main = heatmap_name,
              clustering_distance_rows = dist_method[as.numeric(Distance)],
              clustering_method = link_method[as.numeric(Linkage)],
              annotation_colors = anno_colors,
              show_rownames=names)#,
}
  return(p)
}
#function to check for empty lists
has_empty_list <- function(x) {
  if(is.list(x)) {
    if (length(x)==0) {
      return(TRUE)
    } else {
      return(any(vapply(x, has_empty_list, logical(1))))
    }
  } else {
    return(FALSE)
  }
}

#source of variation output
source_of_variation_op<-function(edata,pData,ids)
{
  # edata<-assay(dds.fc)
  # pData<-colData(dds.fc)
  print(head(pData))
  sov<-as.data.frame(matrix(NA,nrow=length(ids),ncol=2))#ncol(pData)-1
  colnames(sov)<-c("Variation","Percentage")
  for(i in 1:length(ids))#2:ncol(pData))
  {
    print(colnames(pData))
    sov[i,1]<-colnames(pData)[ids[i]]
    sov[i,2]<-source_of_variation_ip(edata,pData,ids[i])
  }
  print(sov)
  return(sov)
} 

source_of_variation_ip<-function(edata,pData,k)
{
  #print(head(colSums(edata)))
  #print(length(edata))
  y_bar<-sum(rowSums(edata))/length(edata)
  #print(y_bar)
  #SSB(Sum of square between groups) and
  #SSW(Sum of squares within each group) calculation
  SSB=0
  SSW=0
  temp<-c()
  print(levels(pData[,k]))
  #k is the annotation variable
  #Loop through each group in k
  for(i in levels(pData[,k]))
  {
    #get which samples belong to a particular group
    id<-which(pData[,k]==i)
    print(id)
    #group mean
    y_i<-mean((rowSums(edata[,id])))
    #print(y_i)
    group<-0
    
    #loop through each sample in a group
    for(j in id)
    {
      p<-sapply(edata[,j],function(idx) (y_i-idx)^2)
      #print(head(p))
      group<-group+sum(p)
      #print("group")
      #print(group)
    }
    SSW<-SSW+group
    temp<-c(temp,(y_i-y_bar)^2)
  }
  print(SSW)
  SSB<-sum(temp)*length(edata)
  print(SSB)
  TSS<-SSB+SSW
  #print(TSS)
  R<-SSB/TSS
  R_Square=R^2
  print("r2")
  #print(R_Square)
  return(R_Square)
}
plotCounts<-function(batch_corrected=NULL,dds, gene, intgroup = "condition", normalized = TRUE, 
                     transform = FALSE, main, xlab = "group", returnData = FALSE, 
                     replaced = FALSE)
{
  #print(dds)
  stopifnot(length(gene) == 1 & (is.character(gene) | (is.numeric(gene) & 
                                                         (gene >= 1 & gene <= nrow(dds)))))
  stopifnot(returnData | all(sapply(intgroup, function(v) is(colData(dds)[[v]], 
                                                             "factor"))))
  if (is.null(sizeFactors(dds)) & is.null(normalizationFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  cnts <- counts(dds, normalized = normalized, replaced = replaced)[gene,]
  print(cnts)
  group <- if (length(intgroup) == 1) {
    colData(dds)[[intgroup]]
  }
  data <- data.frame(count = cnts + 0.5, group = as.integer(group))
  
  if(!is.null(batch_corrected))
  {
    cnts <- batch_corrected
    data <- data.frame(count = cnts, group = as.integer(group))
  }
  if (transform) {
    data$count <- log2(data$count)
    ylab <- expression(log[2] ~ count)
    logxy <- ""
  }
  else {
    ylab <- ifelse(normalized, "normalized count", "count")
    logxy <- "y"
  }
  if (missing(main)) {
    main <- if (is.numeric(gene)) {
      rownames(dds)[gene]
    }
    else {
      gene
    }
  }
  if (returnData)
    return(data.frame(count = data$count, colData(dds)[intgroup]))
  plot(data$group + runif(ncol(dds), -0.05, 0.05), data$count,
       xlim = c(0.5, max(data$group) + 0.5), log = logxy, xaxt = "n",
       xlab = xlab, ylab = ylab, main = main)
  axis(1, at = seq_along(levels(group)), levels(group))
}

gene_counts<-function(dataset,dds.fc,gene_name,log_scale)
{
  geneCounts <- plotCounts(dataset,dds.fc, gene=gene_name, intgroup="condition", returnData=TRUE)
  p<-ggplot(geneCounts, aes(x=condition, y=count,fill=condition)) +
    geom_boxplot()+
    geom_point(position=position_jitter(width=.1,height=0), size=3)+
    ggtitle(paste0("Gene name: ",gene_name)) +
    labs(x="Condition",y="Batch corrected counts") +
    theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=30, hjust=0.5)) +
    theme(axis.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22))+
    theme_bw()
  if(log_scale=="Yes") p+scale_y_continuous(trans='log2')
  else p
}

summary_analysis<-function(result_de,result_tf,result_kegg,result_bp,result_hall,combination)
{
  
    # result_de<-DE_genes()
    # result_tf<-DE_TF()
    # result_bp<-Enriched_BP()[[1]]
    # result_kegg<-Enriched_Kegg()[[1]]
    num <- length(combination)
    res_de<-data.frame(matrix(NA, nrow = num, ncol = 3))
    rownames(res_de)<-lapply(1:num, function(i) {
      combination[[i]]
    })
    res_tf<-data.frame(matrix(NA, nrow = num, ncol = 3))
    rownames(res_tf)<-lapply(1:num, function(i) {
      combination[[i]]
    })
    res_bp<-data.frame(matrix(NA, nrow = num, ncol = 2))
    rownames(res_bp)<-lapply(1:num, function(i) {
      #paste(combination()[[i]][1],' vs ',combination()[[i]][2])
      combination[[i]]
      
    })
    res_kegg<-data.frame(matrix(NA, nrow = num, ncol = 2))
    rownames(res_kegg)<-lapply(1:num, function(i) {
      combination[[i]]
      
    })
    res_hall<-data.frame(matrix(NA, nrow = num, ncol = 2))
    rownames(res_hall)<-lapply(1:num, function(i) {
      combination[[i]]
      
    })
    colnames(res_de)<-c('up-regulated_genes','down-regulated_genes','both')
    colnames(res_tf)<-c('up-regulated_TF','down-regulated_TF','All DE TF')
    colnames(res_bp)<-c('Enriched BP for Up-reg genes','Enriched BP for Down-reg genes')
    colnames(res_kegg)<-c('Enriched Kegg pathways for Up-reg genes','Enriched Kegg pathways for Down-reg genes')
    colnames(res_hall)<-c('Enriched Hallmark gene sets for Up-reg genes','Enriched Hallmark gene sets for Down-reg genes')
    print(res_de)
    print(res_de[1,1])
    print(res_tf)
    print(res_tf[1,1])
    print(res_bp)
    print(res_bp[1,1])
    print(res_kegg)
    print(res_kegg[1,1])
    print(res_hall)
    print(res_hall[1,1])
    for(i in 1:num)
    {
      #de
      print(nrow(as.data.frame(result_de[[i]][1])))
      res_de[i,1]<-nrow(as.data.frame(result_de[[i]][1]))
      res_de[i,2]<-nrow(as.data.frame(result_de[[i]][2]))
      res_de[i,3]<-nrow(as.data.frame(result_de[[i]][3]))
      #tf
      print(nrow(as.data.frame(result_tf[[i]][[1]])))
      res_tf[i,1]<-nrow(as.data.frame(result_tf[[i]][[1]]))
      res_tf[i,2]<-nrow(as.data.frame(result_tf[[i]][[2]]))
      res_tf[i,3]<-nrow(as.data.frame(result_tf[[i]][[3]]))
      
      print(nrow(as.data.frame(result_bp[[i]][[1]])))
      res_bp[i,1]<-nrow(as.data.frame(result_bp[[i]][[1]]))
      res_bp[i,2]<-nrow(as.data.frame(result_bp[[i]][[2]]))
      
      print(nrow(as.data.frame(result_kegg[[i]][[1]])))
      res_kegg[i,1]<-nrow(as.data.frame(result_kegg[[i]][[1]]))
      res_kegg[i,2]<-nrow(as.data.frame(result_kegg[[i]][[2]]))
      
      print(nrow(as.data.frame(result_hall[[i]][[1]])))
      res_hall[i,1]<-nrow(as.data.frame(result_hall[[i]][[1]]))
      res_hall[i,2]<-nrow(as.data.frame(result_hall[[i]][[2]]))
    }
    return(list(res_de,res_tf,res_kegg,res_bp,res_hall))
  
  
}

#all p value plot
p_value_all<-function(result,p_choice,combination)
{
  
  
  res<-as.data.frame(result[[p_choice]][4])
  print(head(res))
  print(head(res$pvalue))
  #i<-as.numeric(input$p_choice)
  
  p<-hist(res$pvalue[res$baseMean > 1],main = paste("Corrected p-value",combination[[p_choice]]),  col = "lavender", xlab = "p-values")
  return(p)
}  
  
  
  
