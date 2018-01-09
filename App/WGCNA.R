library(WGCNA)
source("./functions.r")
WGCNA_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    bsModal("modalwgcna", "Sample clustering", "go_wgcna",size = "large",
            selectInput(ns("Choose_top"), label = h5("Enter input choice"), 
                        choices = list("Consider top 5000 genes" = 1,"Consider the number of top genes specified by user"=2,
                                       "Select genes to be considered"=3,
                                       "Enter a list of genes to be considerd"=4),
                        selected = 1),
             uiOutput(ns("wgcna_gene_options")),
            # uiOutput("wgcna_action"),
            actionButton(ns("go_wgcna_start"),"Start"),
            conditionalPanel(sprintf("input['%s']>0",ns("go_wgcna_start")),
                             helpText("Enter cut off after cluster plot is displayed. 
                            Enter cut off only if there appears to be an outlier and click ok."),
                             uiOutput("cutoff_wgcna_ok"),
                             # actionButton("wgcna_ok","Ok"),
                             textInput(ns("cutoff_wgcna"), "Enter cut-off", value = "", width = NULL, placeholder = NULL),
                             plotOutput(ns("wgcna")),
                             helpText("Enter threshold and close this tab"),
                             textInput(ns("cutoff_soft"), "Enter threshold", value = "", width = NULL, placeholder = NULL),
                             plotOutput(ns("soft") )
                             )
            
            ),
    conditionalPanel(sprintf("input['%s']>0",ns("go_wgcna_start")),
                     actionButton(ns("go_view"),"View soft threshold plot"),
                     selectInput(ns("Choose_merge"), label = h5("Display heatmap of:"), 
                                 choices = list("All modules" = 1,"Modules after merging"=2),
                                 selected = 1),
    textInput(ns("cutoff_wgcna_dendo"),
              "enter cutoff", value = 0.25),
    #choices=nodes, multiple=FALSE, selected = h)
    plotOutput(ns("wgcna_module")),
    plotOutput(ns("heat_wgcna")),
    actionButton(ns("vis_tab"),"Visualize modules as a network"),
    downloadButton(ns('download_network'),'download network of all modules'),
    DT::dataTableOutput(ns("table_wgcna_display")),
    downloadButton(ns('download_wgcna_anova_Table'),'download ANOVA table for selected module'),
    DT::dataTableOutput(ns("filtered_data_wgcna")),
    gene_count_module_UI(ns("module6"))
    ),
    
    bsModal("modal_vis","Visualize modal",ns("vis_tab"),size="large",
            
            DT::dataTableOutput(ns("table_wgcna")),
            uiOutput(ns("network_options")),
            uiOutput(ns("network_hub")),
            uiOutput(ns("network_threshold")),
            visNetworkOutput(ns("plot_input"))
            )
    
  )
}
WGCNA_module<-function(input,output,session,infile,perform_voom,
                       condition,dds.fc,batch_design,conchoice,organism,ok3,combination,anova_table,
                       batch_choice,batch_corrected,normal)
  
{
  #
    datExpr<-infile()
    print(head(datExpr))
    print(typeof(datExpr))
    wgcna_click<-FALSE
  #}
  #choices
  #top 5000 or 3000 or 1000 or genes specified by the user
  observeEvent(input$Choose_top,
               {
                 output$wgcna_gene_options<-renderUI({
                   
                   if(as.numeric(input$Choose_top)==2)
                   {
                     textInput(session$ns("top_num_genes"), "Enter the number of top genes to be considerd", value = "", width = NULL, placeholder = NULL)
                     
                   }
                   else if(as.numeric(input$Choose_top)==3)
                   {
                     
                     selectInput(inputId = session$ns("Genes"),
                                 label = "Select genes",
                                 choices=rownames(datExpr), multiple=TRUE, selectize=FALSE)
                   }
                   else if(as.numeric(input$Choose_top)==4)
                   {
                     textInput(session$ns("top_gene_list"), "Enter a list of genes", value = "", width = NULL, placeholder = NULL)
                   }
                   #else 
                   
                   
                 })
                 output$wgcna_action<-renderUI({
                   actionButton(ns("go_wgcna_start"),"Start")
                 })
               })
  #Computes WGCNA
  wgcna<-reactive({
    #print(unlist(strsplit(input$top_num_genes,"\\s+")[[1]])) remove space or tab probably use  a grep expression
    
    limit<-5000
    gene_list<-NULL
    if(as.numeric(input$Choose_top)==2)
    {
      limit<-as.numeric(input$top_num_genes)
    }
    else if(as.numeric(input$Choose_top)==3)
    {
      gene_list<-input$Genes
    }
    else if(as.numeric(input$Choose_top)==4)
    {
      gene_list<-unlist(strsplit(input$top_gene_list,"\\s+")[[1]])
      print("gene list")
    }
    print("count")
    #print(head(assay(dds.fc()[[1]])))
    datExpr0<-NULL
    if(perform_voom==TRUE) 
      {
      library(limma)
      print("inside wgcna line 105")
      print(head(datExpr))
      design <- model.matrix(~0+colData(datExpr)$condition)
      #Estimate size factors  
      dds.norm=estimateSizeFactors(datExpr)
      print(head(counts(dds.norm, normalized=TRUE)))
      print(design)
      v <- voom(counts=counts(dds.norm, normalized=TRUE), design=design)
      print(v)
       datExpr0<-t(v$E)
       print(head(datExpr0))
    }
    else datExpr0<-t(datExpr)
    print("before")
    print(ncol(datExpr0))
    if(length(gene_list)>0)
    {
      print(gene_list)
      print(typeof(gene_list))
      datExpr0<-datExpr0[,which(colnames(datExpr0) %in% gene_list)]
      limit<-ncol(datExpr0)
    }
    print('after')
    print(ncol(datExpr0))
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    progress$set(message = "Processing Data", value = 0)
    
    n<-2
    print("cut of")
    print(input$cutoff_wgcna)
    if((input$cutoff_wgcna)!="")
    {
      temp<-datExpr0
      sampleTree = hclust(dist(datExpr0), method = "average");
      print("ok")
      #Determine cluster under the line
      clust = cutreeStatic(sampleTree, cutHeight = as.numeric(input$cutoff_wgcna), minSize = 10)
      table(clust)
      # clust 1 contains the samples we want to keep.
      keepSamples = (clust==1)
      datExpr0 = temp[keepSamples, ]
      
    }
    # Increment the progress bar, and update the detail text.
    progress$inc(1/(n), detail = paste("Doing part", 2,"/",2))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    
    dat<-t(datExpr0)
    WGCNA_matrix = t(dat[order(apply(dat,1,mad), decreasing = T)[1:limit],])
    print('wgcna')
    #print(head(WGCNA_matrix))
    # Increment the progress bar, and update the detail text.
    progress$inc(1/(n), detail = paste("Doing part", 2,"/",2))
    # Pause for 0.1 seconds to simulate a long computation.
    Sys.sleep(0.1)
    #print(head(WGCNA_matrix))
    list(datExpr0,WGCNA_matrix)
  })
  #Compute heatmap of WGCNA
  #Package used is WGCNA
  #Please refer to WGCNA module for the algorithm underlying the computation of heatmap
  heat_wgcna<-reactive({
    req(input$go_wgcna_start)
    
    if(input$go_wgcna_start>0)
    {
      WGCNA_matrix<-wgcna()[[2]]
      
      # Create a Progress object
      progress <- shiny::Progress$new()
      # Make sure it closes when we exit this reactive, even if there's an error
      on.exit(progress$close())
      progress$set(message = "Processing Data", value = 0)
      
      n<-2
      METree<-NULL
      if(input$cutoff_soft!="")
      {
        softPower = as.numeric(input$cutoff_soft);
        adjacency = adjacency(WGCNA_matrix, power = softPower);
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency);
        dissTOM = 1-TOM
        #  
        # Call the hierarchical clustering function
        geneTree = hclust(as.dist(dissTOM), method = "average");
        
        # We like large modules, so we set the minimum module size relatively high:
        minModuleSize = 30;
        # Module identification using dynamic tree cut:
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 4, pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize);
        print("dynamicmod")
        #print(table(dynamicMods)) #table
        #print(dynamicMods)
        # Convert numeric lables into colors
        dynamicColors = labels2colors(dynamicMods) #print table
        #print(table(dynamicColors))
        
        # Calculate eigengenes
        MEList = moduleEigengenes(WGCNA_matrix, colors = dynamicColors)
        MEs = MEList$eigengenes
        MEColors= dynamicColors
        print("MES")
        # print(MEs)
        # Calculate dissimilarity of module eigengenes
        MEDiss = 1-cor(MEs);
        # Cluster module eigengenes
        METree = hclust(as.dist(MEDiss), method = "average");
        if(input$Choose_merge==2) #merge modules
        {
          #merge modules
          #merging modules whose expressions are similiar
          
          MEDissThres = 0.25
          if(!is.null(input$cutoff_wgcna_dendo)) MEDissThres = as.numeric(input$cutoff_wgcna_dendo)
          #  # Plot the cut line into the dendrogram
          #  abline(h=MEDissThres, col = "red")
          # Call an automatic merging function
          merge = mergeCloseModules(WGCNA_matrix, dynamicColors, cutHeight = MEDissThres, verbose = 3)
          # The merged module colors
          # mergedColors = merge$colors;
          MEColors = merge$colors;
          # Eigengenes of the new merged modules:
          # mergedMEs = merge$newMEs;
          MEs= merge$newMEs;
          # Calculate dissimilarity of module eigengenes
          mergedMEDiss = 1-cor(merge$newMEs);
          # Cluster module eigengenes
          #merged
          METree = hclust(as.dist(mergedMEDiss), method = "average");
          
          
        }
        # Increment the progress bar, and update the detail text.
        progress$inc(1/(n), detail = paste("Doing part", 1,"/",2))
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        #heatmap
        #dds.fc<-batch_design()[[1]]
        #pheno<-colData(dds.fc)
        
        traitData <- data.frame(matrix(0, ncol=length(unique(condition)), nrow = length(rownames(WGCNA_matrix))),row.names = rownames(WGCNA_matrix))
        colnames(traitData) <- unique(condition)
        #print(traitData)
        for (i in 1:length(colnames(traitData))){traitData[which(condition==colnames(traitData)[i]),colnames(traitData)[i]] <-1}
        #print(head(traitData))
        
        #Define numbers of genes and samples
        nGenes = ncol(WGCNA_matrix);
        nSamples = nrow(WGCNA_matrix);
        # Recalculate MEs with color labels
        #MEs0 = moduleEigengenes(WGCNA_matrix, mergedColors)$eigengenes
        MEs0=MEs#mergedMEs
        MEs = orderMEs(MEs0)
        moduleTraitCor = cor(MEs, traitData, use = "p");
        moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
        # Will display correlations and their p-values
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
        dim(textMatrix) = dim(moduleTraitCor)
        # Increment the progress bar, and update the detail text.
        progress$inc(1/(n), detail = paste("Doing part", 2,"/",2))
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
        list(MEColors,moduleTraitCor,traitData,MEs,textMatrix,TOM,METree)
      }
    }
    
    
  })
#to repoen soft threshold tab
  observeEvent(input$go_view, {
    toggleModal(session, "modalwgcna", toggle = "close")
  })
  
observeEvent(input$go_wgcna_start,
             {
               output$cutoff_wgcna_ok<-renderUI({
                 
                 datExpr<-wgcna()[[1]]
                 sampleTree = hclust(dist(datExpr), method = "average")
                 #########
                 dend2<-as.dendrogram(sampleTree)
                 print(get_nodes_attr(dend2,"height"))
                 nodes<-c(unique(get_nodes_attr(dend2,"height")))+0.25 # node's height
                 nodes<-c(0,nodes)
                 print(nodes)
                 #ggd1 <- as.ggdend(dend2)
                 selectInput(inputId = session$ns("cutoff_wgcna"),
                             label = "Select cutoff",
                             choices=nodes, multiple=FALSE, selected = 0)
               })
               #PLot dendogram 
               output$wgcna<- renderPlot({
                 
                     datExpr<-wgcna()[[1]]
                     sampleTree = hclust(dist(datExpr), method = "average")
                     # h<-0.25
                     # if(as.numeric(input$cutoff_wgcna)!=0) h<-as.numeric(input$cutoff_wgcna)
                     # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
                     # The user should change the dimensions if the window is too large or too small.
                     plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
                          cex.axis = 1.5, cex.main = 2)
                     abline(h=0.25, col = "red")
                   
               })
               
               #display threshold plot
               output$soft<-renderPlot({
                 #Construction of co-expression network
                 #  #similarity measure between gene profiles: biweight midcorrelation
                 WGCNA_matrix<-wgcna()[[2]]
                     
                     # Create a Progress object
                     progress <- shiny::Progress$new()
                     # Make sure it closes when we exit this reactive, even if there's an error
                     on.exit(progress$close())
                     progress$set(message = "Processing Data", value = 0)
                     
                     n<-2
                     s = abs(bicor(WGCNA_matrix))
                     
                     # Increment the progress bar, and update the detail text.
                     progress$inc(1/(n), detail = paste("Doing part", 1,"/",2))
                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     powers = c(c(1:10), seq(from = 12, to=100, by=2))
                     sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
                     
                     # Increment the progress bar, and update the detail text.
                     progress$inc(1/(n), detail = paste("Doing part", 2,"/",2))
                     # Pause for 0.1 seconds to simulate a long computation.
                     Sys.sleep(0.1)
                     
                     plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                          xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
                          type='n', main = paste('Scale independence'));
                     text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
                          labels=powers,cex=1,col='red'); abline(h=0.90,col='red')
                  
               })
               
               #Display heatmap
               output$heat_wgcna<-renderPlot({
                 if(input$cutoff_soft!="" && input$go_wgcna_start>0)
                 {
                   # Display the correlation values within a heatmap plot
                   labeledHeatmap(Matrix = heat_wgcna()[[2]],
                                  xLabels = names(heat_wgcna()[[3]]),
                                  yLabels = names(heat_wgcna()[[4]]),
                                  ySymbols = names(heat_wgcna()[[4]]),
                                  colorLabels = FALSE,
                                  colors = blueWhiteRed(50),
                                  textMatrix = heat_wgcna()[[5]],
                                  setStdMargins = FALSE,
                                  cex.text = 0.5,
                                  zlim = c(-1,1),
                                  main = paste("Module-trait relationships"))
                 }
               })
               output$cutoff_options<-renderUI({
                 if(input$go_wgcna_start>0)
                 {
                   METree<-heat_wgcna()[[7]]
                   if(!is.null(METree))
                   {
                     #library("dendextend")
                     dend2<-as.dendrogram(METree)
                     #print(get_nodes_attr(dend2,"height"))
                     nodes<-c(unique(get_nodes_attr(dend2,"height")))+0.25 # node's height
                     #print(nodes)
                     #ggd1 <- as.ggdend(dend2)
                     h<-0.25
                     if(!is.null(v_pval$click_combo_list)) h<-v_pval$click_combo_list
                     selectInput(inputId = session$ns("cutoff_wgcna_dendo"),
                                 label = "Select cutoff",
                                 choices=nodes, multiple=FALSE, selected = h)
                   }
                 }
                 
               })
               #reactiveValues()
               v_pval <- reactiveValues(
                 click_combo_list = NULL)
               # # Handle clicks on the plot
               observeEvent(input$cutoff_wgcna_dendo,{
                 isolate({
                   print('hey')
                   v_pval$click_combo_list = input$cutoff_wgcna_dendo
                   #print(length(v_pval$click_combo_list)+1)
                   print('ok')
                 })
               })
               output$wgcna_module<-renderPlot({
                 
                   METree<-heat_wgcna()[[7]]
                   if(!is.null(METree))
                   {
                     
                     plot(METree, main = "Clustering of module eigengenes",
                          xlab = "", sub = "")
                     abline(h=as.numeric(input$cutoff_wgcna_dendo), col = "red")#0.25
                     
                   }
                 
                 #else {plotly_empty()}
               })
               #Display computes modules and the number of genes present in each module
               output$table_wgcna <- DT::renderDataTable({
                 if(input$cutoff_soft!="" )
                 {
                   wgcna_click=TRUE
                   moduleColors=heat_wgcna()[[1]]
                   moduleColors<-moduleColors[which(moduleColors!="grey")]
                   df<-as.data.frame(table(moduleColors))
                   colnames(df)<-c("Modules","Number of genes")
                   #print(head(df))
                   datExpr<-wgcna()[[2]]
                   
                   #print(moduleColors)
                   hubs    = chooseTopHubInEachModule(datExpr,moduleColors )
                   print(hubs)
                   #print(cbind(df,hubs))
                   df<-as.data.frame(cbind(df,hubs))
                   colnames(df)[3]<-"Hub Genes"
                   DT::datatable(t(df),class = 'cell-border stripe',
                                 selection = list(target = 'column',mode='single'),
                                 extensions = list('Scroller'=NULL,'Buttons'=NULL),
                                 options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                                buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),
                                                                            text = 'Download only genes on clipboard'))))
                 }
               })
               
               
               #Display computes modules and the number of genes present in each module
               
               output$table_wgcna_display <- DT::renderDataTable({
                 if(input$cutoff_soft!="" )
                 {
                   wgcna_click=TRUE
                   moduleColors=heat_wgcna()[[1]]
                   moduleColors<-moduleColors[which(moduleColors!="grey")]
                   df<-as.data.frame(table(moduleColors))
                   colnames(df)<-c("Modules","Number of genes")
                   #print(head(df))
                   datExpr<-wgcna()[[2]]
                   
                   #print(moduleColors)
                   hubs    = chooseTopHubInEachModule(datExpr,moduleColors )
                   print(hubs)
                   #print(cbind(df,hubs))
                   df<-as.data.frame(cbind(df,hubs))
                   colnames(df)[3]<-"Hub Genes"
                   DT::datatable(t(df),class = 'cell-border stripe',
                                 selection = list(target = 'column',mode='single'),
                                 extensions = list('Scroller'=NULL,'Buttons'=NULL),
                                 options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                                buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),
                                                                            text = 'Download only genes on clipboard'))))
                 }
               })
               
               
               network_cytoscape<- reactive({
                 
                 #Input preparation
                 datExpr<-t(wgcna()[[2]])
                 # Recalculate topological overlap if needed
                 TOM = heat_wgcna()[[6]]#TOMsimilarityFromExpr(datExpr, power = as.numeric(input$cutoff_soft));#heat_wgcna()[[6]]
                 # print(head(TOM))
                 #get modules
                 print('mod')
                 df<-as.data.frame(table(heat_wgcna()[[1]]))
                 colnames(df)<-c("Modules","Number of genes")
                 
                 inp<-c(1:nrow(df))
                 threshold=0.2
                 observeEvent(input$vis_tab,
                              {
                                threshold=0.2
                                inp<-input$table_wgcna_columns_selected
                                req(input$network_threshold)
                                threshold<-as.numeric(input$network_threshold)
                              })
                 # Select modules
                 modules = df$Modules[inp]#c("brown", "red");
                 # print(modules)
                 # Select module probes
                 probes = rownames(datExpr)
                 # print(head(probes))
                 moduleColors=heat_wgcna()[[1]]
                 moduleColors=moduleColors[which(moduleColors!="grey")]
                 
                 group_fc<-list()
                 cyt<-list()
                 edgefile<-NULL
                 nodefile<-NULL
                 file_name<-NULL
                 for(i in modules)
                 {
                   print(i)
                   inModule = is.finite(match(moduleColors, i));
                   print(head(inModule))
                   print(threshold)
                   modProbes = probes[inModule];
                   length(modProbes)
                   # Select the corresponding Topological Overlap
                   modTOM = TOM[inModule, inModule];
                   dimnames(modTOM) = list(modProbes, modProbes)
                   print('modTom')
                   #print(head(modTOM))
                   # Export the network into edge and node list files Cytoscape can read
                   cyt[length(cyt)+1] = exportNetworkToCytoscape(modTOM,
                                                                 edgeFile = paste("CytoscapeInput-edges-", paste(i, collapse="-"), ".txt", sep=""),
                                                                 nodeFile = paste("CytoscapeInput-nodes-", paste(i, collapse="-"), ".txt", sep=""),
                                                                 weighted = TRUE,
                                                                 threshold = threshold,
                                                                 nodeNames = modProbes,
                                                                 #altNodeNames = modGenes,
                                                                 nodeAttr = moduleColors[inModule])
                   #print(cyt)
                   edgefile[length(edgefile)+1]<-paste("CytoscapeInput-edges-", paste(i, collapse="-"), ".txt", sep="")
                   nodefile[length(nodefile)+1]<- paste("CytoscapeInput-nodes-", paste(i, collapse="-"), ".txt", sep="")
                   file_name[length(file_name)+1]<-paste("CytoscapeInput_GFC_", paste(i, collapse="-"), ".xlsx", sep="")
                   print("ok")
                   print(cyt[length(cyt)])
                   if(ok3()>0)
                   {
                     ####calculate group fold change
                     a_tab<-anova_table()[,-c(2,3)]
                     cond<-unique(colData(dds.fc())[,as.numeric(conchoice())])
                     print(cond)
                     c<-colnames(a_tab)
                     print(length(c))
                     temp<-as.vector(c[4:(3+length(cond))])
                     temp2<-as.vector(c[(length(cond)+4):length(c)])
                     print(temp2)
                     # #temp<-as.vector(c[7:length(c)])
                     print(c(temp,c[2],c[3],temp2,c[1]))
                     print('howdy line 592 inside wgcna module')
                     # print(c)
                     # anova <- 
                     #print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
                     a_tab<-a_tab[,c(temp,c[2],c[3],temp2,c[1])]
                     print(head(a_tab))
                     overall_mean_idx<-ncol(a_tab)
                     print(overall_mean_idx)
                     dds.fc<-batch_design()
                     print(colData(dds.fc()))
                     dds<-colData(dds.fc())
                     idx<-c(1:length(unique(dds.fc()$condition)))#ids of mean per condition
                     #get number of comparisons
                     combo<-combination()
                     comb<-length(combo())
                     idx_fc<-c(length(unique(dds.fc()$condition))+3)#foldchange ids
                     if(comb>1)
                     {
                       idx_fc<-seq(length(unique(dds.fc()$condition))+3,(comb*6)+1,6)#foldchange ids
                     }
                     idx_gfc<-c()
                     #3+length(unique(dds))+1+((row-1)*6)
                     library(gtools)
                     for(j in 1:length(unique(dds.fc()$condition)))
                     {
                       GFC <- foldchange(a_tab[,j],a_tab[,overall_mean_idx])#group mean/overall mean
                       a_tab <- cbind(a_tab, round(GFC,3))
                       idx_gfc<-c(idx_gfc,ncol(a_tab))
                       colnames(a_tab)[ncol(a_tab)] <- paste("GFC_",levels(dds.fc()$condition)[j])
                     }
                     print('wgcna at_tab')
                     print(head(a_tab))
                     idx<-c(idx,idx_gfc,idx_fc)
                     #idx<-c(1:length(unique(dds$condition)),overall_mean_idx+1:overall_mean_idx+length(unique(dds$condition)))
                     print(idx)
                     print(head(a_tab[which(rownames(a_tab) %in% modProbes),idx]))
                     #group_fc[[length(group_fc)+1]]<-1
                     #print(group_fc)
                     group_fc[[length(group_fc)+1]]<-a_tab[which(rownames(a_tab) %in% modProbes),idx]
                     print(group_fc)
                     #get the list of transcription factors
                     
                     TF_list<-read.csv("./www/Transcriptome_TFcat.txt", header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
                     group_fc[[length(group_fc)]]$Transcription_factor<-rep('No',nrow(group_fc[[length(group_fc)]]))
                     print(group_fc)
                     if(as.numeric(organism())==1)#humans
                     {
                       idx<-which(rownames(group_fc[[length(group_fc)]]) %in% TF_list$Human)
                       group_fc[[length(group_fc)]]$Transcription_factor[idx]<-rep('Yes',length(idx))
                     }
                     else if (as.numeric(organism())==2)#mouse
                     {
                       idx<-which(rownames(group_fc[[length(group_fc)]]) %in% TF_list$Mouse)
                       group_fc[[length(group_fc)]]$Transcription_factor[idx]<-rep('Yes',length(idx))
                     }
                   }
                 }
                 
                 list(cyt,edgefile,nodefile,group_fc,file_name,modules)
               })
               
               
               #Download the network for the modules selected
                 output$download_network <- downloadHandler(
                   
                   filename =function()
                   {
                     paste("output", "zip", sep=".")
                     #paste(input$plotType_hall,' of Up regulated hallmarkP for ',condition,'.svg')
                     
                   },
                   content = function(file) {
                     
                     
                     fs<-c()
                     df<-as.data.frame(table(heat_wgcna()[[1]]))
                     colnames(df)<-c("Modules","Number of genes")
                     # Select modules
                     modules_num<- ncol(df)
                     network<-network_cytoscape()
                     if(!is.null(network[[1]]))
                     {
                       for(i in 1:modules_num)#input$table_wgcna_columns_selected
                       {
                         #files2zip <- dir(network_cytoscape()[[6]][[i]], full.names = TRUE)
                         #files2zip<-c(files2zip,network_cytoscape()[[2]][[i]],network_cytoscape()[[3]][[i]],network_cytoscape()[[5]][[i]])
                         # edge_data<-network[[1]][[i]]
                         # 
                         # write.xlsx2(edge_data, file=network[[2]][[i]], sheetName = "Sheet1",
                         #             col.names = TRUE, row.names = TRUE, append = FALSE)
                         # 
                         # node_data<-network[[1]][[i]]
                         # print(head(node_data))
                         # write.xlsx2(node_data, file=network[[3]][[i]], sheetName = "Sheet1",
                         #             col.names = TRUE, row.names = TRUE, append = FALSE)
                         fs<-c(fs,network[[2]][[i]],
                               network[[3]][[i]],
                               network[[5]][[i]])
                         #fs<-c(fs,files2zip)
                         #fs<-c(fs,network_cytoscape()[[5]][[i]])
                         #print(head(network))
                         network[[1]][[i]]
                         
                         #print(is.empty(network_cytoscape()[[4]]))
                         print(has_empty_list(network[[4]]))
                         # print((network_cytoscape()[[4]][[i]])!=list())
                         if(!has_empty_list(network[[4]]))#if(!is.null(network_cytoscape()[[4]][[i]]))
                         {
                           print(i)
                           print(network[[4]][[i]])
                           #write.table(network_cytoscape()[[4]][[i]], file=network_cytoscape()[[5]][[i]], row.names = T, col.names=NA, quote = F, sep = "\t")
                           
                           write.xlsx2(network[[4]][[i]],
                                       file=network[[5]][[i]],
                                       sheetName = "Sheet1",
                                       col.names = TRUE, row.names = TRUE, append = FALSE)
                         }
                         
                       }
                     }
                     zip(zipfile=file, files=fs)
                     
                     
                   }
                 )
               
               
               observeEvent(input$table_wgcna_columns_selected,{
                 print("line 507 wgcna")
                 
                 print(input$table_wgcna_columns_selected)
                 
                 
                 
                 output$network_options<-renderUI({
                   if(length(input$table_wgcna_columns_selected)>0)
                   {
                     selectInput(session$ns("network_options"),label="Choice of network",
                                 choices = list("Network constructed usng top hub genes"=1,
                                                "Network constructed using threshold"=2),selected = 2)
                   }
                 })
                 
                 output$network_threshold<-renderUI({
                   req(input$network_options)
                   
                   if(as.numeric(input$network_options)==2)
                   {
                     textInput(session$ns("network_threshold"),
                               "Enter the threshold to be considered to build the network for a module",
                               value = 0.2, width = NULL, placeholder = NULL)
                     
                   }
                 })
                 
                 output$network_hub<-renderUI({
                   req(input$network_options)
                   
                   if(as.numeric(input$network_options)==1)
                   {
                     textInput(session$ns("network_hub"),
                               "Enter the to number of hub genes to be considered to build the network for a module",
                               value = 30, width = NULL, placeholder = NULL)
                     
                   }
                 })
                 print(input$table_wgcna_columns_selected)
                 
                 #install.packages("visNetwork")
                 output$plot_input<-renderVisNetwork({
                   # for(i in 1:length(input$table_wgcna_columns_selected))#input$table_wgcna_columns_selected
                   # {
                   i<-input$table_wgcna_columns_selected
                   req(input$network_threshold)
                       if(input$network_threshold>0)
                       {
                     #edgefile<-network_cytoscape()[[2]][[i]]
                     #nodefile<-network_cytoscape()[[3]][[i]]
                    # print(head(network_cytoscape()))
                     edge_table<-network_cytoscape()[[1]][[1]]
                     node<-unique(edge_table[,1])
                     #print(head(edge))
                     print(head(node))
                     print(head(edge_table$fromNode))
                     #print(edge_table$toNode)
                     library(visNetwork)
                     nodes <- data.frame(id = node, title = paste(node), 
                                         shape = rep("dot", length(node)))
                                         #size = 10:15, color = c("blue"))#, "red"))
                     edges <- data.frame(from = edge_table$fromNode, to = edge_table$toNode)
                     print(head(nodes))
                     print(dim(nodes))
                     print(head(edges))
                     print(dim(edges))
                     visNetwork(nodes, unique(edges)) %>%
                       visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
                     
                     
                   }
                 })
                 
               })   
               observeEvent(input$table_wgcna_display_columns_selected,{
                 if(ok3()>0)
                 {
                 #display anova table
                 row<-input$table_wgcna_display_columns_selected
                   
                       mod<-heat_wgcna()[[1]]
                       modules<-as.data.frame(table(mod))
                       colnames(modules)<-c("Var1","number")
                       
                       print(table(mod))
                       WGCNA_matrix<-wgcna()[[2]]
                       
                       print(modules$Var1[row])
                         #print(head(colnames(WGCNA_matrix)))
                         print(head(mod))
                         print(typeof(mod))
                         #print(which(mod == modules$Var1[row-num] ))
                         idx_w<-which(mod==modules$Var1[row])
                         print(head(idx_w))
                         gene_list<-colnames(WGCNA_matrix)[idx_w]
                         print(head(as.data.frame(gene_list)))
                         
                         
                         #########preparing the anova table in the order as output#########
                         a_tab<-anova_table()[,-c(2,3)]
                         cond<-unique(colData(dds.fc())[,as.numeric(conchoice())])
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
                         print(head(all_genes))
                         anova_genes<-rownames(all_genes)
                         #anova_genes<-rownames(all_genes[which(rownames(all_genes) %in% TFs),])
                         print(head(which(anova_genes %in% gene_list)))
                         print(head(all_genes[which(anova_genes %in% gene_list),]))
                         df<-all_genes[which(anova_genes %in% gene_list),]
                         
              output$filtered_data_wgcna <- DT::renderDataTable({
                         DT::datatable(df,class = 'cell-border stripe',#as.data.frame(gene_list)
                                       selection = list(mode='single',target = 'row'),
                                       extensions = list('Scroller'=NULL,'Buttons'=NULL),
                                       options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                                      buttons = list('copy', list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),
                                                                                  text = 'Download only genes on clipboard'))))
                         
                       
                 
               })
              #download button
              output$download_wgcna_anova_Table <- downloadHandler(
                
                filename = function() 
                {
                  paste('ANOVA genes for ',modules$Var1[row],'.xlsx')
                  
                },
                content = function(file) {
                  nam<-paste('ANOVA genes for ',modules$Var1[row])
                  
                  write.xlsx2(df, file, sheetName = nam,
                              col.names = TRUE, row.names = TRUE, append = FALSE)
                }
                
              )
                 ###plotcounts of gene
                 #display gene expression across conditon
                 observeEvent(input$filtered_data_wgcna_rows_selected,{
                   #get which gene was clicked
                   #get the row  number of the gene clicked
                   print('hey')
                   print(input$filtered_data_wgcna_rows_selected)
                   #print(input$ANOVA_rows_clicked)
                   selected_row <- input$filtered_data_wgcna_rows_selected
                   print(selected_row)
                   #get the normalized data
                   full_data<-normal()
                   print(head(full_data))
                   print('hey')
                   an_gene<-rownames(df)[selected_row]
                   
                   print(an_gene)
                   counts<-as.vector(full_data[which(an_gene %in% rownames(full_data)),])
                   #print(counts)
                   #print(as.factor(colData(dds.fc()[[1]])[,as.numeric(input$conchoice)]))
                   cond<-as.vector(colData(dds.fc())[,as.numeric(conchoice())])
                   #rownames(cond)<-colData(dds.fc()[[1]])[,1]
                   library('data.table')
                   print(length(counts))
                   print(length(cond))
                   df<-data.frame(counts,cond)
                   colnames(df)<-c('count','condition')
                   print(head(df))
                   print(batch_choice())
                   print(as.numeric(batch_choice()))
                   batch<-batch_choice()
                   if(as.numeric(batch)==1) callModule(gene_count_module,"module6",NULL,reactive({dds.fc()}),reactive({an_gene}))
                   else
                   {
                     idx<-which(rownames(batch_corrected()) %in% an_gene)
                     print(idx)
                     callModule(gene_count_module,"module6",batch_corrected()[idx,],reactive({dds.fc()}),reactive({an_gene}))
                   }
                 })
                 }
                   
               })
                # observeEvent(input$network_threshold,
                #  {
                #    req(input$network_threshold)
                #    if(input$network_threshold>0) plot_input()
                #  }
                # )
               
               
             })  
#if() 
return(list(
            #wgcna_start=reactive({input$go_wgcna_start}),
            #soft_thershold=reactive({input$cutoff_soft}),
            modules=reactive({heat_wgcna()[[1]]}),
            WGCNA_matrix=reactive({wgcna()[[2]]})
            #wgcna_click=reactive({wgcna_click}))
       )
)
  
  }
# if(input$ok3>0)
# {
#   ####calculate group fold change
#   a_tab<-anova_table()[,-c(2,3)]
#   cond<-unique(colData(dds.fc()[[1]])[,as.numeric(conchoice)])
#   print(cond)
#   c<-colnames(a_tab)
#   print(length(c))
#   temp<-as.vector(c[4:(3+length(cond))])
#   temp2<-as.vector(c[(length(cond)+4):length(c)])
#   print(temp2)
#   # #temp<-as.vector(c[7:length(c)])
#   print(c(temp,c[2],c[3],temp2,c[1]))
#   print('howdy')
#   # print(c)
#   # anova <- 
#   #print(head(a_tab[,c(temp,c[2],c[3],temp2,c[1])]))
#   a_tab<-a_tab[,c(temp,c[2],c[3],temp2,c[1])]
#   print(head(a_tab))
#   overall_mean_idx<-ncol(a_tab)
#   print(overall_mean_idx)
#   dds.fc<-batch_design()[[1]]
#   print(colData(dds.fc))
#   dds<-colData(dds.fc)
#   idx<-c(1:length(unique(dds$condition)))#ids of mean per condition
#   #get number of comparisons
#   comb<-length(input$combination)
#   idx_fc<-c(length(unique(dds$condition))+3)#foldchange ids
#   if(comb>1)
#   {
#     idx_fc<-seq(length(unique(dds$condition))+3,(comb*6)+1,6)#foldchange ids
#   }
#   idx_gfc<-c()
#   #3+length(unique(dds))+1+((row-1)*6)
#   library(gtools)
#   for(j in 1:length(unique(dds$condition)))
#   {
#     GFC <- foldchange(a_tab[,j],a_tab[,overall_mean_idx])#group mean/overall mean
#     a_tab <- cbind(a_tab, round(GFC,3))
#     idx_gfc<-c(idx_gfc,ncol(a_tab))
#     colnames(a_tab)[ncol(a_tab)] <- paste("GFC_",levels(dds$condition)[j])
#   }
#   print('wgcna at_tab')
#   print(head(a_tab))
#   idx<-c(idx,idx_gfc,idx_fc)
#   #idx<-c(1:length(unique(dds$condition)),overall_mean_idx+1:overall_mean_idx+length(unique(dds$condition)))
#   print(idx)
#   print(head(a_tab[which(rownames(a_tab) %in% modProbes),idx]))
#   #group_fc[[length(group_fc)+1]]<-1
#   #print(group_fc)
#   group_fc[[length(group_fc)+1]]<-a_tab[which(rownames(a_tab) %in% modProbes),idx]
#   print(group_fc)
#   
#   #get the list of transcription factors
#   
#   TF_list<-read.csv("./www/Transcriptome_TFcat.txt", header = TRUE,sep = "\t",check.names = FALSE,quote = "\"")
#   group_fc[[length(group_fc)]]$Transcription_factor<-rep('No',nrow(group_fc[[length(group_fc)]]))
#   print(group_fc)
#   if(as.numeric(input$organism)==1)#humans
#   {
#     idx<-which(rownames(group_fc[[length(group_fc)]]) %in% TF_list$Human)
#     group_fc[[length(group_fc)]]$Transcription_factor[idx]<-rep('Yes',length(idx))
#   }
#   else if (as.numeric(input$organism)==2)#mouse
#   {
#     idx<-which(rownames(group_fc[[length(group_fc)]]) %in% TF_list$Mouse)
#     group_fc[[length(group_fc)]]$Transcription_factor[idx]<-rep('Yes',length(idx))
#   }
# }