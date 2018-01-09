# Shiny-Seq

It is an interactive web based application designed to assist biologists in exploration, visualization, and interpretation of RNA-Seq data.

## Getting Started

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5","tximport",'DESeq2','clusterProfiler',"org.Hs.eg.db","org.Mm.eg.db","org.Mmu.eg.db","sva","limma","geneplotter",'biomaRt',"pcaGoPromoter","pcaGoPromoter.Mm.mm9","pcaGoPromoter.Hs.hg19","pathview")

install.packages("shiny","shinyBS","shinyjs",'RColorBrewer',"stringr",'formula.tools','data.table','fdrtool',"VennDiagram",'colorspace',"xlsx",'svglite',"visNetwork","V8","ggrepel","ReporteRs","ReporteRsjars")

install.packages("gplots",dependencies = TRUE)

devtools::install_github("ropensci/plotly")

devtools::install_github("rstudio/crosstalk",force=TRUE)

devtools::install_github('rstudio/DT')

```


