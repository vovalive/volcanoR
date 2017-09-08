#devtools::install_github('andrewsali/shinycssloaders')
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)

### orgs
org=data.frame(row.names = c('Mus musculus','Rattus norvegicus','Homo sapiens'))
org$code=c('mmu','rno','hsa')
org$db=c("org.Mm.eg.db","org.Rn.eg.db","org.Hs.eg.db")
org=as.matrix(org)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  tags$head(
    tags$style(
      HTML(
        "h3 {color:purple;}
         img {box-shadow: 0 0 1em grey;}
        #dark {color:navy;}"
      )
    )
  ),
  # Application title
  titlePanel("VolcanoR"),
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p('Start with loading your DiffExpression results. We accept tab-separated files with columns: GeneSymbol Pvalue Log2FoldChange
        You can first download the sample',
        a(href = 'https://raw.githubusercontent.com/vovalive/volcanoR/master/res2_GSE70213.txt', 'res.txt'),
        'file for mouse muscular tissue differential expresion results, and then try uploading it.'),
      selectInput("org", "Organism:", 
                  choices=row.names(org)),
      fileInput('file1', 'Upload file',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
      p('when your file is uploaded, you can click button below to get Volcano plot'),
      actionButton("plot", "Generate volcano plot"),
      p('you can adjust tresholds, it will change the look of your volcano plot and enrichment results'),
      sliderInput("fc.tr",
                  "log2 Fold change treshold:",
                  min = 0,
                  max = 8,
                  step=0.01,
                  value = 3),
      sliderInput("p.tr",
                  "-log10 P-value treshold:",
                  min = 0,
                  max = 10,
                  step=0.01,
                  value = 4),
      
      actionButton("go", "Do enrichment test"),
      p('It can take a while, please be patient')
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Start here", 
                 h3('Welcome to VolcanoR web sevice'),
                 h4("Short howto:"),
                 h5("1. Prepare your differential expression results in tab separated format, you can see example",a(href = 'https://raw.githubusercontent.com/vovalive/volcanoR/master/res2_GSE70213.txt', 'HERE')),
                 img(src="samplefile.png", height = 100),
                 h5("2. Upload your file on left panel, and click 'Generate volcano plot' button"),
                 h5("3. Now you can go to the Results tab and see generated plot"),
                 h5("4. You can adjust tresholds and regenerate plot "),
                 h5("5. When you like list of genes you selected, you can do KEGG and GO enrichment test by clicking button"),
                 img(src="sampleplot.png", width = 300)
                 ),
        tabPanel("RESULTS",
          withSpinner(plotOutput("distPlot",height = '700px')),
          p(),
          verbatimTextOutput("pv"),
         # actionButton("copyButton", "Copy!"),
         h4("KEGG enrichment analysis"),
         column(12,withSpinner(dataTableOutput('tablekegg'),proxy.height = '80px')),
         h4("GeneOntology biological process enrichment analysis"),
          column(12,withSpinner(dataTableOutput('tablego'),proxy.height = '80px')),
         h4("Selected genes"),
          column(12,dataTableOutput('table')))
      )
     )
)))

# Define server logic 
server <- shinyServer(function(input, output) {
  #(list=ls(all=T))
  #libs : 
  library(ggplot2)
  library(gplots)
  library(ggrepel)
  library(dplyr)
  library(clipr)
  library("AnnotationDbi")
  library("clusterProfiler")
  library("org.Mm.eg.db")
  library("org.Rn.eg.db")
  library("org.Hs.eg.db")
  
  #options:
  options(shiny.maxRequestSize = 9*1024^2)
  showmax=2000
  
  #load file and detect signiff
  ress=reactive({
    inFile <- input$file1
    r=read.csv(inFile$datapath,sep='\t',header = T)
    colnames(r)=c('name','pval','log2fc')
    fc=input$fc.tr
    pv=input$p.tr
    mutate(r, sig=ifelse((-log10(r$pval) >= pv & abs(r$log2fc) >= fc ), "Sig", "Not Sig"))
  })
  
  #calculate pval and fc
  output$pv  <- renderText({
    paste0("P-value treshold  = ", as.character(10**(-input$p.tr)),"\n","Fold change treshold  = ", as.character(2**(input$fc.tr)),"\nGenes selected:",length(which(ress()$sig == 'Sig')))
  })
  
  volc = eventReactive(input$plot, {
    # generate plot
    b=1  # point size
    fc=input$fc.tr
    pv=input$p.tr
    if(length(which(ress()$sig == 'Sig'))<=showmax)
    {ggplot(data = ress(),aes(x = log2fc,y = -log10(pval))) + geom_point(size=b,aes(col=sig)) + scale_color_manual(values=c("darkgrey","red")) + geom_vline(xintercept = c(-fc,fc)) + geom_hline(yintercept = c(0,pv)) +geom_text_repel(max.iter=10,data=filter(ress(), sig == 'Sig'), aes(label=name))} 
    else 
    {ggplot(data = ress(),aes(x = log2fc,y = -log10(pval))) + geom_point(size=b,aes(col=sig)) + scale_color_manual(values=c("darkgrey","red")) + geom_vline(xintercept = c(-fc,fc)) + geom_hline(yintercept = c(0,pv))}
  })
  
 # volcano plot with labels
  output$distPlot <- renderPlot({
    # generate plot
      volc()
    })
 
  # tanble of signiff genes 
 output$table <- renderDataTable(ress()[which(ress()$sig=='Sig'),],options = list(searching= F))
 
  # KEGG and GO(sloooow) enrichment analysis
  enrich=eventReactive(input$go, {
    or = org[which(row.names(org)==input$org),]
    dbname=as.character(or[2])
    orgname=as.character(or[1])
    r=ress()
    gene = r$name[which(r$sig=='Sig')]
    gene = bitr(gene, fromType="SYMBOL", toType="UNIPROT", OrgDb=dbname)
    eK=enrichKEGG(gene = gene$UNIPROT,keyType = 'uniprot',organism = orgname, pvalueCutoff = 0.05,pAdjustMethod = "BH")
    ego2 <- enrichGO(gene = gene$UNIPROT,OrgDb = dbname, keytype = 'UNIPROT',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.01)
    list(eK@result,ego2@result)
  })
  
  #render KEGG and GO tables
  output$tablekegg=renderDataTable(enrich()[[1]], options = list(pageLength = 5))
  output$tablego=renderDataTable(enrich()[[2]], options = list(pageLength = 5))
  })

# Run the application 
shinyApp(ui = ui, server = server)