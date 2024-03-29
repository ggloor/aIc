#
# This is a Shiny web application. 
#

library(shiny)
library(aIc)

ui <- fluidPage(

  # App title ----
  titlePanel("Compositional test"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(


     
      fileInput("upload", NULL, accept=c(".txt", ".tsv"), buttonLabel = "selex or Upload .tsv ..."), 
# hard coded to 2 lines
      numericInput("z.rem", "max 0 proportion", value = 0.95, min = 0, max=1, step = 0.05),
      tableOutput("files") ,
      numericInput("group_1_size", "group 1 size", value=7, min=3, step=1),
      numericInput("group_2_size", "group 2 size", value=7, min=3, step=1),
#      tableOutput('groups'),
      
      # Input: Selector for test ----
      selectInput("test", "Compositional test:",
                  c("distance dominance" = "dom",
                    "perturbation invariance" = "pert",
                    "scale invariance" = "scale",
                    "correlation coherence:spearman" = "cohere.sp",
                    "correlation coherence:pearson" = "cohere.pe",
                    "data singularity" = "sing")),

      # Input: Selector for normalization ----
      selectInput("norm", "Data normalization:",
                  c("proportion" = "prop",
                    "centred log ratio" = "clr",
                    "intraquartile log ratio" = "iqlr",
                    "low V, high abund log ratio" = "lvha",
                    "edgeR TMM" = "TMM",
                    "edgeR TMMwsp" = "TMMwsp",
                    "DESeq RLE" = "RLE",
                    "entropy" = "H",
                    "none" = "none")),

      # Input: Selector for zero treatment ----
      selectInput("zero", "Zero replacement:",
                  c("prior" = "prior",
                    "Geometric Baysian" = "GBM",
                    "lrSVD (slow)" = "lrSVD",
                   "none" = "none")),

 # Input: Selector for distance measure  ----
      selectInput("distance", "Distance measure:",
                  c("Euclidian (Aitchison if clr)" = "euclidian",
                    "Bray Curtis" = "bray",
                    "Jaccard" = "jaccard",
                    "Jensen-Shannon" = "jsd")),

      # Input: Checkbox for whether outliers should be included ----
      checkboxInput("log", "Take log of transform (prop, TMM, RLE)", F)

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Formatted text for caption ----
      h4(textOutput("caption")),

      # Output: Plot of the data ----
      plotOutput("testPlot"),
      
      
      tableOutput("groups"),
      tableOutput("head")
    )
  )
)

# Define server logic to plot various variables and load data ----
server <- function(input, output, session) {
  
  onStop(function(){ 
   
    # restore to default shiny.maxRequestSize
    options(shiny.maxRequestSize=5*1024^2)
    print("shiny maxRequestSize reset to default")
   
  })
    
  options(shiny.maxRequestSize=10*1024^2)

  # and only exclude outliers if requested
  output$testPlot <- renderPlot({

 up.data <- reactive({
    if(is.null(input$upload)){
      url <- 'https://raw.githubusercontent.com/ggloor/datasets/main/selex.tsv'
      tsv = read.table(url, ,header=T, row.names=1, sep = "\t")
    } else {

    req(input$upload)
    
    #ext <- tools::file_ext(input$upload$name)
    #switch(ext,
     # csv = vroom::vroom(input$upload$datapath, delim = ","),
      tsv = read.table(input$upload$datapath,header=T, row.names=1, sep = "\t")
   #   validate("Invalid file; Please upload a .txt or .tsv file")
   # )
    }
  })
    
  output$head <- renderTable({
    head(up.data(), 2)
  })
  
#  # set groups 
    group <- c(rep('A', input$group_1_size), rep("B", input$group_2_size))
    zero.method = input$zero
    
    if(input$zero=="none") {zero.method = NULL }
     
# perturbation 
    if(input$test=='pert'){
      x <- aIc.perturb(up.data(), norm.method=input$norm, zero.method=zero.method,
        zero.remove=input$z.rem, distance=input$distance, log=input$log, group=group)
      aIc.plot(x)
      if(x$is.perturb == 'Yes'){
        output$caption <- renderText({paste('The data are approximately perturbation invariant with transform ', input$norm,' and distance ', input$distance, ". The maximum observed relative perturbation is: ", x$ol , sep="")})
      } else if(x$is.perturb == 'No') {
        output$caption <- renderText({paste('The data are not perturbation invariant with transform ', input$norm,' and distance ', input$distance, '. The maximum observed relative perturbation is: ', x$ol, ' fold change. Please try the clr transform on this dataset.', sep="")})
      }

# dominance
    } else if (input$test=='dom'){
      x <- aIc.dominant(up.data(), norm.method=input$norm,zero.method=zero.method,
        zero.remove=input$z.rem,  log=input$log, distance=input$distance, group=group)
      aIc.plot(x)
      if(x$is.dominant == 'Yes'){
        output$caption <- renderText({paste('The data are distance dominant with transform ', input$norm, ' and distance ', input$distance, '. The proportion of non-dominant distances in the sub-compositon is: ', 1 - round(x$ol,2), "%.", sep="")})
      } else if(x$is.dominant == 'No') {
        output$caption <- renderText({paste('The data are not distance dominant with transform ', input$norm,' and distance ', input$distance, '. The proportion of non-dominant distances in the sub-compositon is: ', 1 - round(x$ol,2), '%. Please try the clr transform on this dataset.', sep="")})
      }
      
# scale
    } else if (input$test=='scale'){
      x <- aIc.scale(up.data(), norm.method=input$norm, zero.method=zero.method,
        zero.remove=input$z.rem, distance=input$distance, log=input$log, group=group)
      aIc.plot(x)
      if(x$is.scale == 'Yes'){
        output$caption <- renderText({paste('The data are scale invariant with transform ', input$norm,' and distance ', input$distance, ". The proportion of non-scale invariant distances in the sub-compositon is: ", round(x$ol,2), "%.", sep="")})
      } else if(x$is.scale == 'No') {
        output$caption <- renderText({paste('The data are not scale invariant with transform ', input$norm,' and distance ', input$distance, '. The proportion of non-scale consistent distances is: ', round(x$ol,2), ' fold change. Please try the clr transform on this dataset.', sep="")})
      }

# coherence      
    } else if (input$test=='cohere.sp' || input$test=='cohere.pe'){
      if(input$test=='cohere.sp') x <- aIc.coherent(up.data(), norm.method=input$norm,
        zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group,
       cor.test='spearman')
      if(input$test=='cohere.pe') x <- aIc.coherent(up.data(), norm.method=input$norm,
        zero.method=zero.method, zero.remove=input$z.rem,  log=input$log, group=group, 
        cor.test='pearson')
      aIc.plot(x)
      if(x$is.coherent == 'Yes' & input$norm != 'none'){
        output$caption <- renderText({paste('The data are sub-compositionally coherent with ', input$norm,". Network analysis and correlation inference may be OK.", sep="")})
      } else if(x$is.coherent == 'No' & input$norm != 'none') {
        output$caption <- renderText({paste('The data are not sub-compositionally coherent with transform ', input$norm,'. When doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375', sep="")})
      } else if(x$is.coherent == 'Yes' & input$norm == 'none') {
        output$caption <- renderText({paste('The data must be transformed in some way prior to use. You should use some form of compositional association test and should check out the propr R package. In addition, when doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375', sep="")})
      }

    # singularity      
    } else if (input$test=='sing'){
      x <- aIc.singular(up.data(), norm.method=input$norm, zero.method=zero.method,
        zero.remove=input$z.rem,  log=input$log, group=group)
      if(x$is.singular == 'Yes'){
       output$caption <- renderText({paste('The data are singular with transform ', input$norm,'. When doing dimension reduction you need to be aware of compositional effects; see Greenacre and Aitchison 2002 "Biplots of compositional data" JRSA 51:375. You likey will be best served using a compositional analysis approach.', sep="")})
      } else if (x$is.singular == 'No'){
       output$caption <- renderText({paste('The data are not singular with transform ', input$norm,'. Go about your daily business.', sep="")})
      }
    }
  })
}

# Create Shiny app ----
shinyApp(ui, server)