ui <- fluidPage(
  
  # App title ----
  titlePanel("Uploading Files"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose a ucf or svg file list",
                multiple = FALSE,
                accept = c("text/plain")),
      fileInput("Color_File", "Choose a csv file to plot curves",
                multiple = FALSE,
                accept = c("csv")),
      textOutput("contents"),
      tags$hr(),
      selectInput("var", 
                  label = "What kind of curve?",
                  choices = c("Open", 
                              "Closed"),
                  selected = "Closed"),
      actionButton("goButton", "Geodesic distance"),
      p("Click to get a symmetrix matrix of geodesic distance."),
      actionButton("Mean_shape_plotButton", "Mean Shape"),
      p("Click to get mean shape. Always get mean shape first"),
      actionButton("PCA_plotButton", "PCA"),
      p("Click to get PCA plot"),
      actionButton("PCA_variation_plotButton", "PCA Variation"),
      p("Click to get PCA variation along eigen axis 1"),
      actionButton("mds_plotButton", "MDS"),
      p("Click to get MDS plot"),
      actionButton("dendogram_plotButton", "Dendogram"),
      p("Click to get dendogram plot"),
      actionButton("dfa", "Deformation Field (All)"),
      p("Click to get Deformation Field (All) plot"),
      fluidRow(column(3, verbatimTextOutput("value")))
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      titlePanel("Shape Analysis Using Square Root Velocity Function"),
      plotOutput("plots"),
      plotly::plotlyOutput("geo_dist", width = "100%", height="600px"),
      plotOutput("mean_shape",width = "100%", height="600px"),
      plotOutput("PCA_plot",width = "100%", height="600px"),
      plotOutput("PCA_variation_plot",width = "100%", height="600px"),
      plotOutput("mds_plot",width = "100%", height="600px"),
      plotOutput("dendogram_plot",width = "100%", height="600px"),
      plotOutput("dfa_plot",width = "100%", height="600px")


      
      
      
      
      
    )
    
  )
)