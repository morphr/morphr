ui <- fluidPage(
  
  # App title ----
  titlePanel("Uploading Files"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      h3("File Uploading"),
      fileInput("file1", "Choose a ucf or svg file list",
                multiple = FALSE,
                accept = c("text/plain")),
      fileInput("Color_File", "Choose a csv file to plot curves",
                multiple = FALSE,
                accept = c("csv")),
      textOutput("contents"),
      tags$hr(),
      h3("Type of Curve"),
      selectInput("var", 
                  label = "What kind of curve?",
                  choices = c("Open", 
                              "Closed"),
                  selected = "Closed"),
      h3("Heat map of geodesic distance."),
      actionButton("goButton", "Geodesic distance"),
      h3("Plot Mean Shape"),
      actionButton("Mean_shape_plotButton", "Mean Shape"),
      h3("Plot PCA Plot Along Eigen Axis 1 and 2"),
      actionButton("PCA_plotButton", "PCA"),
      h3("Plot PCA Variation Along Eigen Axis 1"),
      actionButton("PCA_variation_plotButton", "PCA Variation"),
      h3("Plot MDS Plot"),
      actionButton("mds_plotButton", "MDS"),
      h3("Plot Dendogram"),
      actionButton("dendogram_plotButton", "Dendogram"),
      h3("Plot Deformation"),
      actionButton("dfa", "Deformation Field (All)"),
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