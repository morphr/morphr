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
      selectInput("var", 
                  label = "What kind of curve?",
                  choices = c("Open", 
                              "Closed"),
                  selected = "Closed"),
      # Horizontal line ----
      tags$hr(),
      actionButton("goButton", "Geodesic distance"),
      p("Click to get a symmetrix matrix of geodesic distance."),
      actionButton("Mean_shape_plotButton", "Mean Shape"),
      p("Click to get mean shape."),
      tags$hr(),
      actionButton("PCA_plotButton", "PCA"),
      p("Click to get PCA plot"),
      actionButton("mds_plotButton", "MDS"),
      p("Click to get MDS plot"),
      actionButton("dendogram_plotButton", "Dendogram"),
      p("Click to get dendogram plot"),
      actionButton("dfa", "Deformation Field (All)"),
      p("Click to get Deformation Field (All) plot"),
      tags$hr(),
      textOutput("contents"),
      tags$hr(),
      textInput("mat_path", label = h3("file path input"), value = "Enter mat file input"),
      textInput("color_path", label = h3("Color code input"), value = "Enter csv input"),
      actionButton("file_ready", "Done"),
      hr(),
      fluidRow(column(3, verbatimTextOutput("value"))),
      tags$hr(),
      p("Uploaded plots"),
      uiOutput("plots")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      titlePanel("Shape Analysis"),
      tableOutput("geo_dist"),
      plotOutput("mean_shape"),
      plotOutput("PCA_plot"),
      plotOutput("mds_plot"),
      plotOutput("dendogram_plot"),
      plotOutput("dfa_plot")
      
      
      
      
      
      
    )
    
  )
)