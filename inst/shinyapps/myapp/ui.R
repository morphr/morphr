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
      textInput("Color_File", label = p("Color File"), value = "Used to plot curves"),
      actionButton("file_ready", "Done"),
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
      titlePanel("Shape Analysis"),
      fluidRow(
        column(3,uiOutput("plots")),
        column(9,plotly::plotlyOutput("geo_dist", width = "100%", height="600px"),
               plotOutput("mean_shape"),
               plotOutput("PCA_plot"),
               plotOutput("mds_plot"),
               plotOutput("dendogram_plot"),
               plotOutput("dfa_plot"))
      )
      #plotly::plotlyOutput("geo_dist", width = "100%", height="600px"),
      #plotOutput("mean_shape"),
      #plotOutput("PCA_plot"),
      #plotOutput("mds_plot"),
      #plotOutput("dendogram_plot"),
      #plotOutput("dfa_plot")
      
      
      
      
      
      
    )
    
  )
)