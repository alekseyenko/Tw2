dashboardPage(skin = "purple",
  

  # Application title
  dashboardHeader( title = HTML(paste("Shiny T", tags$sub("w"), tags$sup(2), sep = ""))),

  dashboardSidebar(
    sidebarMenu(
       textOutput("dataLoaded"),

       menuItem("Dashboard", tabName = "dashboard"),

       actionButton("loadDemoDataBtn", "Load Demo Data"),
       actionButton("clearDataBtn", "Clear Data"),
       
       div(id = "testcontrols",
         selectInput("distanceMethod", "Distance", distance_choices, selected = "jsd"),
         uiOutput("strataFactor"),
         uiOutput("mainFactor") , 
         actionButton("addTestBtn", "Add Test"), 
         textOutput("numTests")
       ),

       conditionalPanel("output.numTests != '0'",
         numericInput("numPermutations", "Number of Permutations",
           value = 999, min = 1, max = 100000),
         actionButton("runTestsBtn", "Run Tests")),
       menuItem(text = "Citation", tabName = "help")
    )),
    
    dashboardBody( 
      tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css')),
      useShinyjs(),
      tabItems(
        tabItem(tabName = "dashboard",

          conditionalPanel("output.dataLoaded != 'data loaded'", 
            fluidRow(box(
                  fileInput("rawCountFile", "Raw Count File"),
                  fileInput("sampleDataFile", "Sample Data File")))),

          conditionalPanel("output.dataLoaded == 'data loaded'", 
            fluidRow(box(title = "Ordination Plot", width = 12, plotOutput("plot"))),

            conditionalPanel("output.numTests != '0'",
              fluidRow(box(width = 12, title = "Tests",
                       actionButton("deleteRowsBtn", "Delete Rows"), DTOutput("testTable")))
           ))),
       tabItem(tabName = "help", "Citation and help not ready")
   )))
