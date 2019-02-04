#'@title view.database
#'@description view database contents made by ProcessHplc package functions
#'@param database character string name or filepath and name of database
#'@export

view.database <- function(my.data=''){

  #1. load libraries
  library('dplyr')
  library('dbplyr')
  library('RSQLite')
  library('shiny')
  library('plotly')

  #load in data and organise columns, add a 'time' column
  db <- src_sqlite(my.data, create=F)

  #start of shiny app section
  app = shinyApp(
    ui <- fluidPage(
      titlePanel("Database Viewer"),

      fluidRow(

        column(12,
               wellPanel(
                 selectInput('table',
                             'Select Table',
                             src_tbls(db),
                             'meta.data')
               )#end of wellPanel
        )),#end of column
        fluidRow(
        column(8,
               dataTableOutput('current.table'))

      )#end of fluid row 1
    ),#end of UI

    server <- function(input, output,session) {

      cur.tab <- reactiveValues()

      observeEvent(input$table,{
        cur.tab$tab <- as.data.frame(tbl(db, input$table))

      })

      output$current.table <- renderDataTable(
        cur.tab$tab
      )

      observeEvent(input$ending,{
        stopApp()
      })

      session$onSessionEnded(function(){
        stopApp()
      })

    }
  )
  runApp(app)
}
