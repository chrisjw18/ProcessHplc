#'@title id.peaks
#'@description identify peaks in output of EXPORT3D.mac for ChemStation hplc data
#'@param my.data character string name or filepath and name of the CSV containing the original data
#'@param nups.451 number of ups to define peak at 451 nm (can be altered if peaks not captured correctly)
#'@param nups.223 number of ups to define peak at 223 nm (can be altered if peaks not captured correctly)
#'@return SQLite database in current working directory updated to output of integrate.peaks() with designated peak identifications determined during interactive peak ID shinyApp.
#'@export

id.peaks <- function(my.data = '', nups.451=10, nups.223=20){

  #access database made in integrate.peaks() for all info needed here
  if(my.data %in% dir()){
    db.nam <- my.data %>% strsplit(.,split='[.]') %>% lapply('[[',1)  %>% unlist
  } else {
    db.nam <- my.data %>% basename %>% strsplit(.,split='[.]') %>% lapply('[[',1)  %>% unlist
  }
  mydb <- src_sqlite(db.nam, create = F)
  print(paste('Opening database *',db.nam,'* in ',getwd(), sep=''))

  #pull out all the objects we need
  dat <- as.data.frame(tbl(mydb,'corrected_abs_data'))
  integrated.223 <- as.data.frame(tbl(mydb, 'integrated.223'))
  integrated.431 <- as.data.frame(tbl(mydb, 'integrated.431'))
  integrated.451 <- as.data.frame(tbl(mydb, 'integrated.451'))
  abs.spec.223.correct <- as.data.frame(tbl(mydb, '223.peaks.abs.spec.correct'))
  abs.spec.431.correct <- as.data.frame(tbl(mydb, '431.peaks.abs.spec.correct'))
  abs.spec.451.correct <- as.data.frame(tbl(mydb, '451.peaks.abs.spec.correct'))
  #lib.meta <- as.data.frame(tbl(mydb, 'library_meta_data'))
  #lib.abs.spec <- as.data.frame(tbl(mydb, 'library_abs_data'))
  #lib.vite.spec <- as.data.frame(tbl(mydb, 'library_vite_abs_data'))

  ids <- reactiveValues(ids.223 = integrated.223, ids.431 = integrated.431, ids.451 = integrated.451)

  app = shinyApp(

    ui <- fluidPage(

  titlePanel('Peak Identification'),

  fluidRow(

    column(4,
      wellPanel(
        selectInput("wave", "Wavelength", choices = c("223" ,"431","451")),
        actionButton("ending", "Finished")
      )#end of wellPanel
      ),#end of column

    column(8,
           plotlyOutput(outputId = "plot1", width = '100%', height = '300px')
           )#end of 2nd column

      ),#end of fluid row 1

  hr(),
  fluidRow(

    column(6,
           wellPanel(
             sliderInput('fit.threshold',
                         'Fit Threshold',
                         min = 0,
                         max = 1,
                         value = 0.5,
                         step = 0.01)
           )#end of wellPanel
    ),#end of column

    column(6,
           sliderInput('fit.norm.threshold',
                       'Fit Norm Threshold',
                       min = 0,
                       max = 30,
                       value = 10,
                       step = 0.2, width = '100%')
    )#end of 2nd column

  ),#end of fluid row 2

  fluidRow(
    column(4,
           wellPanel(
             htmlOutput("peak"),
             verbatimTextOutput('currentID'),
             checkboxInput("checkbox", label = "Amend Peak ID?", value = FALSE),
             htmlOutput("checkbox1"),
             htmlOutput("conditional1"),
             htmlOutput("confirm.choice"),
             checkboxInput("add.trace", label = 'Add trace to plot?', value = FALSE),#this
             htmlOutput("conditional2") #this
           )#end of wellPanel
    ),#end of column

    column(8,
           plotlyOutput(outputId = "plot2", width = '100%', height = '400px')
    )#end of 2nd column
  )
  ),#end of UI call


# Define server logic for random distribution app ----
server <- function(input, output,session) {

  wave <<- reactive({
    switch(input$wave,
           "223" = dat$X223_final_corrected,
           "431" = dat$X431_final_corrected,
           "451" = dat$X451_final_corrected)
  })#end of wave

  this.int <<- reactive({
    switch(input$wave,
           "223" = ids$ids.223,
           "431" = ids$ids.431,
           "451" = ids$ids.451)
  })#end of peak

  this.abs <<- reactive({
    switch(input$wave,
           "223" = abs.spec.223.correct,
           "431" = abs.spec.431.correct,
           "451" = abs.spec.451.correct)
  })#end of this.abs

  this.lib <<- reactive({
    switch(input$wave,
           "223" = lib.vite.spec,
           "431" = lib.abs.spec,
           "451" = lib.abs.spec)
  })#end of this.lib

  this.peak <<- reactive({
    input$peak %>% substr(6, nchar(input$peak)) %>% as.numeric
  })

  #make a list to catch best pigment ID matches....
  my.matches <<- reactiveValues(  )

  output$peak <- renderUI({
    current.int <- this.int()
    selectInput("peak", 'Select Peak', choices = paste('Peak',1:nrow(current.int),sep=' '), selected=character(0))
  })

  plot2 <- reactive({
    abs.2 <<- this.abs()
    lib.2 <<- this.lib()
    peak.2 <<- this.peak()
    my.first <<- 'vitamin.e'
    my.second <<- 'vitamin.e'
    my.third <<- 'vitamin.e'
    my.matches$best <- c(my.first, my.second, my.third)

    plot_ly(x = abs.2$waves, y = abs.2[,peak.2+1] / max(abs.2[,peak.2+1]), type='scatter', mode='line', name = 'Sample') %>%
      add_trace(x = lib.2$wave, y = lib.2$vitamin.e_223_34.57, name = 'Best Match')
  })

  plot3 <- reactive({
    my.abs <<- this.abs()
    my.lib <<- this.lib()
    my.peak <<- this.peak()
    my.int <- ids$ids.431
    #get top three ID choices
    my.first <<- my.int$pigment_id[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.second <<- my.int$pigment_id2[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.third <<- my.int$pigment_id3[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.matches$best <- c(my.first, my.second, my.third)

    ref.1 <<- my.lib[,which(lib.meta$pigments_in_library == my.first)[1]+1]
    ref.2 <<- my.lib[,which(lib.meta$pigments_in_library == my.second)[1]+1]
    ref.3 <<- my.lib[,which(lib.meta$pigments_in_library == my.third)[1]+1]

    plot_ly(x = my.abs$waves, y = my.abs[,my.peak+1] / max(my.abs[,my.peak+1]), type='scatter', mode='line', name='sample') %>%

      add_trace(x = my.lib$wave, y = ref.1, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id[my.peak], '</br> Fit: ', round(my.int$fit[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'orange'), name = '1st best fit') %>%

      add_trace(x = my.lib$wave, y = ref.2, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id2[my.peak], '</br> Fit: ', round(my.int$fit2[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm2[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'green'), name = '2nd best fit') %>%

      add_trace(x = my.lib$wave, y = ref.3, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id3[my.peak], '</br> Fit: ', round(my.int$fit3[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm3[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'red'), name = '3rd best fit')
  })

  plot4 <- reactive({
    my.abs <<- this.abs()
    my.lib <<- this.lib()
    my.peak <<- this.peak()
    my.int <- ids$ids.451

    #get top three ID choices
    my.first <<- my.int$pigment_id[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.second <<- my.int$pigment_id2[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.third <<- my.int$pigment_id3[my.peak] %>% strsplit(.,split='_') %>% lapply('[[',1) %>% unlist
    my.matches$best <- c(my.first, my.second, my.third)

    ref.1 <<- my.lib[,which(lib.meta$pigments_in_library == my.first)[1]+1]
    ref.2 <<- my.lib[,which(lib.meta$pigments_in_library == my.second)[1]+1]
    ref.3 <<- my.lib[,which(lib.meta$pigments_in_library == my.third)[1]+1]

    plot_ly(x = my.abs$waves, y = my.abs[,my.peak+1] / max(my.abs[,my.peak+1]), type='scatter', mode='line', name = 'Sample') %>%

      add_trace( x = my.lib$wave, y = ref.1, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id[my.peak], '</br> Fit: ', round(my.int$fit[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'orange'), name = '1st best fit') %>%

      add_trace( x = my.lib$wave, y = ref.2, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id2[my.peak], '</br> Fit: ', round(my.int$fit2[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm2[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'green'), name = '2nd best fit') %>%

      add_trace( x = my.lib$wave, y = ref.3, text = ~paste('</br> Peak: ', my.peak, '</br> Pigment: ',my.int$pigment_id3[my.peak], '</br> Fit: ', round(my.int$fit3[my.peak],2), '</br> Fit Norm: ', round(my.int$fit_norm3[my.peak],2)), mode = 'line', hoverinfo = 'text', line = list(color = 'red'), name = '3rd best fit')
  })

  # Return the requested graph
  graphInput <- reactive({
    switch(input$wave,
           "223" = plot2(),
           "431" = plot3(),
           "451" = plot4()
    )
  })

  output$plot2 <- renderPlotly({
    if(input$add.trace == FALSE){
      graphInput()
      } else {
        ref.4 <<- my.lib[,which(lib.meta$pigments_in_library == input$trace)[1]+1]
        graphInput() %>% add_trace(x = my.lib$wave, y = ref.4, name = input$trace, mode='line', line=list(color = 'grey'))
      }
  })

  output$checkbox1 <- renderUI({
    xy <- input$checkbox
    if(length(xy) > 0 && xy == TRUE){
    radioButtons("checkbox1", label = 'Select Peak ID', choices = c(
      my.matches$best[1] ,
      my.matches$best[2] ,
      my.matches$best[3] ,
      'unknown' ,
      'other'),selected = character(0))
    }
  })

  output$conditional1 <- renderUI({
    if(input$checkbox == TRUE){
    if(length(input$checkbox1) > 0 && input$checkbox1 == "other"){
      selectInput('other_id', 'Other', choices = lib.meta$pigments_in_library %>% sort, selected = character(0))
    }
    }
  })

  output$conditional2 <- renderUI({
    if(input$add.trace == TRUE){
      selectInput('trace', 'Select Pigment', choices = lib.meta$pigments_in_library %>% sort, selected = character(0))
    }
  })

  output$confirm.choice <- renderUI({
    ff <- input$checkbox
    if(length(ff) > 0 && ff == TRUE){
      checkboxInput("confirm.choice", label = 'Confirm Peak Amendment', value = FALSE)
    }
  })

  #################################### Peak IDs#####################################
  #manualover-ride threshold derived IDs based on user input - Add a button to confirm?
  amend.id <<- reactiveValues()
  observe({
    check.id <<- input$checkbox1
    cond.id <<- input$other_id
    if(!is.null(check.id) && check.id == 'other'){
      amend.id$id <- cond.id
    } else {
      if(!is.null(check.id)){
        amend.id$id <- check.id
      }
    }
  })
  observeEvent(input$confirm.choice, {
    if(!is.null(input$confirm.choice) && input$confirm.choice == TRUE){
      tt <- this.peak()

      if(input$wave == '223'){
        ids$ids.223$final_id[tt] <- amend.id$id
        ids$ids.223$amended_id[tt] <- TRUE
      }
      if(input$wave == '431'){
        ids$ids.413$final_id[tt] <- amend.id$id
        ids$ids.431$amended_id[tt] <- TRUE
      }
      if(input$wave == '451'){
        ids$ids.451$final_id[tt] <- amend.id$id
        ids$ids.451$amended_id[tt] <- TRUE
      }
    }
  })

  #Take fit and fit.norm threshold inputs to select from pigment_id column
  my.fits <<- reactiveValues()
  observeEvent(input$fit.threshold,{
    my.fits$fit.threshold <<- input$fit.threshold
  })
  observeEvent(input$fit.norm.threshold,{
    my.fits$fit.norm.threshold <<- input$fit.norm.threshold
  })

  #calculate which rows of integration tables have IDs above thresholds
  my.223.rows <<- reactive({
    which(ids$ids.223$fit > my.fits$fit.threshold & ids$ids.223$fit_norm > my.fits$fit.norm.threshold)
  })
  my.431.rows <<- reactive({
    which(ids$ids.431$fit > my.fits$fit.threshold & ids$ids.431$fit_norm > my.fits$fit.norm.threshold)
  })
  my.451.rows <<- reactive({
    which(ids$ids.451$fit > my.fits$fit.threshold & ids$ids.451$fit_norm > my.fits$fit.norm.threshold)
  })

  #change the final_ID column of integration tables based on my.rows above
  observe({
    for(i in 1:nrow(ids$ids.223)){
      if(ids$ids.223$amended_id[i] == TRUE){
        ids$ids.223$final_id[i] <- ids$ids.223$final_id[i]
      } else {
        if(length(my.223.rows()) > 0 && i %in% my.223.rows()){
          ids$ids.223$final_id[i] <- ids$ids.223$pigment_id[i]
        } else {
          ids$ids.223$final_id[i] <- 'unknown'
        }
      }}

      for(i in 1:nrow(ids$ids.431)){
        i<<-i
        if(ids$ids.431$amended_id[i] == TRUE){
          ids$ids.431$final_id[i] <- ids$ids.431$final_id[i]
        } else {
          if(length(my.431.rows()) > 0 && i %in% my.431.rows()){
            ids$ids.431$final_id[i] <- ids$ids.431$pigment_id[i]
          } else {
            ids$ids.431$final_id[i] <- 'unknown'
          }
        }}

    for(i in 1:nrow(ids$ids.451)){
      i<<-i
      if(ids$ids.451$amended_id[i] == TRUE){
        next
        #ids$ids.451$final_id[i] <- ids$ids.451$final_id[i]
      } else {
        if(length(my.451.rows()) > 0 && i %in% my.451.rows()){
          ids$ids.451$final_id[i] <- ids$ids.451$pigment_id[i]
        } else {
          ids$ids.451$final_id[i] <- 'unknown'
        }
      }
    }#end of for loop
    updateCheckboxInput(session, 'checkbox', label = "Amend Peak ID?", value = FALSE )
  })

  #produce currentID for verbatum text
  output$currentID <- renderText({
    current.peak <- this.peak()
    current.int <- this.int()
    current.id <<- current.int$final_id[current.peak]
    paste('Current peak ID: ', current.id, sep='')
  })


  output$plot1 <- renderPlotly({
    current.int <- this.int()
    int.451 <- ids$ids.451
    if('chla' %in% int.451$final_id){
      chla.rt <- int.451$rt_real[int.451$final_id == 'chla'][1]
      time.dif <- ref$rent_time[ref$pigment == 'chla'] - chla.rt
      ref$rent_time <- ref$rent_time - time.dif
    }
    known <- subset(current.int, final_id!='unknown')
    known.nam <<- known$final_id
    if(length(known.nam)<1){known.nam <<-''} #this can definitely be imporved...
    unknown <- subset(current.int, final_id == 'unknown')
    p <- plot_ly(x = dat$time, y = wave(), type='scatter', mode='line') %>%
      add_annotations(x = current.int$rt_real, y = current.int$height, text = current.int$peak) %>%
         add_markers(x = known$rt_real, y = known$height, type='scatter', mode='marker', hoverinfo='text', text = known.nam) %>%
            add_markers(x = as.numeric(ref$rent_time), y = rep(0,nrow(ref)), xend = ref$rent_time, yend = rep(max(ids$ids.451$height)*1.2, nrow(ref)), hoverinfo='text', text = ref$pigment) %>%
      layout(showlegend = FALSE)
  })

  #reset inputs on peak / wave change
  observeEvent(input$peak,{
    updateCheckboxInput(session, 'checkbox', label = "Amend Peak ID?", value = FALSE )
    updateSelectInput(session, 'other_id', 'Other', choices = lib.meta$pigments_in_library %>% sort, selected = character(0))
    updateCheckboxInput(session, 'confirm.choice', label = 'Confirm Peak Amendment', value = FALSE)
    updateCheckboxInput(session, 'add.trace', label = 'Add trace to plot?', value = FALSE)
    updateSelectInput(session, 'trace', 'Select Pigment', choices = lib.meta$pigments_in_library %>% sort, selected = character(0))
  })
  observeEvent(input$wave,{
    updateCheckboxInput(session, 'checkbox', label = "Amend Peak ID?", value = FALSE )
    updateCheckboxInput(session, 'confirm.choice', label = 'Confirm Peak Amendment', value = FALSE)
    updateCheckboxInput(session, 'add.trace', label = 'Add trace to plot?', value = FALSE)
    updateSelectInput(session, 'trace', 'Select Pigment', choices = lib.meta$pigments_in_library %>% sort, selected = character(0))
  })

  observeEvent(input$ending, {
    stopApp()
  })
  session$onSessionEnded(function() {
    stopApp()
  })
}

)#end of shinAPP call
runApp(app)

#get back info
integrated.223 <- isolate(ids$ids.223)
integrated.431 <- isolate(ids$ids.431)
integrated.451 <- isolate(ids$ids.451)
copy_to(dest = mydb, df = integrated.223, name = "integrated.223", overwrite = T, temporary=F)
copy_to(dest = mydb, df = integrated.431, name = "integrated.431", overwrite = T, temporary=F)
copy_to(dest = mydb, df = integrated.451, name = "integrated.451", overwrite = T, temporary=F)


#write out meta.data table for this sample
meta.info<-c('Sample_Name','Peaks_Integrated_Per_Wavelength', 'Peaks_IDd_Per_Wavelength', 'Processing_Settings')
processing.settings<-paste('nups.451 = ',nups.451, ';nups.223 = ', nups.223,';min.peak.height = ', min.height, ':fit.threshold = ', isolate(my.fits$fit.threshold),';fit.norm.threshold = ', isolate(my.fits$fit.norm.threshold), sep='')
meta.data<-data.frame(attributes=meta.info, values=c(db.nam,paste(nrow(integrated.223), nrow(integrated.431), nrow(integrated.451), sep=','),
paste(length(which(grepl(pattern = 'unknown',x = integrated.223$final_id)==F)), length(which(grepl(pattern = 'unknown',x = integrated.431$final_id)==F)),length(which(grepl(pattern = 'unknown',x = integrated.451$final_id)==F)),sep=','),processing.settings))
copy_to(dest = mydb, df = meta.data, name = "meta.data", overwrite = T, temporary = F)

###############################################################################################
print('all finished')

}#end of identify.peaks function

#things to complete:
# new tab showing integration tables
# option to plot new library spectrum on lower graph
