#'@title integarte.peaks
#'@description identify peaks in output of EXPORT3D.mac for ChemStation hplc data
#'@param my.data character string name or filepath and name of the CSV containing the data
#'@param blank character string name or filepath and name of a single representative blank sample CSV
#'@param nups.451 number of ups to define peak at 451 nm (can be altered if peaks not captured correctly)
#'@param nups.223 number of ups to define peak at 223 nm (can be altered if peaks not captured correctly)
#'@return SQLite database containing the following tables:
#'    1. meta.data: meta data table containing sample information, including information on extraction and filtration volumes and cell counts associated with the samples.
#'    2. vite_blank_data: containing associated blank samples vitamin E concentrations
#'    3. mobile_phase_reference_spectra: full spectra for blank samples used in blanking process
#'    4. library_meta_data: meta data table for reference library showing pigment names, wavelength observed, retention times
#'    5. library_abs_data: full spectra for all pigments in the library
#'    6. library_vite_abs_data: spectra for vitamin E internal standard library
#'    7. raw_abs_data: raw EXPORT3D.mac output from ChemStation
#'    8. corrected_abs_data: 223 nm, 431 nm, and 451 nm absorbance data, corrected for background (mobile phase reference) and 751 nm background absorbance
#'    9.etc
#'@export

integrate.peaks <- function(my.data='', blank = NULL, nups.451=10, nups.223=20){

  #1. load libraries
  library(alsace)
  library(Peaks)
  library.dynam('Peaks', 'Peaks', lib.loc=NULL)
  library('dplyr')
  library('dbplyr')
  library('RSQLite')
  library('pracma')
  library('shiny')
  library('plotly')
  #source('./R/integrate.peaks.R')
  #source('./R/identify.peaks.R')
  source('./R/mobile_reference_construction.R')

  #load in data and organise columns, add a 'time' column
  df1 <- read.csv(my.data, skipNul = T, header=F, sep=',', fileEncoding="UCS-2LE")
  df1[1,1] <- 'time'
  colnames(df1) <- paste('X',df1[1,],sep='')
  colnames(df1)[1] <- 'time'
  df1 <- df1[2:nrow(df1),]
  df1$time <- as.numeric(df1$time)


  #3. create an SQLite data base to store all the outputs of the script below (nice and easy to read back in). This will grepl a name from the database name. It will check first to see if database is in another directory (and thus needs file path removed)
  if(my.data %in% dir()){
    db.nam <- my.data %>% strsplit(.,split='[.]') %>% lapply('[[',1)  %>% unlist
  } else {
    db.nam <- my.data %>% basename %>% strsplit(.,split='[.]') %>% lapply('[[',1)  %>% unlist
  }
  mydb <- src_sqlite(db.nam, create = T)
  print(paste('Creating database *',db.nam,'* in ',getwd(), sep=''))

  #add in raw_abs data to database
  copy_to(mydb, df1, temporary=F, name='raw_abs_data', overwrite=T)

  #extract time and wavelength info
  this.time <- df1$time
  wave <- as.numeric(gsub(pattern = 'X',replacement = '',x = colnames(df1)[2:ncol(df1)]))


  #subset out my wavelengths of interest for chromatograms
  my.waves <- c(222,430,450)+1
  my.ref.wave <- df1$X751
  my.cols <- paste('X',my.waves,sep='')
  dat <- data.frame(time = df1[,'time'])

  #deduct the reference 751nm from each of my wave lengths of interest
  for(i in 1:length(my.cols)){
    dat[, my.cols[i]] <- df1[,my.cols[i]] - my.ref.wave
  }

  #read in mobile phase reference, add to my database and blank correct wavelengths of interest.
  if(is.null(blank)){
    blank <- blank_r1_example #read in example if blank not provided by user
  } else {
    blank <- read.csv(blank, skipNul = T, header=F, sep=',', fileEncoding="UCS-2LE")
  }

  #run mobile_reference_construction function on blank data
  blank <- mobile_reference_construction(blank)
  #fix to catch difference (1 row) between blank dataset and sample data sets
  blank <- blank[1:nrow(dat), ]

  #add blank to main database
  copy_to(mydb, name='mobile_phase_reference_spectra', df = blank, temporary=F, overwrite=T)

  #subtract the blank from each wavelength
  dat$X223_mobile_ref_corrected <- dat$X223 - blank$X223
  dat$X431_mobile_ref_corrected <- dat$X431 - blank$X431
  dat$X451_mobile_ref_corrected <- dat$X451 - blank$X451


  #calculate 'spectrum' background - given blank deduction above, this does very little
  dat$X223_final_corrected <- dat$X223_mobile_ref_corrected -
    SpectrumBackground(dat$X223_mobile_ref_corrected)
  dat$X431_final_corrected <- dat$X431_mobile_ref_corrected -
    SpectrumBackground(dat$X431_mobile_ref_corrected)
  dat$X451_final_corrected <- dat$X451_mobile_ref_corrected -
    SpectrumBackground(dat$X451_mobile_ref_corrected)

  #add dat into database
  copy_to(mydb, dat, temporary=F, name='corrected_abs_data', overwrite=T)


  #From here we need to be interactive to minimum peak height.......

  #identify peaks in each chromatograph.
  peaks.223 <- pracma::findpeaks(dat$X223_final_corrected, nups = nups.223) %>% as.data.frame
  peaks.431 <- pracma::findpeaks(dat$X431_final_corrected, nups = nups.451) %>% as.data.frame
  peaks.451 <- pracma::findpeaks(dat$X451_final_corrected, nups = nups.451) %>% as.data.frame
  colnames(peaks.223) <- c('height','max_loc','start','end')
  colnames(peaks.431) <- c('height','max_loc','start','end')
  colnames(peaks.451) <- c('height','max_loc','start','end')

  #re-track peak position to retention time
  peaks.223$rt_real <- dat$time[peaks.223[,2]]
  peaks.431$rt_real <- dat$time[peaks.431[,2]]
  peaks.451$rt_real <- dat$time[peaks.451[,2]]

  #get peak information (height / area) using alsace package, order the table, add peak number, write tables to database
  #get peak info
  integrated.223 <- alsace::fitpeaks(dat$X223_final_corrected, peaks.223[,2]) %>% as.data.frame
  integrated.431 <- alsace::fitpeaks(dat$X431_final_corrected, peaks.431[,2]) %>% as.data.frame
  integrated.451 <- alsace::fitpeaks(dat$X451_final_corrected, peaks.451[,2]) %>% as.data.frame

  #add in real time retention times
  integrated.223$rt_real <- peaks.223$rt_real
  integrated.431$rt_real <- peaks.431$rt_real #throws the error if no peaks found, so drop number to something stupid to get a least one peak!
  integrated.451$rt_real <- peaks.451$rt_real

  #order according to rt_real
  integrated.223 <- integrated.223[order(integrated.223$rt_real, decreasing = F),]
  integrated.431 <- integrated.431[order(integrated.431$rt_real, decreasing = F),]
  integrated.451 <- integrated.451[order(integrated.451$rt_real, decreasing = F),]

  tmp.223 <- integrated.223
  tmp.431 <- integrated.431
  tmp.451 <- integrated.451


  #start of shiny app section
  app = shinyApp(
    ui <- fluidPage(
      titlePanel("Select minimum peak height"),

      fluidRow(

        column(12,
               wellPanel(
                 sliderInput('min.451.height',
                             'Min 451 nm Height',
                             min = 0,
                             max = round(max(dat$X451_final_corrected)*0.2,0),
                             value = 10,
                             step = 0.1)
               )#end of wellPanel
        )#end of column
      ),#end of fluid row 1

      fluidRow(
        column(12,
               plotlyOutput('plot', height = "500"))
      ),

      fluidRow(
        column(4,
               actionButton("ending", "Finished"))
      )
     ),

    server <- function(input, output,session) {

      #use inputs from server to define new limits for graph...
      #subset based on the 'height' column in the integration tables (set to -Inf in main function...)
      output$plot <- renderPlotly({

        tmp.451 <- integrated.451[integrated.451$height > input$min.451.height,]
        tmp.451$peak <- 1:nrow(tmp.451)
        min.height <<- input$min.451.height
        plot_ly(x = dat$time, y = dat$X451_final_corrected, type='scatter', mode = 'line') %>%
          add_annotations(x = tmp.451$rt_real, y = tmp.451$height, text = tmp.451$peak, font = list(color='orange'))
      })

      observeEvent(input$ending,{
        stopApp()
      })

      session$onSessionEnded(function(){
        stopApp()
      })

    }
  )
  runApp(app)

  #should now have min.height value to use for subsetting integration tables...
  integrated.223 <- integrated.223[integrated.223$height > min.height, ]
  integrated.431 <- integrated.431[integrated.431$height > min.height, ]
  integrated.451 <- integrated.451[integrated.451$height > min.height, ]

  #add peak number column
  integrated.223$peak <- 1:nrow(integrated.223)
  integrated.431$peak <- 1:nrow(integrated.431)
  integrated.451$peak <- 1:nrow(integrated.451)


  #make an abs spectra for each peak integrated per wavelength processed from 300 - 700 nm, write to data frames, identify peaks in these for printing to abs spec plots, write to database, add info to meta.data

  #223 data
  abs.spec.223 <- data.frame(waves=wave)
  abs.spec.223.correct <- data.frame(waves=wave)

  for(i in 1:nrow(integrated.223)){
    peak.no <- integrated.223$peak[i]
    rt <- integrated.223$rt_real[i]
    abs.spec.223[,i+1] <- df1[which(df1$time == rt), 2:ncol(df1)] %>% as.numeric #causing a problem here!!!
    colnames(abs.spec.223)[i+1] <- paste('peak_',peak.no,sep='')
    blank.row <- which(abs(blank$time - rt) == min(abs(blank$time -rt)))
    if(length(blank.row) > 1){
      blank.row <- blank.row[1]
    }
    blank.dat <- blank[blank.row, 2:ncol(blank)] %>% as.numeric
    abs.spec.223.correct[,i+1] <- as.numeric(df1[which(df1$time==rt), 2:ncol(df1)]) - as.numeric(blank.dat)
    colnames(abs.spec.223.correct)[i+1] <- paste('peak_corrected_', peak.no, sep='')
  }

  #write to db
  copy_to(mydb, abs.spec.223, temporary=F, name='223.peaks.abs.spec', overwrite=T)
  copy_to(mydb, abs.spec.223.correct, temporary=F, name='223.peaks.abs.spec.correct', overwrite=T)

  #integrate main peaks in abs spec (TOTAL Range for 223 (as is UV signal)) and add into integrated.223 table
  for(i in 2:ncol(abs.spec.223.correct)){
    x <- SpectrumSearch(abs.spec.223.correct[,i], sigma=3.0, background=F, markov=F,threshold = 40)$pos
    integrated.223$abs_max[i-1] <- paste0(abs.spec.223.correct$waves[x], collapse=',')
  }

  #431 data
  abs.spec.431 <- data.frame(waves=wave)
  abs.spec.431.correct <- data.frame(waves=wave)
  for(i in 1:nrow(integrated.431)){
    peak.no <- integrated.431$peak[i]
    rt <- integrated.431$rt_real[i]
    abs.spec.431[,i+1] <- as.numeric(df1[which(df1$time==rt), 2:ncol(df1)])
    colnames(abs.spec.431)[i+1] <- paste('peak_',peak.no,sep='')

    blank.row <- which(abs(blank$time-rt)==min(abs(blank$time-rt)))
    if(length(blank.row) > 1){
      blank.row <- blank.row[1]
    }
    blank.dat <- blank[blank.row, 2:ncol(blank)]
    abs.spec.431.correct[,i+1] <- as.numeric(df1[which(df1$time==rt), 2:ncol(df1)]) - as.numeric(blank.dat)
    colnames(abs.spec.431.correct)[i+1] <- paste('peak_corrected_', peak.no, sep='')
  }

  #cut down to 350 - 700 region
  abs.spec.431 <- subset(abs.spec.431, waves > 349 & waves < 701)
  abs.spec.431.correct <- subset(abs.spec.431.correct, waves > 349 & waves < 701)
  copy_to(mydb, abs.spec.431, temporary=F, name='431.peaks.abs.spec', overwrite=T)
  copy_to(mydb, abs.spec.431.correct, temporary=F, name='431.peaks.abs.spec.correct', overwrite=T)

  #add in abs spec integration for 431 into integrated.431
  for(i in 2:ncol(abs.spec.431.correct)){
    x <- SpectrumSearch(abs.spec.431.correct[,i], sigma=3.0, background=F, markov=F,threshold = 40)$pos
    integrated.431$abs_max[i-1] <- paste0(abs.spec.431.correct$waves[x], collapse=',')
  }

  #451 data
  abs.spec.451 <- data.frame(waves=wave)
  abs.spec.451.correct <- data.frame(waves=wave)
  for(i in 1:nrow(integrated.451)){
    peak.no <- integrated.451$peak[i]
    rt <- integrated.451$rt_real[i]
    abs.spec.451[,i+1] <- as.numeric(df1[which(df1$time==rt), 2:ncol(df1)])
    colnames(abs.spec.451)[i+1] <- paste('peak_',peak.no,sep='')

    blank.row <- which(abs(blank$time-rt)==min(abs(blank$time-rt)))
    if(length(blank.row) > 1){
      blank.row <- blank.row[1]
    }
    blank.dat <- blank[blank.row, 2:ncol(blank)]
    abs.spec.451.correct[,i+1] <- as.numeric(df1[which(df1$time==rt), 2:ncol(df1)]) - as.numeric(blank.dat)
    colnames(abs.spec.451.correct)[i+1] <- paste('peak_corrected_', peak.no, sep='')
  }

  #cut down to 300 - 700 region
  abs.spec.451 <- subset(abs.spec.451, waves > 349 & waves < 701)
  abs.spec.451.correct <- subset(abs.spec.451.correct, waves > 349 & waves < 701)
  copy_to(mydb, abs.spec.451, temporary=F, name='451.peaks.abs.spec', overwrite=T)
  copy_to(mydb, abs.spec.451.correct, temporary=F, name='451.peaks.abs.spec.correct', overwrite=T)

  #add in abs spec integration for 451 into integrated.451 table
  for(i in 2:ncol(abs.spec.451.correct)){
    x <- SpectrumSearch(abs.spec.451.correct[,i], sigma=3.0, background=F, markov=F,threshold = 40)$pos
    integrated.451$abs_max[i-1] <- paste0(abs.spec.451.correct$waves[x], collapse=',')
  }


  ###############################################################
  #Draw data from library that is incorporated into the package - loaded in package
  # lib.meta <- as.data.frame(tbl(library, "library_meta_data"))
  # lib.abs.spec<-as.data.frame(tbl(library, "library_abs_data"))
  # lib.abs.spec<-subset(lib.abs.spec, wave > 349)
  # lib.vite.spec<-as.data.frame(tbl(library, "vite_abs_data"))

  #add these into sample database
  copy_to(mydb, name='library_meta_data', df = lib.meta, temporary=F, overwrite=T)
  copy_to(mydb, name='library_abs_data', df = lib.abs.spec, temporary=F, overwrite=T)
  copy_to(mydb, name='library_vite_abs_data', df = lib.vite.spec, temporary=F, overwrite=T)

  ###############################################################

  #Peak Identification
  #for each peak, take, normalise, compare, compute statistics, write outputs to integrated.X df, then re-write to mydb.

  #223 data - only comparing to lib.vite.spec which has only x1 spec in it anyway!
  integrated.223$pigment_id <- NA
  integrated.223$fit <- NA
  integrated.223$time_diff <- NA
  integrated.223$fit_norm <- NA

  for(i in 1:nrow(integrated.223)){
    #get peak name
    peak <- integrated.223$peak[i]
    peak.nam <- paste('peak_corrected',peak,sep='_')

    #get peak rt info
    rt <- integrated.223$rt_real[i]

    #get out of abs.spec table and normalise to max value
    y <- abs.spec.223.correct[,peak.nam]
    y.norm <- y/max(y)

    if(is.nan(y.norm[1])){
      y.norm <- y
    }

    #compare to vite library
    x <- lib.vite.spec[,2]
    my.sum <- summary(lm(y ~ x))$adj.r.squared
    if(is.nan(my.sum)){
      integrated.223$fit[i] <- 0
    } else {
      integrated.223$fit[i] <- my.sum
    }
    rm(my.sum)
    #compare retention times
    integrated.223$time_diff[i] <- rt - lib.meta$retention_times[lib.meta$pigments_in_library=='vitamin.e']

    #add fit norm value
    integrated.223$fit_norm[i] <- integrated.223$fit[i]*(1/abs(integrated.223$time_diff[i]))

    #add in vite name
    integrated.223$pigment_id[i] <- 'vitamin.e'

  }#end of i loop

  #431 data
  integrated.431$pigment_id<-NA
  integrated.431$fit<-NA
  integrated.431$time_diff<-NA
  integrated.431$fit_norm<-NA
  integrated.431$pigment_id2<-NA
  integrated.431$fit2<-NA
  integrated.431$time_diff2<-NA
  integrated.431$fit_norm2<-NA
  integrated.431$pigment_id3<-NA
  integrated.431$fit3<-NA
  integrated.431$time_diff3<-NA
  integrated.431$fit_norm3<-NA

  for(i in 1:nrow(integrated.431)){
    #get peak name
    peak<-integrated.431$peak[i]
    peak.nam<-paste('peak_corrected',peak,sep='_')

    #get peak rt info
    rt<-integrated.431$rt_real[i]

    #get out of abs.spec table and normalise to max value
    y<-abs.spec.431.correct[,peak.nam]
    y.norm<-y/max(y)
    if(is.nan(y.norm[1])){
      y.norm<-y
    }

    #get library and loop through library entries and compare
    active.lib.meta<-lib.meta[lib.meta$wave_length_observed==451,] #add fits here to start with...
    active.lib.meta$fit<-NA
    active.lib.meta$time_diff<-NA
    active.lib.meta$fit_norm<-NA
    for(j in 1:nrow(active.lib.meta)){
      x<-lib.abs.spec[,j+1]
      my.info<-summary(lm(y.norm ~ x))$adj.r.squared #changed this to y.norm as I think it should be!!
      if(is.nan(my.info)){
        active.lib.meta$fit[j]<-0
      } else {
        active.lib.meta$fit[j]<-my.info
      }
      rm(my.info)
    }#end of j loop

    #compare retention times
    active.lib.meta$time_diff<-rt-active.lib.meta$retention_times

    #Compute fit_norm, i.e. fit R2 multiplied by 1/time_diff - NB if time_diff = 0 1/0 will give infinity
    active.lib.meta$fit_norm<-active.lib.meta$fit*(1/abs(active.lib.meta$time_diff))

    #order fit results based on fit column
    active.lib.meta<-active.lib.meta[order(active.lib.meta$fit, decreasing=T),]
    #then re-order the top three rows of this based on fit_norm
    active.lib.meta<-active.lib.meta[1:5,]
    active.lib.meta<-active.lib.meta[order(active.lib.meta$fit_norm, decreasing=T),]

    #add outcomes to integrated.451
    integrated.431$pigment_id[i]<-active.lib.meta$pigments_in_library[1]
    integrated.431$fit[i]<-active.lib.meta$fit[1]
    integrated.431$time_diff[i]<-active.lib.meta$time_diff[1]
    integrated.431$fit_norm[i]<-active.lib.meta$fit_norm[1]

    integrated.431$pigment_id2[i]<-active.lib.meta$pigments_in_library[2]
    integrated.431$fit2[i]<-active.lib.meta$fit[2]
    integrated.431$time_diff2[i]<-active.lib.meta$time_diff[2]
    integrated.431$fit_norm2[i]<-active.lib.meta$fit_norm[2]

    integrated.431$pigment_id3[i]<-active.lib.meta$pigments_in_library[3]
    integrated.431$fit3[i]<-active.lib.meta$fit[3]
    integrated.431$time_diff3[i]<-active.lib.meta$time_diff[3]
    integrated.431$fit_norm3[i]<-active.lib.meta$fit_norm[3]

  }#end of i loop


  #451 data
  integrated.451$pigment_id<-NA
  integrated.451$fit<-NA
  integrated.451$time_diff<-NA
  integrated.451$fit_norm<-NA
  integrated.451$pigment_id2<-NA
  integrated.451$fit2<-NA
  integrated.451$time_diff2<-NA
  integrated.451$fit_norm2<-NA
  integrated.451$pigment_id3<-NA
  integrated.451$fit3<-NA
  integrated.451$time_diff3<-NA
  integrated.451$fit_norm3<-NA

  for(i in 1:nrow(integrated.451)){
    #get peak name
    peak<-integrated.451$peak[i]
    peak.nam<-paste('peak_corrected',peak,sep='_')

    #get peak rt info
    rt<-integrated.451$rt_real[i]

    #get out of abs.spec table and normalise to max value
    y<-abs.spec.451.correct[,peak.nam]
    y.norm<-y/max(y)

    if(is.nan(y.norm[1])){
      y.norm<-y
    }

    #get library and loop through library entries and compare
    active.lib.meta<-lib.meta[lib.meta$wave_length_observed==451,] #add fits here to start with...
    active.lib.meta$fit<-NA
    active.lib.meta$time_diff<-NA
    active.lib.meta$fit_norm<-NA
    for(j in 1:nrow(active.lib.meta)){
      x<-lib.abs.spec[,j+1]
      my.info<-summary(lm(y.norm ~ x))$adj.r.squared #changed this to y.norm as I think it should be!!
      if(is.nan(my.info)){
        active.lib.meta$fit[j]<-0
      } else {
        active.lib.meta$fit[j]<-my.info
      }
      rm(my.info)
    }#end of j loop

    #compare retention times
    active.lib.meta$time_diff<-rt-active.lib.meta$retention_times

    #Compute fit_norm, i.e. fit R2 multiplied by 1/time_diff - NB if time_diff = 0 1/0 will give infinity
    active.lib.meta$fit_norm<-active.lib.meta$fit*(1/abs(active.lib.meta$time_diff))


    #order fit results based on fit
    active.lib.meta<-active.lib.meta[order(active.lib.meta$fit, decreasing=T),]
    #then re-order the top three rows of this based on fit_norm
    active.lib.meta<-active.lib.meta[1:5,]
    active.lib.meta<-active.lib.meta[order(active.lib.meta$fit_norm, decreasing=T),]

    #add outcomes to integrated.451
    integrated.451$pigment_id[i]<-paste(active.lib.meta$pigments_in_library[1])
    integrated.451$fit[i]<-active.lib.meta$fit[1]
    integrated.451$time_diff[i]<-active.lib.meta$time_diff[1]
    integrated.451$fit_norm[i]<-active.lib.meta$fit_norm[1]

    integrated.451$pigment_id2[i]<-paste(active.lib.meta$pigments_in_library[2])
    integrated.451$fit2[i]<-active.lib.meta$fit[2]
    integrated.451$time_diff2[i]<-active.lib.meta$time_diff[2]
    integrated.451$fit_norm2[i]<-active.lib.meta$fit_norm[2]

    integrated.451$pigment_id3[i]<-paste(active.lib.meta$pigments_in_library[3])
    integrated.451$fit3[i]<-active.lib.meta$fit[3]
    integrated.451$time_diff3[i]<-active.lib.meta$time_diff[3]
    integrated.451$fit_norm3[i]<-active.lib.meta$fit_norm[3]

  }#end of i loop



  #Clean up integration tables based on user defined fit.threshold and fit.norm.threshold
  integrated.223$final_id <- 'unknown'
  integrated.431$final_id <- 'unknown'
  integrated.451$final_id <- 'unknown'
  integrated.223$amended_id <- F
  integrated.431$amended_id <- F
  integrated.451$amended_id <- F

  copy_to(mydb, integrated.223, temporary = F, name = 'integrated.223', overwrite=T)
  copy_to(mydb, integrated.431, temporary = F, name = 'integrated.431', overwrite=T)
  copy_to(mydb, integrated.451, temporary = F, name = 'integrated.451', overwrite=T)
  copy_to(mydb, abs.spec.223.correct, temporary=F, name='223.peaks.abs.spec.correct', overwrite=T)
  copy_to(mydb, abs.spec.431.correct, temporary=F, name='431.peaks.abs.spec.correct', overwrite=T)
  copy_to(mydb, abs.spec.451.correct, temporary=F, name='451.peaks.abs.spec.correct', overwrite=T)

}
