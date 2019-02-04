####################################
#Create mobile phases reference for use in hplc_processing function using blank sample information

mobile_reference_construction<-function(blank){

  library('dplyr')
  library('dbplyr')
  library('RSQLite')

  df1 <- blank
  df1[1,1] <- 'time'
  colnames(df1) <- paste('X',df1[1,],sep='')
  colnames(df1)[1] <- 'time'
  df1 <- df1[2:nrow(df1),]
  df1$time <- as.numeric(df1$time)

  #take blank, section out 34 - 35 mins when vite E peak comes through, and use regression between 1 min before and 1 min after data to fill in the gap for all wavelengths, then put into library database and use to correct both library and sample data...
  blank <- df1
  blank.correct.1 <- df1[df1$time > 32 & df1$time < 33,]
  blank.correct.2 <- df1[df1$time > 35 & df1$time < 36,]
  blank.correct <- rbind(blank.correct.1, blank.correct.2) #use this to predict the period between 34 and 35 minutes for each wavelength

  new.time <- df1$time[df1$time > 32 & df1$time < 36]
  modelled.dat <- data.frame(time=new.time)
  modelled.dat[,colnames(df1)[2:ncol(df1)]] <- NA

  #loop through wavelengths, model the gap and split data into modelled.dat
  for(i in 2:ncol(blank.correct)){
    x <- blank.correct$time
    y <- blank.correct[,i]
    m1 <- lm(y~x)
    new.dat <- data.frame(x=new.time)
    my.pred <- as.numeric(predict(object=m1,newdata = new.dat ))
    modelled.dat[,i] <- my.pred
  }
  #cut down modelled.dat to only > 33 and < 35 mins
  modelled.dat <- modelled.dat[modelled.dat$time > 33 & modelled.dat$time < 35,]

  #paste into blank in right time rows..
  blank[which(blank$time %in% modelled.dat$time),] <- modelled.dat

  #put into own data base
  return(blank)

}#end of function

