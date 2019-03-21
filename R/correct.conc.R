#'@title correct.conc
#'@description calculate pigment concetrations using std coeffients
#'@param sample.database character string name or filepath and name of DATABASE name (or path to database if not in working directory) output from id.peaks function
#'@param std.coefs character string name or filepath and name of a Standard Curve RSQLite database
#'@param blank.vite value of the AVERAGE blank vitamin e peak area (need to integrate peaks for your blank samples and use view.database function to manually get and make average)
#'@param extraction.vol.ml volume of the extraction solvent used for pigment extractions - in MILLILITRES
#'@param filtration.vol.l volume of the amout of sample filtered for extraction - in LITRES
#'@return 'correct.conc' table within main sample database with concentrations for pigments a standard is available for in ug per litre, or the ratio of pigment area to chlorophyll area for those pigments a standard is not available for.
#'@export


correct.conc <- function(sample.database = NULL, std.coefs = NULL, blank.vite = NULL, extraction.vol.ml = NULL, filtration.vol.l = NULL){

  library('dplyr')
  library('dbplyr')
  library('RSQLite')


  if(is.null(sample.database) | is.null(std.coefs) | is.null(blank.vite)| is.null(extraction.vol.ml)| is.null(filtration.vol.l)){

    print('Not all requried information input')

  } else {

  #read in std coefs from RSQLite database
  stds <- src_sqlite(std.coefs,create = F)
  std.coefs <- as.data.frame(tbl(stds, 'all_coefficients'))

  #read in sample database
  db1 <- src_sqlite(sample.database,create = F)

  #load int.223 and get sample vite conc
  int.223 <- as.data.frame(tbl(db1, 'integrated.223'))
  samp.vite <- int.223$area[int.223$final_id=='vitamin.e']

  #take in blank vite average area from user input
  blank.vite <- blank.vite

  #load in extraction volume
  ex.vol <- extraction.vol.ml

  #load in filtration volume
  fil.vol <- filtration.vol.l

  #load int.431 and get 431 chla area (if exists)
  int.431 <- as.data.frame(tbl(db1, 'integrated.431'))
  row.431 <- which(int.431$final_id == 'chla')
  if(length(row.431) > 0 ){
    chla.area <- sum(int.431$area[row.431]) #incase of multiple entries
    chla.time <- int.431$rt_real[row.431][1] # do not input into int451, input into correct.conc only!
  } else {
    chla.area <- NA
    chla.time <- NA
  }

  #catch for if we have no pigments in the file
  if(is.na(chla.area)){
    print(paste(files[i], ':no pigments present'))
  } else {
    #aggregate int.431 to make sure chla data gets put into 1 row
    dat.431 <- aggregate(area ~ final_id, data=int.431, FUN=sum)

    #load in int.451
    int.451 <- as.data.frame(tbl(db1, 'integrated.451'))

    #make correct.conc data frame for calculations
    correct.conc <- data.frame(final_id=int.451$final_id,area=int.451$area)

    #sum incase we have multiple entries of the same pigment
    correct.conc <- aggregate(area ~ final_id, data=correct.conc, FUN=sum)

    #switch chla (or not) for dat.431 chla
    if('chla' %in% correct.conc$final_id){
      correct.conc$area[correct.conc$final_id=='chla'] <- dat.431$area[dat.431$final_id=='chla']
    } else {
      correct.conc[nrow(correct.conc)+1,] <- c('chla', dat.431$area[dat.431$final_id=='chla'])
    }

    #for those samples we have a standard for, calculate actual concentration, for those we dont, input same conc in this column, i.e. area * 1
    correct.conc$std.area <- NA
    for(i in 1:nrow(correct.conc)){
      pig <- correct.conc$final_id[i]
      if(pig %in% std.coefs$pigment){
        m <- std.coefs$m[std.coefs$pigment==pig]
        c <- std.coefs$c[std.coefs$pigment==pig]
        correct.conc$std.area[i] <- ((m) * as.numeric(correct.conc$area[i])) + (c)
      } else {
        correct.conc$std.area[i] <- correct.conc$area[i]
      }
    }


    #correct std conc for vite (blank vite / sample vite) and extraction / filtration volumes
    correct.conc$std.area <- as.numeric(correct.conc$std.area)
    correct.conc$area <- as.numeric(correct.conc$area)
    correct.conc$final.conc.ug.l <- correct.conc$std.area * (blank.vite / samp.vite) * (ex.vol / fil.vol)

    #make ratios to chla for those pigments with actual concentrations, and non.true ratios for those without, based on conc and/or area
    correct.conc$true.ratio <- NA
    correct.conc$area.ratio <- NA
    chla.conc <- correct.conc$final.conc.ug.l[which(grepl('chla',correct.conc$final_id))]
    chla.area <- correct.conc$area[which(grepl('chla',correct.conc$final_id))]* (blank.vite / samp.vite)*(ex.vol / fil.vol) #this is chla area with only extract/filter vols corrections applied
    #loop through and make ratios
    for(i in 1:nrow(correct.conc)){
      if(correct.conc$area[i]==correct.conc$std.area[i]){
        correct.conc$area.ratio[i]<-correct.conc$final.conc.ug.l[i]/chla.area
      }  else {
        correct.conc$true.ratio[i]<-correct.conc$final.conc.ug.l[i]/chla.conc
      }
    }

    #copy info into original database
    copy_to(db1, correct.conc, name='correct.conc', temporary=F, overwrite=T)
    copy_to(db1, std.coefs, name='standard.coefficients', temporary=F, overwrite=T)

    met <- as.data.frame(tbl(db1, 'meta.data'))
    print(paste('Concentrations calculate for',met$values[met$attributes=='Sample_Name']))
    return(correct.conc)
  }#end of input else catch
  }
}
