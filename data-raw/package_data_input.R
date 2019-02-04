#data input to make available to ProcessHplc package
library(RSQLite)
library(magrittr)

#example blank data
blank_r1_example <- read.csv('data-raw/blank_r1_example.csv', skipNul = T, header=F, sep=',', fileEncoding="UCS-2LE")
#blank_r2_example <- read.csv('data-raw/blank_r2_example.csv', skipNul = T, header=F, sep=',', fileEncoding="UCS-2LE")
#blank_r3_example <- read.csv('data-raw/blank_r3_example.csv', skipNul = T, header=F, sep=',', fileEncoding="UCS-2LE")

devtools::use_data(blank_r1_example, overwrite = T)
#devtools::use_data(blank_r2_example, overwrite = T)
#devtools::use_data(blank_r3_example, overwrite = T)

#example mobile phase reference
mobile_phase_reference_example <- src_sqlite('data-raw/mobile_phase_reference_example', create=F) %>% tbl('mobile_phase_reference_spectra') %>% as.data.frame
devtools::use_data(mobile_phase_reference_example, overwrite = T)

#pigment library
library <- src_sqlite('data-raw/library', create=F)
lib.meta <- as.data.frame(tbl(library, "library_meta_data"))
lib.abs.spec<-as.data.frame(tbl(library, "library_abs_data"))
lib.abs.spec<-subset(lib.abs.spec, wave > 349)
lib.vite.spec<-as.data.frame(tbl(library, "vite_abs_data"))
devtools::use_data(library, overwrite = T)
devtools::use_data(lib.meta, overwrite = T)
devtools::use_data(lib.abs.spec, overwrite = T)
devtools::use_data(lib.vite.spec, overwrite = T)


#pigment reference retention time lines...
ref <- read.csv('data-raw/pigment_retention_times.csv')
devtools::use_data(ref, overwrite=T)
