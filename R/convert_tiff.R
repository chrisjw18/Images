#convert images from png to TIFF

library(magick)
library(magrittr)
library(utils)

my.dir <- '~/Dropbox/2018 Greenland Samples/260718_surface_ICE_samples_26-50/'

my.folders <- dir(my.dir)
for(j in 1:length(my.folders)){
  my.files <- list.files(paste(my.dir, my.folders[j], sep=''), full.names = T)
  for(i in 1:length(my.files)){
    orig <- image_read(my.files[i])
    nam <- my.files[i] %>% strsplit(., split = '.tif') %>% lapply('[[',1) %>% unlist
    image_write(orig, path = paste(nam, '.png',sep=''), format = 'png')
  }
}


