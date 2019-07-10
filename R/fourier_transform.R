#attempt at discrete fourier transform of an image to remove the haemocytometer grid, as per...
#https://mathematica.stackexchange.com/questions/88494/how-do-i-remove-grid-from-this-photo

library(EBImage)
library(CRImage)
library(magrittr)
library(magick)
library(spatialfil) #has sobel edge detection...masks EBImage display()
display <- EBImage::display
library(raster)
library(stats)
library(imager)

#get list of images to work on
my.images <- list.files('./common/rob_2018_fixed/', full.names = T)


#read in tiff image with  EBImage
img <- my.images[1] %>% readImage()
display(img)

#gray scale image?
img_gry <- EBImage::channel(img[,,1:3], "gray")
display(img_gry < 0.2)

fft <- fft(img, inverse = T)
fft <- fft(img_gry, inverse = T)
display(fft)
