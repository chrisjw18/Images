# Chris script on G's computer to segment algal images

#install bioconductor and use bioLite() to install EBImage package
source("https://bioconductor.org/biocLite.R")
biocLite('CRImage') #use this line of code to install specific parts of bioconductor
install.packages('magick')

library(EBImage)
library(CRImage)
library(magrittr)
library(magick)
library(spatialfil) #has sobel edge detection...masks EBImage display()
display <- EBImage::display

input.dir <- './common/unfixed/'
output.dir <- './common/unfixed_processed/'

my.images <- list.files(input.dir)

#pdf(file = paste(output.dir, 'test1.pdf', sep=''))

for(i in 1:length(my.images)){

#load as bmp with magick
img <- paste(input.dir, my.images[i], sep='') %>% image_read

#convert to tiff and write out as same name
img %>% image_write('./common/tmp.tiff', format='tiff')

#read back in with EBImage
img <- './common/tmp.tiff' %>% readImage()
#display(img)

#resize image
img_re <- img %>% resize(500,500)

#get the array
dat <- imageData(img_re)

#test the different filters
sob_original <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'sobel')) %>% Image
sob_original  %>% display #this seems to give really good results - need to use this as an 'edge layer'

#threshold sob > mean
sob <- sob_original > 0.5 #gives T/F, turn binary
sob[sob==FALSE] <- 0
sob[sob==TRUE] <- 1
display(sob)
sob <- dilate(sob, kern= makeBrush(size=5, shape='disc')) #can be too much for close cells...
display(sob)

#fill the complete regions
sob <- fillHull(sob)
display(sob,all=T)

#remove small regions?
segs <- sob %>% bwlabel
#table(segs)
#display(colorLabels(segs))
paintObjects(segs[,,1], img_re,opac = 0.8, col='black') %>%  display(method = 'raster')
#num.to.remov <- #need to remove objects that are smaller than threshold using table info but dont know how to subset the table///

table(segs)
#can we have a look at some features...
fea_int <- computeFeatures.basic(x = segs[,,1], img_re) %>% as.data.frame #gives the pixel intensities (mean, sd, mad, quantile)
fea_shape <-computeFeatures.shape(x = segs[,,1], img_re) %>% as.data.frame #gives the object features (area, perimeter, mean radius, sd, max radius, min radius; all in pixel values)
fea_moment <- computeFeatures.moment(x = segs[,,1], img_re) %>% as.data.frame #gives center mass ox x, center mass of y, elliptical fit of major acxis, elliptical eccentricity, object angle
fea_haralick <- computeFeatures.haralick(x = segs[,,1], img_re) %>% as.data.frame #computes features that quantify pixel texture according to Haralicks original paper...R. M. Haralick, K Shanmugam and Its'Hak Deinstein (1979). Textural Features for Image Classification. IEEE Transactions on Systems, Man and Cybernetics.

#can i use the center of mass from fea_moment as seeds in a propagate routine, with segs [binary] image as the mask. WHich image to use as the image - perhaps the ORIGINAL sobel edge image [sob_original]?
my.seeds <- array(0, dim=c(500,500))
for(j in 1:nrow(fea_moment)){
  my.seeds[round(fea_moment$m.cx[j],0), round(fea_moment$m.cy[j],0)] == j
}
test <- propagate(x = sob_original[,,1], seeds = my.seeds, mask = segs[,,1]) #this didnt work yet.....keep trying.
display(test)

}#end of i loop
dev.off()
