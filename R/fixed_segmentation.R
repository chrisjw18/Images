#image segmentation pipeline for those images Robb has been producing from fixed samples

#example images are stored in common>rob_2018_fixed

library(EBImage)
library(CRImage)
library(magrittr)
library(magick)
library(spatialfil) #has sobel edge detection...masks EBImage display()
display <- EBImage::display
library(raster)

#get list of images to work on
my.images <- list.files('./common/rob_2018_fixed/', full.names = T)


#read in tiff image with  EBImage
#img <- my.images[1] %>% readImage()
img <- readImage('./common/rob_2018_fixed/test_blobs.png')
display(img)

#resize image (why?)
img_re <- test %>% resize(1000,1000)

#gray scale image?
img_gry <- EBImage::channel(img_re[,,1:3], "gray")
display(img_gry < 0.2) #this is a really good way to threshold!! The cells are dark compared to both the grid and the background!!

#get the array
dat <- imageData(img_re)

#test the different filters
sob_original <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'sobel')) %>% Image
sob_original  %>% display #this seems to give really good results - need to use this as an 'edge layer'

#maybe need to separate out the grids into a different image to deduct later?
my.grd <- sob_original[,,1] #> 0.9
my.grd[my.grd < 0.9] <- 0
display(my.grd)
my.grd <- closing(my.grd, kern= makeBrush(size=17, shape='disc'))
display(my.grd)
#the below switches the binary image so we can multiply out (as zeros) the grid area from other images
my.grd[my.grd == 0 ] <- 100
my.grd[my.grd < 100] <- 0
my.grd[my.grd == 100] <- 1

#threshold sob > mean
sob <- sob_original > 0.2 #gives T/F, turn binary
sob[sob==FALSE] <- 0
sob[sob==TRUE] <- 1
display(sob)
sob <- closing(sob, kern= makeBrush(size=7, shape='disc')) #can be too much for close cells...
display(sob)

#remove grid here? prior to applying bwlabel?
sob[,,1] <- sob[,,1] * my.grd
display(sob)

#fill the complete regions - doesnt work due to the grids in the fixed images...skip for now
sob <- fillHull(sob)
display(sob,all=T)

#remove small regions?
segs <- sob %>% bwlabel
table(segs)
display(colorLabels(segs))
paintObjects(segs[,,1], img_re,opac = 0.8, col='black') %>%  display()
#num.to.remov <- #need to remove objects that are smaller than threshold using table info but dont know how to subset the table///
x <- table(segs)
my.shapes <- 0:(length(x)-1)
my.sizes <-c()
for(i in 1:length(x)){
  my.sizes[i] <- x[i][[1]]
}
to.correct <- my.sizes > 2000
#loop through shape numbers and change to 0 if size associated with shape number is above our threshold above (currently 500)
for(i in 1:length(my.shapes)){
  if(to.correct[i] == TRUE){
    segs[segs == my.shapes[i]] <- 0
  }
}
segs <- segs %>% bwlabel
display(colorLabels(segs))
paintObjects(segs[,,1], img_re[,,1],opac = 0.9, col='purple') %>%  display()



table(segs)
#can we have a look at some features...
fea_int <- computeFeatures.basic(x = segs[,,1], img_re) %>% as.data.frame #gives the pixel intensities (mean, sd, mad, quantile)
fea_shape <-computeFeatures.shape(x = segs[,,1], img_re) %>% as.data.frame #gives the object features (area, perimeter, mean radius, sd, max radius, min radius; all in pixel values)
fea_moment <- computeFeatures.moment(x = segs[,,1], img_re) %>% as.data.frame #gives center mass ox x, center mass of y, elliptical fit of major acxis, elliptical eccentricity, object angle
fea_haralick <- computeFeatures.haralick(x = segs[,,1], img_re) %>% as.data.frame #computes features that quantify pixel texture according to Haralicks original paper...R. M. Haralick, K Shanmugam and Its'Hak Deinstein (1979). Textural Features for Image Classification. IEEE Transactions on Systems, Man and Cybernetics.

#can i use the center of mass from fea_moment as seeds in a propagate routine, with segs [binary] image as the mask. WHich image to use as the image - perhaps the ORIGINAL sobel edge image [sob_original]?
my.seeds <- array(0, dim=c(1000,1000))
for(j in 1:nrow(fea_moment)){
  my.seeds[round(fea_moment$m.cx[j],0), round(fea_moment$m.cy[j],0)] == j
}
test <- propagate(x = sob_original[,,1], seeds = my.seeds, mask = segs[,,1]) #this didnt work yet.....keep trying.
display(test)
