#work using imager package
#working through:http://dahtah.github.io/imager/canny.html
# and here for pix sets https://cran.r-project.org/web/packages/imager/vignettes/pixsets.html
library(imager)
library(magrittr)
library(magick)
library(purrr)
library(dplyr)
library(EBImage)

#get list of images to work on
my.images <- list.files('./common/rob_2018_fixed', pattern='png', full.names = T)

# #convert to png
# for(i in 1:length(my.images)){
#   img <- image_read(my.images[i])
#   img <- image_convert(img, "png")
#   nam <- my.images[i] %>% strsplit(., '[.]') %>% lapply('[[',2) %>% unlist
#   nam <- paste('.', nam, '.png', sep='')
#   image_write(img, path = nam, format = 'png')
# }


#read in as png using imager package
img <- load.image(my.images[7])
my.chan <- channels(img, drop=F)
img <- my.chan[1:3] %>% imappend("c")

#convert img to grayscale
img <- grayscale(img) #%>% isoblur(3)
plot(img)

#illumination correction using linear model.
d <- as.data.frame(img)
##Subsamble, fit a linear model
m <- sample_n(d,1e4) %>% lm(value ~ x*y,data=.)
##Correct by removing the trend
im.c <- img - predict(m,d)

#edge detection
#imgradient(img, "xy") %>% enorm %>% plot()
#calculate image gradient along x and y
gr <- imgradient(im.c, "xy")
#plot(gr, layout = 'row')
#the gradient has two components (the image derivatives along x and y) and is stored as an image list with two components
#we can compute the gradient magnitude via
mag <- with(gr, sqrt(x^2+y^2))
#plot(mag)
#threshold(mag, approx = F) %>% plot

#this gives us just the grid outline - works very well
thr <- threshold(mag, approx = F)
#plot(thr)
#we can then use morphological dilation to grow the area so that it encapsulates all associated pixels
thr <- grow(thr, 15) #remember we could do this individually per channel if that is better?
#plot(thr)

#now i want to use threshold image of mag that has been grown (showing highlighted grid) to mask pixels to NA values in original image, then use inpaint to fill these in with surrounding colour, effectively removing the grid!

#plot original with highlighted grid
#plot(im.c)
#highlight(thr)
#or use colourise to highlight on plot
#colorise(img, thr, 'red', alpha=.5) %>% plot

#to make thr the reciprocal can all 1 - thr (i.e. changes to TRUES of the grid to FALSE and the reverse, so that we can use it as a mask for our original image.)
#make intermediate reciprolca mask, the set TRUE value to be something negative
recip <- 1 - thr
grd.rm <- img * recip
#plot(grd.rm)
grd.rm[grd.rm == 0] <- NA
#plot(grd.rm)

test <- inpaint(grd.rm, 100) # perhaps could be dynamic
plot(test)


#Now i need to start again with this image and see if we can subset out the algae!!
#calculate image gradient along x and y
grad <- imgradient(test, "xy")
#plot(grad, layout = 'row')
new.mag <- with(grad, sqrt(x^2+y^2))
#plot(new.mag)
new.thr <- threshold(new.mag, thr = '99.5%', approx = F) #might need to be dynamic
new.thr <- clean(new.thr, 2)
new.thr <- grow(new.thr, 15) #should be dynamic - definitely
new.thr <- shrink(new.thr, 5)
plot(im.c)
highlight(new.thr)

#can i convert to raster for nicer plotting in shiny, perhaps as a leafmap?
gg <- img %>% as.array %>% Image # works a treat, can just punt to array!!
hh <- new.thr %>% as.array %>% Image
segs <- hh %>% bwlabel
paintObjects(segs, gg,opac = 0.8, col='black') %>%  display()

#write out test and read in the EBImage to process like before...
save.image(test, file = './common/rob_2018_fixed/test.png', quality = 1)
save.image(as.cimg(new.thr), file = './common/rob_2018_fixed/test_blobs.png', quality = 1)

#when i have used as similar approch to select the algae in the image, I could use split_connected on the px to separate these into individual images...can also use grow to expand a pixset or shrink to contract a pix set! and bbox(pixset) to compute the bounding box for a whole set or if split into individual particles, for an individual!


#switching to EBImage
img <- readImage(my.images[7])
blobs <- readImage('./common/rob_2018_fixed/test_blobs.png')

segs <- blobs %>% bwlabel
#table(segs)
display(colorLabels(segs))
paintObjects(segs, img[,,1:3],opac = 0.8, col='black') %>%  display()

#can we have a look at some features...
fea_int <- computeFeatures.basic(x = segs, img) %>% as.data.frame #gives the pixel intensities (mean, sd, mad, quantile)
fea_shape <-computeFeatures.shape(x = segs, img) %>% as.data.frame #gives the object features (area, perimeter, mean radius, sd, max radius, min radius; all in pixel values)
fea_moment <- computeFeatures.moment(x = segs, img) %>% as.data.frame #gives center mass ox x, center mass of y, elliptical fit of major acxis, elliptical eccentricity, object angle
fea_haralick <- computeFeatures.haralick(x = segs, img) %>% as.data.frame #computes features that quantify pixel texture according to Haralicks original paper...R. M. Haralick, K Shanmugam and Its'Hak Deinstein (1979). Textural Features for Image Classification. IEEE Transactions on Systems, Man and Cybernetics.

#lets have a look at each of our objects...
tt <- stackObjects(segs, img) #can then click through these in the viewer! Might want to add a  buffer to aid in identification...NB there is a shiny app commad for the EBImage viewer!!
display(tt)

#which of segs are edge effects..
#2, 11, 27 (smaller),  29, 33 (though has a cell attached), 34, 35 (smaller), 38 + 40 (v small), 46 (small), 47, 50,

#fea int features
plot(1:nrow(fea_int), fea_int$b.mean) #some have mean pixel intensity of 1 = only white pixels = not good
#to remove = which(fea_int$b.mean == 1) - this gets rid of 32, 35, 38, 40, 46, 52 - a lot of the smaller non cell ones!
plot(1:nrow(fea_int), fea_int$b.sd) #same objects have zero SD = only white pixels!
plot(1:nrow(fea_int), fea_int$b.mad) # and zero mad
plot(1:nrow(fea_int), fea_int$b.q095) #etx same for all intensity metrics - so good way to remove these initially

#fea shape features - perhaps a dynamic slider for some of these (radius max perhaps?)
plot(1:nrow(fea_shape), fea_shape$s.area) #most high ones not of interest, but no real clear distinction
plot(1:nrow(fea_shape), fea_shape$s.perimeter) # > 200 pixels could work, but we lose debris bundle
plot(1:nrow(fea_shape), fea_shape$s.radius.mean) # again > 30 here
plot(1:nrow(fea_shape), fea_shape$s.radius.sd) # > 20 = 2, 11, 33, 34, 47 - all to remove
plot(1:nrow(fea_shape), fea_shape$s.radius.min) # no pattern
plot(1:nrow(fea_shape), fea_shape$s.radius.max) # gets several > 50 all to remove

#fea moment features
plot(1:nrow(fea_moment), fea_moment$m.majoraxis) #>150 gives good - can we calc SD and take those over eg 1 x sd? (gives 131)?
which(fea_moment$m.majoraxis > sd(fea_moment$m.majoraxis)) #all caught are to remove - could be applied across more boradly
plot(1:nrow(fea_moment), fea_moment$m.eccentricity)#no pattern at all
plot(1:nrow(fea_moment), fea_moment$m.theta)# no pattern at all

#fea_haralick features
plot(1:nrow(fea_haralick), fea_haralick$h.asm.s1) #anything that == 1 is bollocks
plot(1:nrow(fea_haralick), fea_haralick$h.con.s1) #same ones == 0 on this
plot(1:nrow(fea_haralick), fea_haralick$h.cor.s1) #same ones == 0
plot(1:nrow(fea_haralick), fea_haralick$h.var.s1) #no discernable pattern
plot(1:nrow(fea_haralick), fea_haralick$h.idm.s1) #not helpful
plot(1:nrow(fea_haralick), fea_haralick$h.sav.s1) #not helpful
plot(1:nrow(fea_haralick), fea_haralick$h.sen.s1) #not helpful
plot(1:nrow(fea_haralick), fea_haralick$h.ent.s1) #not helpful
plot(1:nrow(fea_haralick), fea_haralick$h.asm.s2) #rest are about the same, not so great.

###features to remove
#based on fea_int
fi.rm <- which(fea_int$b.mean == 1 & fea_int$b.sd == 0) #this is reliable for edge effects
my.rm <- c(fi.rm, fs.rm) %>% unique %>% sort # not bad only misses a couple of upright, non white ones

segs <- rmObjects(segs, my.rm)
paintObjects(segs, img[,,1:3],opac = 0.8, col='black') %>%  display()


#2, 11, 27 (smaller),  29, 33 (though has a cell attached), 34, 35 (smaller), 38 + 40 (v small), 46 (small), 47, 50,

###########Using features to separate out non-cell like objects before making training dataset...

#haralick names
#asm - angular second moment
#con - contrast
#cor - correlation
#var - variance
#idm - inverse difference moment
#sav - sum average
#sva - sum variance
#sen - sum entropy
#ent - entropy
#dva - difference variance
#den - difference entropy
#f12 -
#f13

#can I chuck these in a PCA (or machine learning algorithm / classification algorithm) to classify the cells
dat <- cbind(fea_shape, fea_haralick[1:13])
pp1 <- prcomp(formula(paste('~',paste0(names(dat), collapse='+'))), data = dat, scale = T)
summary(pp1) #good amount of variance explained by PC3 (96%)
biplot(pp1)
scores <- pp1$x %>% as.data.frame
#obvious cluster at left of PC1 (< -0.2)
#drop them from the segs and re-plot to see who disappears - it was the straight ones on righ tof image - not got rid of right ones - add in more to PCA?
segs <- rmObjects(segs, which(scores$PC1 < - 0.2))
paintObjects(segs, img[,,1:3],opac = 0.8, col='black') %>%  display()
#and some < -5 PC2
segs <- rmObjects(segs, which(scores$PC2 < - 5))
paintObjects(segs, img[,,1:3],opac = 0.8, col='black') %>%  display()


x <- stack()
for(i in 1:(dim(f)[4])){
  print(i)
  #plot(h)
  h <- getFrame(f, i = i, type='render')[,,1]
  writeImage(h, files = 'tmp.tiff', type='tiff')
  k <- raster('tmp.tiff')
  plot(k)
  names(k) <- i
  x <- stack(x, k)
}
