# Chris script on G's computer to segment algal images

#install bioconductor and use bioLite() to install EBImage package
source("https://bioconductor.org/biocLite.R")
biocLite('CRImage') #use this line of code to install specific parts of bioconductor
install.packages('pixmap')

library(EBImage)
library(CRImage)
library(magrittr)
#library(readbitmap)
#library(pixmap)

#load example image - EBImage does not seem to like BMP files!
img <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A001 - 20180723_144843.bmp' %>% read.bitmap()

img %>% pixmapRGB %>% plot #cant seem to get this img array to an EBImage format..

#This is the same image, but manually converted to TIFF using preview.
img2 <- './common/microscope_images_unfixed/220718_5am_1/A001 - 20180723_144843.bmp' %>% readImage

#to make neagtive, max(img2)-img2
#to increase brightness img2 + 0.3 (Addition)
#to increase contrast img2 * 2 (multiplication)
#to apply gamma correction img2 ^ 0.5 (exponentiation)
#can crop and threshold images with standard matrix operations
display(img2 > 0.36) #here for example, this treshold cuts out everything but the ice algae!
display(img2)

#filtering - the makeBrush function generates brushes of various sizes and shapes that can be used as structuring elements. There are many options in the help here some for binary images some for grayscale. Morphological functions position the center of the structuring element over each pixel in the inputimage. If pixels are at the edge of an image, part of the neighbourhood defined by the structuring element may extend past the border image, in which case a value is assigned to these undefined pixels as if the image was padded with additional rows and columns.

#makeBrush generates a 2D matrix continaing the desired brush; size = size of brush in pixels; shape = shape of brush, can be box, dice, diamond, guassian or line; step = logical whether brush is binary; sigma = optional numeric containing the SD of the guassian shape (defaul = 0.3)


#then use filter() with the brush - filter2 is a 2D convolution filter (filters an image using the fast 2D FTT convolution product.
w <- makeBrush(size=31, shape='gaussian', sigma=5)
img_flo <- filter2(img2, w)
img_flo %>% display() #this has made a blurry image

#linear filtering is used to perform low-pass filtering (to blur images, remove oise),
x = readImage(system.file("images", "sample-color.png", package="EBImage"))
display(x, title='Sample')

## Low-pass disc-shaped filter
f = makeBrush(21, shape='disc', step=FALSE)
display(f, title='Disc filter')
f = f/sum(f)
y = filter2(x, f)
display(y, title='Filtered image')

## Low-pass filter with linear padded boundary
y = filter2(x, f, boundary=c(0,.5,1))
display(y, title='Filtered image with linear padded boundary')

#and high-pass filtering (to detect edges, sharpen images).
## High-pass Laplacian filter
la = matrix(1, nc=3, nr=3)
la[2,2] = -7.8
y = filter2(x, la)
display(y, title='Filtered image')

## High-pass Laplacian filter with replicated boundary
y = filter2(x, la, boundary='replicate')
display(y, title='Filtered image with replicated boundary')


#test the high-pass Laplacian filter on the ice algal image
img_hp <- img2 %>% filter2(la)
img_hp %>% display
display(img_hp*0.5)

#make image grayscale, then binary classify using otsu thresholding (clustering-based image thresholding, looks for bimodal distribution of image inteisities)
img_g <- img2 %>% channel(mode = 'gray')
img_g %>% display
thres <- img_g %>% otsu
img_thres <- img_g > thres
img_thres %>% display
display(combine(img2, img_thres))
#adaptive thresholding - threshold allowed to be different in different regions of the image, to account for e.g uneven illumination or stray signal from nearby bright objects. It compares each pxels intensity to the background determined from a local neighbourhood. This can be achieved by comparing the image to its smoothed version, where the filtering window is bigger than the typical size of objects we want to capture. Thres function uses rectangular
img2_thres <- img2 %>% thresh(w=20, h=20) #this does not do well, it does not get the whoel cell, more like the edge of the cell!
img2_thres %>% display
hist(img2_thres)
display(img2_thres==1)

#perhaps otsu thresholding after low-pass filtering - yes, this works the best so far.
f = makeBrush(21, shape='disc', step=FALSE)
f = f/sum(f)
y = filter2(img2, f)
display(y, title='Filtered image')
img_g <- y %>% channel(mode = 'gray')
img_g %>% display
thres <- img_g %>% otsu
img_thres <- img_g > thres
img_thres %>% display

img_thres %>% distmap %>% watershed %>% display()
img2 %>% display

#watershed approach; watershed transformation treats a grayscale image as a topographic relief, or heightmap. Objects that standout of the background are identified and separated by flooding an inverted source image.
img_w <- img_g %>% distmap
img_w %>% display
img_g %>% display

#applying my own mask visually...
nmask <- img2 < 0.36 #I can certainly apply this threshold here to create a T/F mask for training a neural network - right?!
nmask %>% display
nmask <- opening(nmask, makeBrush(5, shape='disc')) #this was good to get rid of particle
nmask %>% display
nmask <- fillHull(nmask)
nmask %>% display

#none-touching connected objects can be segmented using bwlabel (see watershed/propagate for segmenting touching objects)
nlabs <- bwlabel(nmask)
table(nlabs) #here is just 0 and 1, but would be e.g. 0 = bg, 1, 2, 3, 4, if more cells...
#to display we tend to normalise to the 0,1 range expected by the display function = different shades of grayscale for each object....
nlabs %>% normalize %>% display


#try on anothe image
img3 <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A002 - 20180723_144859.tiff' %>% readImage()
img3 %>% display
img3_g <- img3 %>% channel(mode='gray')
img3_g %>% hist
img3_g <- img3_g < 0.5
display(img3_g)
img3_g <- opening(img3_g, makeBrush(5, shape='gaussian')) #this was good to get rid of particle
img3_g %>% display
img3_g <- fillHull(img3_g)
img3_g %>% display
img3_g_labs <- img3_g %>% bwlabel()
display(colorLabels(img3_g_labs))

#and another
img4 <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A004 - 20180723_144927.tiff' %>% readImage()
img4 <- img4 * 1.8 #increase contrast
img4 %>% display

img4_g <- img4 %>% channel(mode='gray')
img4_g %>%  display()
img4_g <- img4_g < 0.9
img4_g %>% display()
img4_g <- opening(img4_g, makeBrush(10, shape='disc'))
img4_g <- fillHull(img4_g)
img4_g %>% display()
img4_g <- img4_g %>% bwlabel

#can i use img4_g as the seeds for
seed <- img4_g
display(seed)
cells <- img4 %>% channel(mode='gray')
ctmask <- opening(cells<0.9, makeBrush(5, shape='disc'))
ctmask %>% display()
cmask <- propagate(cells, seeds = seed, mask = ctmask, lambda = 1e-12)
segmented <- paintObjects(cmask, cells, col='red')
display(segmented)
segmented <- paintObjects(seed, segmented, col='red')
display(segmented)


#quick test with a fixed image!
img <- '../Image_segmentation/fixed_images/image.tif' %>% readImage()
#display(img*1.5)
display(img)
imgg <- img %>% channel(mode='gray')
display(imgg)
imgg <- imgg < 0.45
imgg <- opening(imgg, makeBrush(10, shape='disc'))
display(imgg)
imgg <- fillHull(imgg)
imgg %>% display
imgg_labs <- imgg %>% bwlabel()

seed <- img %>% channel(mode='gray')
seed <- seed < 0.45
seed <- opening(seed, makeBrush(10, shape='disc'))
display(seed) #this needs to be more selective

cells <- img %>% channel(mode='gray')
ctmask <- opening(cells<0.5, makeBrush(10, shape='disc'))
ctmask %>% display() #make sure mask has sufficient breadth to allow expansion to outter of cells
cmask <- propagate(cells, seeds = seed, mask = ctmask, lambda = 1e-4)
segmented <- paintObjects(cmask, cells, col='red')
display(segmented)
segmented <- paintObjects(seed, segmented, col='red')
display(segmented)


library('spatialfil') #has sobel edge detection...
display <- EBImage::display
# • gaussian for Gaussian kernel
# • LoG for Laplacian of Gaussian kernel
# • sharpen for 3x3 convolution matrix for sharpening edges
# • laplacian for a 3x3 convolution matrix that enhances the edges
# • emboss for a 3x3 kernel that draws edges as embossed image
# • sobel gives one of the two 3x3 matrices needed to apply the Sobel filter - for this one, #Sobel convolution kernel returns the possibility to detect edges in a more sofisticated way, #the convKernel function returns only one of the two matrices needed to apply the filter. The #second one is calculated by transposing the returned matrix in the other needed one. i.e. it #makes the horizontal matrix and we trnaspose this for the vertical??

img_re <- resize(img2, 300, 300) #requires a square image!
plot(img_re)
dat <- imageData(img_re) #get the array from the image

#test the different filters
gau <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'gaussian'))
gau %>% Image %>% display
log <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'LoG'))
log %>% Image %>% display
sharp <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'sharpen'))
sharp %>% Image %>% display
lap <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'laplacian'))
lap %>% Image %>% display
emb <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'emboss'))
emb %>% Image %>% display
sob <- applyFilter(dat, kernel = convKernel(sigma = 1.4, k = 'sobel'))
sob %>% Image %>% display #this seems to give really good results

#test sobel with other images
i3 <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A004 - 20180723_144927.tiff' %>% readImage()
i3 %>% display
i3 %>% resize(500,500) %>% imageData %>% applyFilter(kernel = convKernel(k='sobel')) %>% Image %>% display

i4 <- '../Image_segmentation/fixed_images/image.tif' %>% readImage()
i4 %>% display()
i4 %<>% resize(500,500) %>% imageData %>% applyFilter(kernel = convKernel(k='sobel')) %>% Image
i4 %>% display
display(i4>1)

i5 <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A005 - 20180723_144936.tiff' %>% readImage() %>% channel(mode='gray')
i5 %>% display()
i5 %<>% resize(500,500) %>% imageData %>% applyFilter(kernel = convKernel(k='sobel')) %>% Image
i5 %>% display


i5 <- i5 > 0.5 #will need a threshold to implement the sobel filter.
#i5 <- fillHull(i5)
display(i5)
#erode and dilate
ii<- erode(i5, kern = makeBrush(1, 'disc'))
ii %>% display
iii <- erode(ii, kern = makeBrush(1, 'disc'))
iii %>% display
ii <- dilate(ii, kern = makeBrush(3, 'disc'))
ii %>% display()
ii <- fillHull(ii)
ii %>% display
#found a new package called CRImage in bioconducter that is to classify cells in biological images!!
source("https://bioconductor.org/biocLite.R")
biocLite('CRImage')
library(CRImage)
img2 <- '../Image_segmentation/microscope_images_unfixed/220718_5am_1/A001 - 20180723_144843.tiff' %>% readImage
plot(img2)

#see how the segmentImage function goes - not bad but not great!
test <- segmentImage(image=img2, threshold = 'otsu', numWindows = 1, maxShape = 40, greyscaleImage = 2)

test$segmentedImage %>% Image %>% plot
test$features %>% nrow

fea <- test$features %>% as.data.frame

img3 %>% display


#also need to look into the simpleITK package....
install.packages('devtools')
devtools::install_github("SimpleITK/SimpleITKRInstaller")
