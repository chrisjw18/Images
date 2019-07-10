#shiny app version
library(shiny)
library(shinydashboard)
library(shinyAce)
library(shinyFiles)
library(imager)
library(magrittr)
library(magick)
library(purrr)
library(dplyr)
library(EBImage) #*** This is a doctured version of the package from install_github("jkh1/EBImage") that allows locating pixel position when clicking. Devtools will need to be installed inorder to install from github directly.
library(abind)
library(snow)
library(mapview)
library(leaflet)
library(leaflet.extras)
library(rgdal)
library(sp)
library(raster)
library(DT)
library(stringr)
library(data.table)
library(shinyBS)
library(RColorBrewer)


#this increases uplaod file size to 30MB
options(shiny.maxRequestSize = 30*1024^2)
leafletProj <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +towgs84=0,0,0,0,0,0,0 +units=m +nadgrids=@null +wktext +no_defs"#(re) define projections that will be used in the polgon fitting within the app
#define polygon projection
poly.proj <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

ui <- dashboardPage(

  dashboardHeader(title = "Williamson Cells"),

  dashboardSidebar(
    #dynamic sidebarMenu driven by folder heirachry
    sidebarMenu(id = 'sidebarmenu',
                menuItem("Segment Cells", tabName = "segmentCells",  icon = icon("group", lib="font-awesome")),
                menuItem("Cells ID", tabName = "cellId", icon = icon("check-circle", lib = "font-awesome")),
                menuItem("Summary Statistics", tabName = "summary", icon = icon("check-circle", lib = "font-awesome")),
                conditionalPanel("input.sidebarmenu === 'segmentCells'",
                                 fileInput('fileIn', '1. Upload PNG Image file',
                                           multiple = F, accept = ".png"),
                                 sliderInput('thresBar', '2. Magnitude Threshold',
                                             97, 100, value = 99.5, step = 0.1),
                                 sliderInput('growBar', '3. Grow Threshold',
                                             0, 20, value = 10, step = 1),
                                 sliderInput('shrinkBar', '4. Shrink Threshold',
                                             0, 20, value = 5, step = 1),
                                 radioButtons('scale_input', label = '6. Record scale', choices = c('Y','N'), selected = 'Y', inline = T),
                                 actionButton('reset_scale', 'Reset scale', width = '200px'),
                                 textOutput('scale_text', inline = T)
                ),
                conditionalPanel("input.sidebarmenu === 'cellId'",
                                 selectInput('id', 'Species',
                                             choices = c('Non algal', 'm1', 'm2', 'm3', 'm4',
                                                         'a1', 'a2', 'a3', 'a4', 'a5', 'a6',' a7', 'a8',
                                                         'a9', 'a10', 'a11', 'a12', 'c1','c2','c3','c4',
                                                         'snow', 'bluegreen','other'), selected = 'm1'),
                                 radioButtons('rec', 'Record Input', choices = c('Y','N'), selected = 'Y', width = '200px', inline = T),
                                 actionButton('comp', 'Compute Summary', width = '200px')
                ),
                conditionalPanel("input.sidebarmenu === 'summary'",
                                 shinyDirButton('wd', 'Select WD', title = 'Select WD'),
                                 actionButton('export', 'Export Data', width = '200px')
                )
    )
  ),#end of dashboardSidebar

  dashboardBody(
    tabItems(
      tabItem(
        tabName = 'segmentCells',
        displayOutput('img1', width = '100%', height = '600px')
      ),
      tabItem(tabName = 'cellId',
              displayOutput('img2', width = '100%', height = '600px')
      ),
      tabItem(tabName = 'summary',
              fluidRow(
                column(12,
                       renderText('Summary statistics'),
                       verbatimTextOutput('summary_text')
                ),
                column(12,
                       plotOutput('img3', width = '100%', height = '600px')
                ),
                column(12,
                       plotOutput('sp_abundance', 'Relative species abundance', width = '100%', height = "200px")
                ),
                column(12,
                       dataTableOutput('summary_table', width = '100%'))
              ))
    )#tabItems
  )#dashboardBody

)#end of UI


server <- function(input, output, session) {

  app.data <<- reactiveValues(data = NULL)

  #establish catcher dataframes needed later.
  observe({
    app.data$scale_input <- data.frame(x = numeric())
    app.data$cell_data <- data.frame(id = character(), x = numeric(), y = numeric(), blob = numeric())
  })

  #read in file; do initial cleaning of image; do initial cell segmentation
  observeEvent(input$fileIn, {
    #get file details - also makes this reactive on input$fileIn
    my.file <- input$fileIn
    app.data$file_data <- my.file
    my.nam <- my.file$name %>% as.character %>% strsplit(., split='[.]') %>% lapply('[[',1) %>% unlist
    app.data$file_name <- my.nam
    app.data$data_path <- my.file$datapath

    #load the data as imageR image
    my.img <- load.image(my.file$datapath)
    if(dim(my.img)[4] > 3){
      my.chan <- channels(my.img, drop=F)
      my.img <- my.chan[1:3] %>% imappend("c")
    }
    app.data$imager_img <- my.img

    #load copy in EBImage format
    orig_img <- readImage(my.file$datapath)

    #make gray scale image
    gry_img <- grayscale(my.img)

    #illumination correction using linear model.
    d <- as.data.frame(gry_img)
    m <- sample_n(d,1e4) %>% lm(value ~ x*y,data=.)
    im.c <- gry_img - predict(m,d)

    #edge detection
    #calculate image gradient along x and y
    gr <- imgradient(im.c, "xy")
    #compute the gradient magnitude via
    mag <- with(gr, sqrt(x^2+y^2))
    #threshold magnutide automatically
    thr <- threshold(mag, approx = F)
    #grow to encapsulate all pixels - dynamic?
    thr <- grow(thr, 15) #remember we could do this individually per channel if that is better?
    #make intermediate reciprocal mask
    recip <- 1 - thr
    #make grid to remove
    grd.rm <- gry_img * recip
    #set pixels to NA for inpaint function
    grd.rm[grd.rm == 0] <- NA
    #inpaint haemo grid
    grd_rm_img <- inpaint(grd.rm, 100) # perhaps could be dynamic

    #push to app.data
    app.data$original_img <- orig_img
    app.data$gray_img <- gry_img
    app.data$grid_removed_img <- grd_rm_img

  })


  #observeEvent(c(input$thresBar, input$growBar, input$shrinkBar),{
  observeEvent(input$thresBar, {
    if(!is.null(app.data$grid_removed_img)){
      #print('test1')
      app.data$thresBar <- input$thresBar
      grad <- imgradient(app.data$grid_removed_img, "xy")
      new.mag <- with(grad, sqrt(x^2+y^2))
      my.thres <-  input$thresBar %>% as.numeric
      #print(my.thres)
      my.thres <- paste(my.thres, '%', sep='')
      new.thr <- threshold(new.mag, thr = my.thres,approx = F) #dynamic
      new.thr <- clean(new.thr, 2)
      grow.thres <- input$growBar
      new.thr <- grow(new.thr, as.numeric(input$growBar)) #dynamic
      new.thr <- shrink(new.thr, as.numeric(input$shrinkBar)) #dynamic
      app.data$final_total_blobs_imager <- new.thr #Imager - single pixel set of total blobs (previously my_blobs)


      #extract blobs as EBImage product
      segs <- new.thr %>% as.cimg %>% as.array
      segs <- adrop(segs, drop = c(3,4)) %>% Image #EBImage
      segs <- segs %>% bwlabel #EBImage
      all_blobs <- stackObjects(segs, app.data$original_img) #EBImage

      #calculate summary statistics per blob
      #pixel intensities (mean, sd, mad, quantiles)
      app.data$final_fea_int <- computeFeatures.basic(x = segs, app.data$original_img) %>% as.data.frame
      #shape features (area, perimeter, mean radius, sd, max radius, min radius) - all in pixel values
      app.data$final_fea_shape <- computeFeatures.shape(x = segs, app.data$original_img) %>% as.data.frame
      #location features (center mass x/y, elliptical fit of major axis, elliptical eccentricity, object angle)
      app.data$final_fea_moment <- computeFeatures.moment(x = segs, app.data$original_img) %>% as.data.frame
      #haralick texture parameters
      app.data$final_fea_haralick <- computeFeatures.haralick(x = segs, app.data$original_img) %>% as.data.frame

      #identify which

      #push all to app.data
      app.data$final_total_blobs_ebimage_single_img <- segs #EBImage; single image with all blobs (previously total_blobs)
      #app.data$all_blobs <- all_blobs #EBImage; image stack of separated blobs - dont need this until later...

      #new columns for catching cell ids later
      app.data$all_data <- app.data$final_fea_shape
      app.data$all_data$blob <- NA #this will be sanity check right matches are made
      app.data$all_data$id <- NA

    }

  })


  #record pixel x location from click position when input$scale radio button is set to Y
  observeEvent(input$pixelPosition,{
    if(input$scale_input == 'Y'){
      new.x <- data.frame(x = input$pixelPosition[1])
      if(nrow(app.data$scale_input) < 2){
        app.data$scale_input <- rbind(app.data$scale_input, new.x)
      }
    }
  })
  #calculate pixelPosition when 2 values present in dataframe
  observe({
    input$pixelPosition
    if(nrow(app.data$scale_input) == 2){
      pixel.length <- diff(c(app.data$scale_input$x[1],app.data$scale_input$x[2]))
      app.data$pixel_length_um <- 250 / pixel.length #*********THIS ASSUMES 1 SMALL SQUARE is used for the scale-setting, i.e. 250 um in length!!******#######
    }
  })
  #reset scale dataframe in the case of error
  observeEvent(input$reset_scale,{
    app.data$scale_input <- data.frame(x = numeric())
    app.data$pixel_length_um <- NULL
  })

  #write pixel scale to menu bar...
  output$scale_text <- renderText({
    if(!is.null(app.data$pixel_length_um)){
      paste('~~~~~~ One pixel = ', round(as.numeric(app.data$pixel_length_um), 3), ' µm', sep = '')
    }
  })

  #Highlight cells image - first tab
  output$img1 <- renderDisplay({ #EBImage derived from Imager product
    if(!is.null(app.data$final_total_blobs_imager)){
      #print('test')
      my.img <- isolate(app.data$original_img) #EBImage, 3 channel RGB
      my.blobs <- isolate(app.data$final_total_blobs_imager) %>% as.array
      my.blobs <- adrop(my.blobs, drop = c(3,4)) %>% Image
      my.blobs <- my.blobs %>% bwlabel #Imager --> EBImage
      paintObjects(my.blobs, my.img[,,1:3], col='#ff00ff', #for EBImage product
                   thick = T, opac = c(0.5)) %>%  EBImage::display()
    }

  })

  #Highlight cells image - second tab
  output$img2 <- renderDisplay({ #EBImage derived from Imager product
    input$fileIn
    if(!is.null(app.data$final_total_blobs_imager)){
      my.img <- app.data$original_img #EBImage, 3 channel RGB
      my.blobs <- app.data$final_total_blobs_imager %>% as.array
      my.blobs <- adrop(my.blobs, drop = c(3,4)) %>% Image
      my.blobs <- my.blobs %>% bwlabel #Imager --> EBImage
      paintObjects(my.blobs, my.img[,,1:3], col='#ff00ff', #for EBImage product
                   thick = T, opac = c(0.5)) %>%  EBImage::display()
    }
  })


  #capture particle coordinates and IDs in the table, identify which blob they correspond to
  #paste ID into final_fea_shape table - now called app.data$all_data
  observeEvent(input$pixelPosition, {
    if(input$rec == 'Y'){
      current.blob <- app.data$final_total_blobs_ebimage_single_img[input$pixelPosition[1],input$pixelPosition[2]] %>% as.numeric
      new.row <- data.frame(id = input$id,
                            x = input$pixelPosition[1],
                            y = input$pixelPosition[2],
                            blob = current.blob)
      app.data$cell_data <- rbind(app.data$cell_data, new.row)
      app.data$all_data[current.blob, 'blob'] <- current.blob
      app.data$all_data[current.blob, 'id'] <- input$id
    }
  })

  #render a new plot that shows which cells have been identified
  output$img3 <- renderPlot({ #EBImage derived from Imager product
    if(!is.null(app.data$final_total_blobs_imager)){
      img <- app.data$original_img
      blobs <- app.data$final_total_blobs_ebimage_single_img
      blob_info <- app.data$all_data
      labels = blob_info$id
      x <- app.data$final_fea_moment$m.cx
      y <- app.data$final_fea_moment$m.cy
      res <-paintObjects(blobs, img[,,1:3], col=c('#ff00ff'), #for EBImage product
                         thick = T, opac = c(0.5))
      plot(res)
      text(x = x, y = y, label = labels, adj = c(0,2), col = 'white')
      #paintObjects(blobs, img[,,1:3], col=c('#ff00ff'), #for EBImage product
      #            thick = T, opac = c(0.5)) %>%  display(method = 'browser') %>% text(x = x, y = y, label = labels, adj = c(0,2), col = 'white')
    }
  })


  #produce summary info to be displayed in a third tab
  observeEvent(input$comp, {

    #catch for no cells present
    if(nrow(app.data$cell_data) == 0){
      dat <- app.data$all_data[1,]
      dat$species <- 'no cells'
      dat$blob <- 'no cells'
      dat$id <- 'no cells'
      dat$filament_length <- 0
      dat$biovolume_um3 <- 0

      total.cells <- 0
      total.cells.per.ml <- 0
    } else {

    #produce a data table
    dat <- app.data$all_data %>% as.data.table() #use the shape info (in pixel values) as the starting dataframe

    #normalise this data from pixels to um scale **have not normalised area or perimeter**
    dat$s.radius.mean <- dat$s.radius.mean * app.data$pixel_length_um
    dat$s.radius.sd <- dat$s.radius.sd * app.data$pixel_length_um
    dat$s.radius.min <- dat$s.radius.min * app.data$pixel_length_um
    dat$s.radius.max <- dat$s.radius.max * app.data$pixel_length_um

    #extract species name as first character
    dat$species <- substr(dat$id, 0,1)
    #extract any digits within character string (i.e. cell numbers of each filament)
    dat$filament_length <- str_extract(dat$id, "\\d") %>% as.numeric
    #input fil length of 1 cell is snow or blue green algae selected
    dat[, filament_length := ifelse(dat$species == 's' | dat$species == 'b', 1, dat$filament_length)]
    #calclualte pixel based biovolume (pi / 4 * d2 * h) from x2 radius.min (d), and x2 radius.max (h) / filament_length as
    ###****Rob to manually check this automated biovolume approximation using images in imageJ ***###
    dat[, biovolume_um3 := (pi / 4) * ((s.radius.min*2) * (s.radius.min*2)) * ((s.radius.max*2) / filament_length)]

    #total cells counted in the 1 large square
    total.cells <- sum(dat$filament_length, na.rm = T)
    total.cells.per.ml <- sum(dat$filament_length, na.rm = T) / 0.002
    }

    #produce a summary table per species recorded (from a NA removed dat...)
    dat.clean <- na.omit(dat)
    lev <- dat.clean$species %>% factor %>% levels
    if(length(lev) > 1){
      sum <- dat.clean[,.(total.count = sum(filament_length, na.rm = T), cells.per.ml = sum(filament_length, na.rm = T)/0.002,
                          #average filament length per sp
                          fil.av = mean(filament_length, na.rm = T), fil.sd = sd(filament_length, na.rm = T),
                          #relative species proportions
                          sp.proportions = sum(filament_length, na.rm = T) / total.cells * 100,
                          #biovolume calculated from pixel radius.min and radius.max (divided by filament length!)
                          biovolume.um3.av = mean(biovolume_um3, na.rm = T), biovolume.um3.sd = sd(biovolume_um3, na.rm = T)
      ), by = species]# group by species. - dont think this works when just one species- revert to normal code...
    } else {

      sum <- data.frame(species = as.character(lev),
                        total.count = sum(dat.clean$filament_length),
                        cells.per.ml = sum(dat.clean$filament_length) / 0.002,
                        fil.av = mean(dat.clean$filament_length),
                        fil.sd = sd(dat.clean$filament_length),
                        sp.proportions = sum(dat.clean$filament_length) / total.cells * 100,
                        biovolume.um3.av = mean(dat.clean$biovolume_um3),
                        biovolume.um3.sd = sd(dat.clean$biovolume_um3))
    }

    #push to app.data
    app.data$all_data <- dat
    app.data$summary_data <- sum
    app.data$total.cells.per.ml <- total.cells.per.ml
    app.data$total.cells.counted <- total.cells

  })

  #render summary table
  output$summary_table <- renderDataTable({
    if(nrow(app.data$cell_data) > 0 ){

    f <- app.data$summary_data
    f[,2:ncol(f)] <- round(f[,2:ncol(f)], 2)
    f
    }
  })

  #render species proportions bar chart
  output$sp_abundance <- renderPlot({
    if(nrow(app.data$cell_data) > 0){
    f <- app.data$summary_data %>% as.data.frame
    my.cols <- brewer.pal(nrow(f), 'Dark2')
    par(mar=c(4,3,1,1), mgp=c(1.8,0.6,0), las=1, tck=-0.01, lwd = 2, lty = 1, oma = c(0,0,0,0))
    p1 <- barplot(as.matrix(f$sp.proportions), beside = F, horiz = T, xlab = 'Relative Abundance', ylab = '', col =my.cols, lwd = 2)
    legend('top', legend = f$species, fill = my.cols, horiz = T, bty = 'n', inset = -0.03, xpd = T)
    }
  })

  #render output summary_text for summary panel
  output$summary_text <- renderText({
    paste0('Sample name:', app.data$file_name,'\nTotal cells counted:', app.data$total.cells.counted,'\nTotal cells per ml:', app.data$total.cells.per.ml )
  })

  #get working directory to save files in when 'wd' button clicked
  observeEvent(input$wd,{
    shinyDirChoose(input, 'wd',roots=c(wd='~/Dropbox/'))
    app.data$file_path <- parseDirPath(roots=c(wd='~/Dropbox/'), input$wd)
  })

  #export essentials once export button clicked and close the app
  observeEvent(input$export,{
    if(!is.null(app.data$file_path)){

      #get file names and wd path
      samp.nam <- app.data$file_name
      path <- app.data$file_path

      #write.csv of all data
      write.csv(app.data$all_data, file = paste(as.character(path), '/', paste(samp.nam, '_all_data.csv', sep=''), sep=''), row.names =F)

      #write.csv of summary data
      write.csv(app.data$summary_data, file = paste(as.character(path), '/', paste(samp.nam, '_summary_data.csv', sep=''), sep=''), row.names =F)

      #writeout out EBImage Image and associated blob layer
      writeImage(app.data$final_total_blobs_ebimage_single_img, files = paste(as.character(path), '/', paste(samp.nam, '_blob_image.png', sep=''), sep=''),
                 type = 'png')

      #write out annotated image from summary pane as pdf
      pdf(file = paste(as.character(path), '/', paste(samp.nam, '_summary_image.pdf', sep=''), sep=''))
      img <- app.data$original_img
      blobs <- app.data$final_total_blobs_ebimage_single_img
      blob_info <- app.data$all_data
      labels = blob_info$id
      x <- app.data$final_fea_moment$m.cx
      y <- app.data$final_fea_moment$m.cy
      res <-paintObjects(blobs, img[,,1:3], col=c('#ff00ff'), #for EBImage product
                         thick = T, opac = c(0.5))
      plot(res)
      text(x = x, y = y, label = labels, adj = c(0,2), col = 'white')
      dev.off()

      #write out settings used for cell segmentation
      my.set <- data.frame(sample = samp.nam,
                           mag_threshold = input$thresBar,
                           grow_threshold = input$growBar,
                           shrink_threshold = input$shrinkBar,
                           scale = app.data$pixel_length_um)
      write.csv(my.set, file = paste(as.character(path), '/', paste(samp.nam, '_all_settings.csv', sep=''), sep=''), row.names =F)


    }
  })



}

shinyApp(ui, server)

#https://stackoverflow.com/questions/55616335/identify-which-frame-is-currently-displayed-in-ebimage-shiny-display

