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
library(EBImage)
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
                                             99, 100, value = 99.5, step = 0.1),
                                 sliderInput('growBar', '3. Grow Threshold',
                                             0, 20, value = 10, step = 1),
                                 sliderInput('shrinkBar', '4. Shrink Threshold',
                                             0, 20, value = 5, step = 1),
                                 uiOutput('sizeFraction'),
                                 textInput('scale_input', label = '6. Input x1,x2 values', width = '200px'),
                                 actionButton('scale', '7. Calc Scale', width = '200px'),
                                 textOutput('scale_text', inline = T)
                                 ),
                conditionalPanel("input.sidebarmenu === 'cellId'",
                                 dataTableOutput('cells_data')
                )
                )
  ),#end of dashboardSidebar

  dashboardBody(
    tabItems(
      tabItem(
        tabName = 'segmentCells',
        displayOutput('img1', width = '100%', height = '800px')
        # bsModal('modalExample', 'Set scale by drawing transect across one small square', 'scale', size = 'large',
        #         leafletOutput("leafmap", height = '600px'))
      ),
      tabItem(tabName = 'cellId',
              fluidRow(
                column(4,
                       #displayOutput('img2',width = '100%', height = '300px')
                       #leafletOutput('img2', width = '100%', height = '200px')
                       plotOutput('img2', height = '200px')
                       ),
                column(8,
                       actionButton('back', 'Back', width = '100%'),
                       actionButton('forward', 'Forward', width = '100%'),
                       textOutput('speciesId'),
                       textOutput('currentFrame'),
                       radioButtons('id', '', choices = c(
                         'Non algal', 'm1', 'm2', 'm3', 'm4',
                         'a1', 'a2', 'a3', 'a4', 'a5', 'a6',' a7', 'a8',
                         'a9', 'a10', 'a11', 'a12', 'c1','c2','c3','c4','snow', 'bluegreen','other'
                       ), inline = T, selected = 'Non algal')
                       )
              ),
              br(),
              fluidRow(
                column(12,
                       plotOutput('img3',width = '100%', height = '500px')) # this is slow!!
              ),
              fluidRow(
                column(12,
                       actionButton('comp', 'Compute Summary', width = '100%'))
              )
              ),
      tabItem(tabName = 'summary',
              fluidRow(
                column(12,
                       renderText('Summary statistics'),
                       verbatimTextOutput('summary_text')
                       ),
                column(12,
                       plotOutput('sp_abundance', 'Relative species abundance', width = '100%', height = "200px")
                       ),
                column(12,
                       dataTableOutput('summary_table', width = '100%')),
                column(2,
                       shinyDirButton('wd', 'Select WD', title = 'Select WD')
                       ),
                column(2,
                       actionButton('export', 'Export Data'))
              ))
    )#tabItems
  )#dashboardBody

)#end of UI


server <- function(input, output, session) {

  app.data <<- reactiveValues(data = NULL)

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

  #dynamically threshold for algae based on inputs
  observe({
    app.data$grid_removed_img
    if(!is.null(app.data$grid_removed_img)){
      grad <- imgradient(app.data$grid_removed_img, "xy")
      new.mag <- with(grad, sqrt(x^2+y^2))
      my.thres <- input$thresBar %>% as.numeric
      my.thres <- paste(my.thres, '%', sep='')
      new.thr <- threshold(new.mag, thr = my.thres,approx = F) #dynamic
      new.thr <- clean(new.thr, 2)
      grow.thres <- input$growBar
      new.thr <- grow(new.thr, as.numeric(input$growBar)) #dynamic
      new.thr <- shrink(new.thr, as.numeric(input$shrinkBar)) #dynamic
      app.data$initial_total_blobs_imager <- new.thr #Imager - single pixel set of total blobs (previously my_blobs)


      #extract blobs as EBImage product
      segs <- new.thr %>% as.cimg %>% as.array
      segs <- adrop(segs, drop = c(3,4)) %>% Image #EBImage
      segs <- segs %>% bwlabel #EBImage
      all_blobs <- stackObjects(segs, app.data$original_img) #EBImage

      #calculate summary statistics per blob
      #pixel intensities (mean, sd, mad, quantiles)
      app.data$initial_fea_int <- computeFeatures.basic(x = segs, app.data$original_img) %>% as.data.frame
      #shape features (area, perimeter, mean radius, sd, max radius, min radius) - all in pixel values
      app.data$initial_fea_shape <- computeFeatures.shape(x = segs, app.data$original_img) %>% as.data.frame
      #location features (center mass x/y, elliptical fit of major axis, elliptical eccentricity, object angle)
      app.data$initial_fea_moment <- computeFeatures.moment(x = segs, app.data$original_img) %>% as.data.frame
      #haralick texture parameters
      app.data$initial_fea_haralick <- computeFeatures.haralick(x = segs, app.data$original_img) %>% as.data.frame

      #identify which

      #push all to app.data
      app.data$initial_total_blobs_ebimage_single_img <- segs #EBImage; single image with all blobs (previously total_blobs)
      #app.data$all_blobs <- all_blobs #EBImage; image stack of separated blobs - dont need this until later...

    }

  })

  output$sizeFraction <- renderUI({
    if(!is.null(app.data$initial_fea_shape)){
      sliderInput('sizeRange', '5. Size Threshold',
                  min(app.data$initial_fea_shape$s.area),
                  (max(app.data$initial_fea_shape$s.area)/8),
                  step = 10,
                  value = min(app.data$initial_fea_shape$s.area))
    } else {
     # 'Waiting for image to load'
    }
  })

  #dynamically reduce blobs based on area threshold
  observeEvent( input$sizeRange, {
    if(!is.null(app.data$initial_fea_shape)){

    #get area threshold
    area.thres <- input$sizeRange
    app.data$area_thres <- area.thres

    #get range of areas (from EBImage initial_fea_shape)
    initial.areas <- app.data$initial_fea_shape$s.area

    #which to drop?
    to.rm <- which(initial.areas < area.thres)

    #remove from initial pix set to create new pixset for subsequent analysis
    f <- app.data$initial_total_blobs_imager # want to remove from this object
    lab <- label(f) #ImageR image labelled like bwlabel in EBImage
    if(length(to.rm) > 0){
    for(i in 1:length(to.rm)){
      lab[lab == to.rm[i]] <- 0
    }
    }
    mod.f <- lab > 0
    app.data$final_total_blobs_imager <- mod.f

    #make subsequent objects and push to app.data
    #extract blobs as EBImage product
    segs.mod <- mod.f %>% as.cimg %>% as.array
    segs.mod <- adrop(segs.mod, drop = c(3,4)) %>% Image #EBImage
    segs.mod <- segs.mod %>% bwlabel #EBImage
    all_blobs_mod <<- stackObjects(segs.mod, app.data$original_img) #EBImage

    #calculate summary statistics per blob
    #pixel intensities (mean, sd, mad, quantiles)
    app.data$final_fea_int <- computeFeatures.basic(x = segs.mod, app.data$original_img) %>% as.data.frame
    #shape features (area, perimeter, mean radius, sd, max radius, min radius) - all in pixel values
    app.data$final_fea_shape <-computeFeatures.shape(x = segs.mod, app.data$original_img) %>% as.data.frame
    #location features (center mass x/y, elliptical fit of major axis, elliptical eccentricity, object angle)
    app.data$final_fea_moment <- computeFeatures.moment(x = segs.mod, app.data$original_img) %>% as.data.frame
    #haralick texture parameters
    app.data$final_fea_haralick <- computeFeatures.haralick(x = segs.mod, app.data$original_img) %>% as.data.frame

    #push all to app.data
    app.data$final_total_blobs_ebimage_single_img <- segs.mod #EBImage; single image with all blobs (previously total_blobs)
    app.data$final_total_blobs_ebimage_stacked_img <- all_blobs_mod #EBImage; image stack of separated blobs - dont need this until later...

    }
  })


  #render pop-up plot of raster image for drawing of transect to get pixel scale
  # output$leafmap <- renderLeaflet({
  #   my.rast <- getFrame(app.data$original_img, i = 1, type='render')[,,1] %>% #EBImage --> raster
  #     as.array %>% raster
  #   extent(my.rast) <- c(0, dim(my.rast)[1], 0, dim(my.rast)[2])
  #   proj4string(my.rast) <- CRS(leafletProj)
  #   leaflet() %>%
  #     addRasterImage(my.rast, project = T) %>%
  #     addDrawToolbar( polylineOptions = T, circleOptions = F, rectangleOptions = F,
  #                     markerOptions = F, circleMarkerOptions = F, singleFeature = T,
  #                     polygonOptions = drawPolygonOptions(repeatMode = F), editOptions = editToolbarOptions()
  #     )
  # })

  #calculate pixel scale in um from input and write to app.data list
  # observeEvent(input$leafmap_draw_new_feature, {
  #   input$leafmap_draw_new_feature
  #   poly.coords <- input$leafmap_draw_new_feature$geometry$coordinates %>% unlist
  #   poly.mat <- data.frame(x = poly.coords[c(TRUE, FALSE)], y = poly.coords[c(FALSE, TRUE)])
  #   xy <- SpatialPoints(poly.mat)
  #   proj4string(xy) <- poly.proj
  #   xy <- as.data.frame(spTransform(xy, leafletProj))
  #   pixel.length <- diff(xy$x)
  #   pixel.um <- 250 / pixel.length #*********THIS ASSUMES 1 SMALL SQUARE is used for the scale-setting, i.e. 250 um in length!!******#######
  #   app.data$pixel_length_um <- pixel.um
  #
  # })

  #calculate pixel scale in um from text input$scale_input when input$scale action button is pressed
  observeEvent(input$scale, {
    my.text <- input$scale_input
    x1 <- my.text %>% as.character %>% strsplit(., split=',') %>% lapply('[[',1) %>% unlist %>% as.numeric
    x2 <- my.text %>% as.character %>% strsplit(., split=',') %>% lapply('[[',2) %>% unlist %>% as.numeric
    pixel.length <- diff(c(x1, x2))
    pixel.um <- 250 / pixel.length #*********THIS ASSUMES 1 SMALL SQUARE is used for the scale-setting, i.e. 250 um in length!!******#######
    app.data$pixel_length_um <- pixel.um
  })

  #write pixel scale to menu bar...
   output$scale_text <- renderText({
     if(!is.null(app.data$pixel_length_um)){
     paste('~~~~~~ One pixel = ', round(as.numeric(app.data$pixel_length_um), 3), ' µm', sep = '')
     }
   })

  #Highlight cells image
  output$img1 <- renderDisplay({ #EBImage derived from Imager product
    input$fileIn
    input$sizeRange
    if(!is.null(app.data$final_total_blobs_imager)){
      my.img <- app.data$original_img #EBImage, 3 channel RGB
      my.blobs <- app.data$final_total_blobs_imager %>% as.array
      my.blobs <- adrop(my.blobs, drop = c(3,4)) %>% Image
      my.blobs <- my.blobs %>% bwlabel #Imager --> EBImage
      paintObjects(my.blobs, my.img[,,1:3], col='#ff00ff', #for EBImage product
                   thick = T, opac = c(0.5)) %>%  display()
    }
  })

  #set up raster stack for plotting indivivudal blobs in cellId tab
  observe({
    if(!is.null(app.data$final_total_blobs_ebimage_stacked_img)){
    f <- app.data$final_total_blobs_ebimage_stacked_img #EBImage; image stack of separated blobs
    app.data$cells_data <- data.frame(blob = 1:dim(f)[4], id = NA)
    x <- stack()
    for(i in 1:(dim(f)[4])){
      h <- getFrame(f, i = i, type='render')[,,1] %>% #EBImage function --> raster
        as.array %>% raster
      names(h) <- i
      x <- stack(x, h)
    }
    app.data$raster_stack_blobs <- x
    app.data$current_layer <- 1
    }

  })

  observeEvent(input$forward, {
    if(app.data$current_layer < dim(app.data$final_total_blobs_ebimage_stacked_img)[4]){
      app.data$current_layer <- app.data$current_layer + 1
    }
  })

  # observeEvent(input$id, {
  #   if(app.data$current_layer < dim(app.data$final_total_blobs_ebimage_stacked_img)[4]){
  #     app.data$current_layer <- app.data$current_layer + 1
  #   }
  # })

  observeEvent(input$back,{
    if(app.data$current_layer > 1){
      app.data$current_layer <- app.data$current_layer - 1
    }
  })

  output$img2 <- renderPlot({
    if(input$sidebarmenu == 'cellId'){
      current_layer_to_plot <- raster(app.data$raster_stack_blobs,
                                      layer = app.data$current_layer)
      #ff <- app.data$final_total_blobs_ebimage_stacked_img
      #current.blob <- app.data$current_layer

      par(mar=c(0,0,0,0), oma=c(0,0,0,0))
      plot(current_layer_to_plot, legend = F)
      #plot(ff, current.blob)
    }
  })

  output$currentFrame <- renderText({
    paste('Current Frame: ', app.data$current_layer, ' of ',
          dim(app.data$final_total_blobs_ebimage_stacked_img)[4], sep='')
  })

  output$speciesId <- renderText({
    paste('Current ID:', input$id)
  })

  #render main display in cell ID showing only cell of interest highlighted
  #** this works but is VERY slow to render!!**
  output$img3 <- renderPlot({
    if(!is.null(app.data$final_total_blobs_imager)){

      current.blob <- app.data$current_layer

    #the below works but is still a bit slow....how about as rasters.
     img <- app.data$imager_img %>% grayscale #much quicker if plotted in grayscale!
     tot.blobs <- isolate(app.data$final_total_blobs_imager)
     sp <- split_connected(tot.blobs)
     my_xy <<- where(sp[[current.blob]]) %>% dplyr::summarise(mx=mean(x),my=mean(y)) %>% round() #added round
     my_x <<- my_xy[1] %>% as.numeric
     my_y <<- my_xy[2] %>% as.numeric
     my_xlim <<- c((my_x - dim(img)[1]* 0.1), (my_x + dim(img)[1]* 0.1))
     #if(my_xlim[1] < 0){my_xlim[1] <- 0}
     #if(my_xlim[2] > dim(img)[1]){my_xlim[1] <- dim(img)[1]}
     my_ylim <<- c((my_y + dim(img)[2]* 0.1),(my_y - dim(img)[2]* 0.1))
     #if(my_ylim[1] > dim(img)[2]){ my_ylim[1] <- dim(img)[2]}
     #if(my_ylim[2] < 0){my_ylim[2] <- 0}

     par(mar=c(0,0,0,0))
     plot(img, xlim=my_xlim, ylim=my_ylim, xaxt='n', yaxt='n') #plotting imager product
     imager::highlight(sp[[current.blob]], lwd = 3) #imager product

    }

  })

  #capture currently selected species id in cells_data dataframe
  observeEvent(c(input$id, input$back, input$forward), {
    if(!is.null(app.data$current_layer)){
    input$back
    input$forward
    app.data$cells_data[app.data$current_layer, 'id'] <- input$id

    }
  })

  #render the dataframe to view
  output$cells_data <- renderDataTable({
    if(!is.null(app.data$current_layer)){
      datatable(app.data$cells_data, rownames = F,
                filter = 'none', options = list(
        autoWidth = TRUE, dom = 't', paging = F)) %>%
          formatStyle(c('blob','id'),  color = 'red',
                      backgroundColor = 'orange', fontWeight = 'bold')

      }
  })

  #update radio labels if id already assigned or if skipping back/forward
  observeEvent(c(input$back, input$forward), {

    if(!is.null(app.data$current_layer)){

    if(!is.na(app.data$cells_data$id[app.data$current_layer])){
    updateRadioButtons(session, inputId = 'id', label = '', choices = c(
      'Non algal', 'm1', 'm2', 'm3', 'm4',
      'a1', 'a2', 'a3', 'a4', 'a5', 'a6',' a7', 'a8',
      'a9', 'a10', 'a11', 'a12', 'c1','c2','c3','c4','snow', 'bluegreen','other'
    ), inline = T, selected = app.data$cells_data$id[app.data$current_layer])
    } else {
      updateRadioButtons(session, inputId = 'id', label = '', choices = c(
        'Non algal', 'm1', 'm2', 'm3', 'm4',
        'a1', 'a2', 'a3', 'a4', 'a5', 'a6',' a7', 'a8',
        'a9', 'a10', 'a11', 'a12', 'c1','c2','c3','c4','snow', 'bluegreen','other'
      ), inline = T, selected = 'Non algal')
    }
    }
  })

  #capture data table species IDs and images of cells and their delination on the main image (might want to build CNN with total main image!)
  #produce summary info to be displayed in a third tab
  observeEvent(input$comp, {

    #name raster slices in app.data$raster_stack_blobs - this also numbers their occurrences by default..
    names(app.data$raster_stack_blobs) <- app.data$cells_data$id

    #produce a data table
    dat <- app.data$final_fea_shape %>% as.data.table() #use the shape info (in pixel values) as the starting dataframe

    #normalise this data from pixels to um scale **have not normalised area or perimeter**
    dat$s.radius.mean <- dat$s.radius.mean * app.data$pixel_length_um
    dat$s.radius.sd <- dat$s.radius.sd * app.data$pixel_length_um
    dat$s.radius.min <- dat$s.radius.min * app.data$pixel_length_um
    dat$s.radius.max <- dat$s.radius.max * app.data$pixel_length_um

    #get particle number (essentially same as row number)
    dat$particle <- app.data$cells_data$blob
    #get user input ID
    dat$id <- app.data$cells_data$id
    #remove non algal particles
    dat <- dat[-which(dat$id == 'Non algal'),]
    #extract species name as first character
    dat$species <- substr(dat$id, 0,1)
    #extract any digits within character string (i.e. cell numbers of each filament)
    dat$filament_length <- str_extract(dat$id, "\\d") %>% as.numeric
    #input fil length of 1 cell is snow or blue green algae selected
    dat[, filament_length := ifelse(dat$species == 's' | dat$species == 'b', 1, dat$filament_length)]
    #input fil length of 0, i.e. no cells, if 'other' was selected as cell ID
    #dat[, filament_length := ifelse(dat$species == 'o', 0, dat$filament_length)]
    #calclualte pixel based biovolume (pi / 4 * d2 * h) from x2 radius.min (d), and x2 radius.max (h) / filament_length as
    ###****Rob to manually check this automated biovolume approximation using images in imageJ ***###
    dat[, biovolume_um3 := (pi / 4) * ((s.radius.min*2) * (s.radius.min*2)) * ((s.radius.max*2) / filament_length)]

    #total cells counted in the 1 large square
    total.cells <- sum(dat$filament_length, na.rm = T)
    total.cells.per.ml <- sum(dat$filament_length) / 0.002

    #produce a summary table per species recorded
    #haemo is usually 9 squares of 1mm by 1mm with 0.2 mm depth, which gives 0.002 cm3 volume, i.e. 0.002 ml so divide cells in 1 large sq. by this = cells per ml
    sum <- dat[,.(total.count = sum(filament_length, na.rm = T), cells.per.ml = sum(filament_length, na.rm = T)/0.002,
                  #average filament length per sp
                  fil.av = mean(filament_length, na.rm = T), fil.sd = sd(filament_length, na.rm = T),
                  #relative species proportions
                  sp.proportions = sum(filament_length, na.rm = T) / total.cells * 100,
                  #biovolume calculated from pixel radius.min and radius.max (divided by filament length!)
                  biovolume.um3.av = mean(biovolume_um3, na.rm = T), biovolume.um3.sd = sd(biovolume_um3, na.rm = T)
                  ), dat$species]# group by species.

    app.data$all_data <- dat
    app.data$summary_data <- sum
    app.data$total.cells.per.ml <- total.cells.per.ml
    app.data$total.cells.counted <- total.cells

  })

  #render summary table
  output$summary_table <- renderDataTable({
    f <- app.data$summary_data
    f[,2:ncol(f)] <- round(f[,2:ncol(f)], 2)
    f
  })

  #render species proportions bar chart
  output$sp_abundance <- renderPlot({
  f <- isolate(app.data$summary_data) %>% as.data.frame
  my.cols <- brewer.pal(nrow(f), 'Dark2')
  par(mar=c(4,3,1,1), mgp=c(1.8,0.6,0), las=1, tck=-0.01, lwd = 2, lty = 1, oma = c(0,0,0,0))
  p1 <- barplot(as.matrix(f$sp.proportions), beside = F, horiz = T, xlab = 'Relative Abundance', ylab = '', col =my.cols, lwd = 2)
  legend('top', legend = f$dat, fill = my.cols, horiz = T, bty = 'n', inset = -0.03, xpd = T)
  })

  #render output summary_text for summary panel
  output$summary_text <- renderText({
    paste0('Sample name:', app.data$file_name,'\nTotal cells counted:', app.data$total.cells.counted,'\nTotal cells per ml:', app.data$total.cells.per.ml )
  })

  #get working directory to save files in when 'wd' button clicked
  observeEvent(input$wd,{
    shinyDirChoose(input, 'wd',roots=c(wd='~/'))
    app.data$file_path <- parseDirPath(roots=c(wd='~/'), input$wd)
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

      #writeout out EBImage stack and associated identification

      #writeout annotated single image?

      #writeout settings used for cell segmentation.

    }
  })

  #rv <- reactiveValues(pixelPosition = NULL)

  #retrieve pixel x/y from clicks
  observeEvent(input$pixelPosition, {
    current.pixel <- input$pixelPosition
    current.species <- input$species

    print(current.pixel)
  })

}

shinyApp(ui, server)

#https://stackoverflow.com/questions/55616335/identify-which-frame-is-currently-displayed-in-ebimage-shiny-display
