#processing 2018 counts for Ted
library(magrittr)
library(data.table)

#################################################
#Pull all raw counts into one data frame, with a 'large square' counter column, and a sample name column (Should be 75 samples in total)
my.dirs <- list('~/Dropbox/2018 Greenland Samples/260718_surface_ice_samples_1-25/', '~/Dropbox/2018 Greenland Samples/260718_surface_ice_samples_2-1_2-25/', '~/Dropbox/2018 Greenland Samples/260718_surface_ICE_samples_26-50/')
#make example data frame to catch data
dat <- read.csv('~/Dropbox/2018 Greenland Samples/260718_surface_ice_samples_1-25/1/260718_Surf_Ice_sample_3_bak_all_data.csv')
dat <- dat[-(1:nrow(dat)),]
dat$large_sq <- numeric()
dat$sample <- character()

for(i in 1:length(my.dirs)){

  setwd(my.dirs[[i]])
  my.folders <- list.dirs(dir())

  for(j in 1:length(my.folders)){
    setwd(my.folders[[j]])
    my.files <- list.files(pattern = '_all_data.csv')

    tmp.dat <- dat[-(1:nrow(dat)),]
    tmp.dat$large_sq <- numeric()

    if(length(my.files > 0)){

    for(l in 1:length(my.files)){
      dd <- read.csv(my.files[l])
      dd <- na.omit(dd)
      dd$large_sq <- l

      tmp.dat <- rbind(tmp.dat, dd)

    }#my.files L

    tmp.dat$sample <- getwd()
    setwd(my.dirs[[i]])

    } else {
      tmp.dat[1, ] <- 0
      tmp.dat$sample <- getwd()
      setwd(my.dirs[[i]])

    }

    dat <- rbind(dat, tmp.dat)
  }#j loop folders
}# i loop directories

##################################################

#write out
write.csv(dat, file = '~/Dropbox/2018 Greenland Samples/all_75_all_data.csv', row.names = F)


##################################################
#Compute summary per sample

dat <- read.csv('~/Dropbox/2018 Greenland Samples/all_75_all_data.csv')
dat$sample <- factor(dat$sample)
lev <- levels(dat$sample)# should be 75

total.sum <- data.frame("species" = numeric() ,   "total.count"  = numeric() ,    "cells.per.ml"  = numeric() ,   "fil.av" = numeric() ,    "fil.sd"    = numeric() , "sp.proportions" = numeric() ,  "biovolume.um3.av" = numeric() , "biovolume.um3.sd"= numeric(), total.cells.per.ml = numeric(), sample = character() )

#for each sample, we now assume all counts were taken over NINE large squares, i.e. one total haemo grid...
for(i in 1:length(lev)){
  dat.clean <- subset(dat, sample == lev[i]) %>% as.data.table
  if(is.na(dat.clean$species)){dat.clean$species <- 'no cells'}
  total.cells <- sum(dat.clean$filament_length)
  total.cells.per.ml <- sum(dat.clean$filament_length) / 0.018 #this is 0.002 for 1 large sq. * 9

  sp.lev <- dat.clean$species %>% factor %>% levels

  if(length(sp.lev) > 1){
    sum <- dat.clean[,.(total.count = sum(filament_length, na.rm = T), cells.per.ml = sum(filament_length, na.rm = T)/0.018,
                        #average filament length per sp
                        fil.av = mean(filament_length, na.rm = T), fil.sd = sd(filament_length, na.rm = T),
                        #relative species proportions
                        sp.proportions = sum(filament_length, na.rm = T) / total.cells * 100,
                        #biovolume calculated from pixel radius.min and radius.max (divided by filament length!)
                        biovolume.um3.av = mean(biovolume_um3, na.rm = T), biovolume.um3.sd = sd(biovolume_um3, na.rm = T)
    ), by = species]# group by species. - dont think this works when just one species- revert to normal code...
  } else {

    sum <- data.frame(species = as.character(sp.lev),
                      total.count = sum(dat.clean$filament_length),
                      cells.per.ml = sum(dat.clean$filament_length) / 0.018,
                      fil.av = mean(dat.clean$filament_length),
                      fil.sd = sd(dat.clean$filament_length),
                      sp.proportions = sum(dat.clean$filament_length) / total.cells * 100,
                      biovolume.um3.av = mean(dat.clean$biovolume_um3),
                      biovolume.um3.sd = sd(dat.clean$biovolume_um3))
  }
  sum$total.cells.per.ml <- NA
  sum$total.cells.per.ml[1] <- total.cells.per.ml
  sum$sample <- lev[i]

  total.sum <- rbind(total.sum, sum)
}

write.csv(total.sum, file = '~/Dropbox/2018 Greenland Samples/all_75_sum_data.csv', row.names = F)


##################################################
#Cut down to just cells per ml per sample and do some summary for Ted
my.sum <- total.sum[!is.na(total.sum$total.cells.per.ml),] %>% as.data.frame

range(my.sum$total.cells.per.ml)
hist(my.sum$total.cells.per.ml, breaks = 75)

plot(my.sum$total.cells.per.ml, log = 'y')
median(my.sum$total.cells.per.ml)
mean(my.sum$total.cells.per.ml)

clean <- 625
clean.low <- 625 - 381
clean.high <- 625 + 381

light <- 4.73 * 10^3
light.low <- 4.73 * 10^3 - 2.57 * 10^3
light.high <- 4.73 * 10^3 + 2.57* 10^3

heavy <- 2.9 * 10^4
heavy.low <- 2.9 * 10^4 - 2.01 * 10^4
heavy.high <- 2.9 * 10^4 + 2.01 * 10^4

my.clean <- subset(my.sum, total.cells.per.ml < clean.high & total.cells.per.ml > clean.low)
nrow(my.clean) / nrow(my.sum) * 100 #22.6%
my.light <- subset(my.sum, total.cells.per.ml < light.high & total.cells.per.ml > light.low)
nrow(my.light) / nrow(my.sum) * 100 #18.8%
my.heavy <- subset(my.sum, total.cells.per.ml < heavy.high & total.cells.per.ml > heavy.low)
nrow(my.heavy) / nrow(my.sum) * 100 #29.3
#this sums to 70%

plot(x = 1:75, my.sum$total.cells.per.ml, ylab = 'cells.per.ml', xlab = 'Sample Index')
rect(0, clean.low, 75, clean.high, col='blue')
rect(0, light.low, 75, light.high, col = 'green')
rect(0, heavy.low, 75, heavy.high, col = 'red')
points(my.sum$total.cells.per.ml, pch = 19)

my.clean <- subset(my.sum, total.cells.per.ml < clean)
nrow(my.clean) / nrow(my.sum) * 100 #20%
my.light <- subset(my.sum, total.cells.per.ml < heavy & total.cells.per.ml > clean)
nrow(my.light) / nrow(my.sum) * 100 #72%
my.heavy <- subset(my.sum, total.cells.per.ml > heavy)
nrow(my.heavy) / nrow(my.sum) * 100 #8%

write.csv(my.sum, file = '~/Dropbox/2018 Greenland Samples/all_75_small_sum.csv', row.names = F)
