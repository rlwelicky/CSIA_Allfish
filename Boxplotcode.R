### Packages ###
#install.packages(c("maps","mapdata","maptools","mapproj"))
#install.packages("PBSmapping")
#install.packages("ggplot2")
#install.packages("digest") #sometimes, I have to quit R, delete the digest directory from the library folder, then reinstall
#devtools::install_github("dkahle/ggmap")
#install.packages("RgoogleMaps")
#install.packages('rgdal')
#install.packages("colorspace")

library(nord)       #color palette  
library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(rgdal)
library(colorspace)
library(png) #used to get teh slhouettes
library(inauguration)

summ<-read.csv("SLSumm.csv", header=TRUE)


pdf(file="Samplesizefigure1.pdf", width=8.5, height=11)


mat <- matrix(data = c(1,2,3,4,
                        5, 6,7), ncol = 7) #the multi-panel component of the map is assigned using 
                                                          #layout, you could use par instead. Right now it is set up for a 3X2

par(oma=c(2,2,2,2),mar=c(5, 4, 4, 2), family = 'sans')  #adjusting margins
lay <- layout(mat=mat,  heights=c(0.5,0.5,0.5,0.5,0.5), widths=c(70,70,70,70,70)) #this assigns the width/height of each panel if you want to adjust it
#layout.show(6) check what your layout looks like



pdf(file="Samplesizefigure.pdf", width=11, height=5)

mat <- matrix(data = c(1,2,3, 4,5), nrow=1) #the multi-panel component of the map is assigned using 
                                                          #layout, you could use par instead. Right now it is set up for a 3X2

par(oma=c(1,2,0,0),mar=c(6, 3, 2, 2), family = 'sans')  #adjusting margins
lay <- layout(mat=mat,  heights=c(1,1), widths=c(1,1,1,1,1)) #this assigns the width/height of each panel if you want to adjust it
#layout.show(6) check what your layout looks like




library(png)
solepic <- readPNG("solepic.png") #this imports the silhouettes as a png
hakepic<- readPNG("hakepic.png")
pollockpic <- readPNG("pollockpic.png")
herringpic <- readPNG("herringpic.png")
rockpic <- readPNG("rockpic.png")




barplot(solen$fishcount, ylab = "", xlab="", col="#5445b1", 
     names.arg = solen$decade, ylim = c(0,12), main = "English Sole", cex.lab=2, cex.axis = 2, 
      cex.names = 2, cex.main=2, las=2)
rasterImage(solepic, 0.5,12,3,10.5, interpolate = TRUE) #this adds the png to the plot -- you may need to play around with the numbers because it can skew teh dimensions depending on the png shape
mtext('Sample Size', side = 2, line = 0, outer = FALSE, at = NA,
      adj = NA, padj =0, cex = 1.1)

barplot(rockn$fishcount, ylab = "", xlab="", col= "#0077b6", 
     names.arg = rockn$decade,ylim = c(0,12), main = "Rockfish", cex.lab=2, 
      cex.axis = 2, cex.names = 2, cex.main=2, las=2)
rasterImage(rockpic, 0.5,12,3,10.5, interpolate = TRUE)


barplot(haken$fishcount, ylab = "", xlab="", col= "#f3c483", 
     names.arg = haken$decade, ylim = c(0,12), main = "Pacific Hake", cex.lab=1, cex.axis = 2, 
      cex.names = 2, cex.main=2, las=2)
rasterImage(hakepic, 0.5,12,3,10.5, interpolate = TRUE)
mtext('Decade', side = 1, line = 0, outer = FALSE, at = NA,
      adj = 0, padj = NA, cex = 1.1)


barplot(pollockn$fishcount, ylab = "", xlab="", col= "#5c1a33", 
        names.arg = pollockn$decade, ylim = c(0,12), main = "Walleye Pollock", cex.lab=1, cex.axis = 2, 
      cex.names = 2, cex.main=2, las=2)
rasterImage(pollockpic, 0.5,12,3,10.5, interpolate = TRUE)


barplot(herringn$fishcount, ylab = "", xlab="", col= "#cd3341",
       names.arg = herringn$decade, ylim = c(0,14), main = "Pacific Herring", cex.lab=2, 
      cex.axis = 2, cex.names = 2, cex.main=2, las=2)
rasterImage(herringpic, 0.5,12.5,3,12, interpolate = TRUE)



barplot(solen$fishcount, ylab = '', xlab="", col="#5445b1", 
        names.arg = solen$decade, ylim = c(0,12), main = "English Sole", cex.lab=1.5, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(solepic, 2,10.5,8.5,11.5, interpolate = TRUE) #this adds the png to the plot -- you may need to play around with the numbers because it can skew teh dimensions depending on the png shape
text(5.25,10, label = 'SL = 129 ± 20 (mm)', cex=1.1)
text(5.25,9.5, label = 'n = 42', cex=1.1)
mtext('Number of fish sampled', side = 2, line = 0, outer = TRUE, at = NA,
      adj = NA, padj = NA, cex = 1.1)


barplot(rockn$fishcount, ylab = "", xlab="", col= "#0077b6", 
        names.arg = rockn$decade,ylim = c(0,12), main = "Rockfish", cex.lab=1.5, 
        cex.axis = 1.25, cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(rockpic, 2.5,10.5,9.5,11.5, interpolate = TRUE)
text(6,10, label = 'SL = 191 ± 36 (mm)', cex=1.1)
text(6,9.5, label = 'n = 36', cex=1.1)

barplot(haken$fishcount, ylab = "", xlab="", col= "#f3c483", 
        names.arg = haken$decade, ylim = c(0,14), main = "Pacific Hake", cex.lab=1, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(hakepic,1,12.5,7.5,13.5, interpolate = TRUE)
text(4.25,12, label = 'SL = 180 ± 66 (mm)', cex=1.1)
text(4.25,11.5, label = 'n = 42', cex=1.1)
mtext('Decade', side = 1, line = 5, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.1)

barplot(pollockn$fishcount, ylab = "", xlab="", col= "#5c1a33", 
         names.arg = pollockn$decade, ylim = c(0,14), main = "Walleye Pollock", cex.lab=1, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(pollockpic, 2.5,12.5,9,13.5, interpolate = TRUE)
text(5.75,12, label = 'SL = 133 ± 56 (mm)', cex=1.1)
text(5.75,11.5, label = 'n = 50', cex=1.1)

barplot(herringn$fishcount, ylab = "", xlab="", col= "#cd3341",
        names.arg = herringn$decade, ylim = c(0,20), main = "Pacific Herring", cex.lab=1.5, 
        cex.axis = 1.25, cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(herringpic, 2.5,18.5,9,19.25, interpolate = TRUE)
text(5.75,17, label = 'SL = 141 ± 33 (mm)', cex=1.1)
text(5.75,16.25, label = 'n = 52', cex=1.1)


dev.off()

