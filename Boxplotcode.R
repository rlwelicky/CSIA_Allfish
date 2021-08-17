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


nord.aur<-nord(palette="aurora", 5) #assigning colors from the nord package
nord.fro<-nord(palette="frost", 1) # adding blue to the aurora package from the frost package
col <- c(nord.aur, nord.fro)

pdf(file="Samplesizefigure.pdf", width=10, height=5)

mat <- matrix(data = c(1,2,3, 4,5), nrow=1) #the multi-panel component of the map is assigned using 
                                                          #layout, you could use par instead. Right now it is set up for a 3X2

par(oma=c(1,0,0,0),mar=c(6, 3, 2, 2), family = 'sans')  #adjusting margins
lay <- layout(mat=mat,  heights=c(1,1), widths=c(1,1,1)) #this assigns the width/height of each panel if you want to adjust it
#layout.show(6) check what your layout looks like



library(png)
solepic <- readPNG("solepic.png") #this imports the silhouettes as a png
hakepic<- readPNG("hakepic.png")
pollockpic <- readPNG("pollockpic.png")
herringpic <- readPNG("herringpic.png")
rockpic <- readPNG("rockpic.png")

x2<- 7

barplot(solen$fishcount, ylab = "", xlab="", col=col[5], 
        names.arg = solen$decade, ylim = c(0,12), main = "English Sole", cex.lab=1.5, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(solepic, 0.5,11.5,x2,10.5, interpolate = TRUE) #this adds the png to the plot -- you may need to play around with the numbers because it can skew teh dimensions depending on the png shape

barplot(rockn$fishcount, ylab = "", xlab="", col= col[6], 
        names.arg = rockn$decade,ylim = c(0,12), main = "Rockfish", cex.lab=1.5, 
        cex.axis = 1.25, cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(rockpic, 0.5,11.5,x2,10.5, interpolate = TRUE)


barplot(haken$fishcount, ylab = "", xlab="", col= col[4], 
        names.arg = haken$decade, ylim = c(0,12), main = "Pacific Hake", cex.lab=1, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(hakepic, 0.5,11.5,x2,10.5, interpolate = TRUE)

mtext('Decade', side = 1, line = 5, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.1)

barplot(pollockn$fishcount, ylab = "", xlab="", col= col[3], 
         names.arg = pollockn$decade, ylim = c(0,12), main = "Walleye Pollock", cex.lab=1, cex.axis = 1.25, 
        cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(pollockpic, 0.5,11.5,x2,10.5, interpolate = TRUE)


barplot(herringn$fishcount, ylab = "", xlab="", col= col[2],
        names.arg = herringn$decade, ylim = c(0,14), main = "Pacific Herring", cex.lab=1.5, 
        cex.axis = 1.25, cex.names = 1.25, cex.main=1.25, las=2)
rasterImage(herringpic, 0.5,12.5,x2,13.25, interpolate = TRUE)


#barplot(EG$Total, ylab = "", xlab="", col= col[1], 
        #ylim = c(0,20), main = "F. Eastern Stock", cex.lab=1, cex.axis = 1.25, 
        #cex.names = 1.25, cex.main=1.25, las=2, names.arg = EG$Decade)
#rasterImage(sl.eg, 3,13,5.5,18, interpolate = TRUE)

dev.off()

