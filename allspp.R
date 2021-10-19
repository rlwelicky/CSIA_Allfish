#lets input each spp file
sole<-read.csv("Sole_CSIA_withmetadata.csv", header=TRUE)
hake<-read.csv('Hake_CSIA_withmetadata.csv', header=TRUE) #one fish was manually removed because it was from the coast.
herring<-read.csv("Herring(with Pollock)_CSIA_withmetadata.csv")
herring$lat<-as.numeric(herring$lat)
is.numeric(herring$lat)
pollock<-read.csv("Pollock_CSIA_withmetadata.csv")
pollock$lag0_temp<-as.numeric(pollock$lag0_temp)
rock<-read.csv("Rock_CSIA_withmetadata.csv")
reruns<-read.csv("Reruns_CSIA_withmetadata.csv")
reruns$lat<-as.numeric(reruns$lat)             
library(tidyverse)
library(dplyr)
library(broom)
library(DHARMa)
library(spdep)
library(glmmTMB)
library(multcomp)
library(lmtest)
library(car)
library(ggplot2)
library(cowplot)
library(ggpubr)
#now lets combine the sp. files into one datatable
speciescombined<-bind_rows(hake,herring, pollock, sole,rock,reruns)
#I need to make a column of final values for each AA for indiv analysis. # I included this values_fn statement because if there are duplicate samples, we don't know which best, so we will take the average of the duplicate.
allspp<-speciescombined%>%
  filter(!is.na(c(finalvalue)))%>%
  pivot_wider(id_cols = c(hostid, hostsp, sl, tl, lat, long, year, decade, lag0_temp, lag1_temp, lag2_temp, lag3_temp, lag12_temp), names_from = AA, values_from = finalvalue)
duplicated(allspp$hostid)#there are a few duplicated but its ok, bc one of the duplicate doensn't have glu and phe so it won't be used anyways
#I need to make one variable called site using lat and long
allspp$latjitt<-jitter(allspp$lat, factor=0.1, amount=NULL)
allspp$longjitt<-jitter(allspp$long, factor=0.1, amount=0)
allspp$site=paste(allspp$lat,allspp$long, na.rm = TRUE)
#c need to reate a column of trophic position calculations

allspp$tp=(((allspp$glu-allspp$phe-3.4)/7.6)+1)
allspp$diff=(allspp$glu-allspp$phe)


a<-allspp%>%
  group_by(hostsp, na.rm = TRUE)%>% 
   summarize(
            averagediff =mean(diff,na.rm=TRUE),
            averagesl = mean(sl, na.rm = TRUE),
            sdsl = sd(sl, na.rm = TRUE),
            mediansl = median(sl, na.rm = TRUE),
            fishcount = n_distinct(hostid, na.rm = TRUE),
      
            averagetp = mean (tp, na.rm = TRUE))


#Are the data spatially autocorrelated?
#spatial checks for glutamic acid, this would be the same as phe or tp, since all these values are from 1 fish from the same site; One fish didn't have a lat/long, so I removed it to run the check using a unique dataframe for the check
#the Moran's I test cannot handle NAs, so first I'm removing all NAs from this 'mock' dataset
allsppspatial<- allspp %>%
  filter(!is.na(latjitt)) %>%
  filter(!is.na(longjitt))%>%
  filter(!is.na(tp)) %>%
  filter(!is.na(year))

spatial.allspp<-glm(tp~year, data=allsppspatial)
simspatial.allspp<-simulateResiduals(fittedModel = spatial.allspp)
spatialtest.allspp<-testSpatialAutocorrelation(simulationOutput = simspatial.allspp,  x = allsppspatial$longjitt, y = allsppspatial$latjitt)
spatialtest.allspp #Yay, no spatial corr


#temporal checks for glutamic acid, his would be the same as phe or tp, since all these values are from 1 fish from the same site

allspp$timegroup<-allspp$year
temporal.glu<-glm(glu~year, data= allspp)
temporaltest.glu<-dwtest(temporal.glu, order.by = NULL, alternative = "two.sided", exact = FALSE, tol = 1e-10) 
temporaltest.glu #yay, no temporal corr
dwtest(temporal.glu)




#Analyses for only Sole

sole_only<-allspp%>%
  filter(hostsp == "sole")


#what's the median sl needed for figures?
sole_only%>%
  summarise(
    mediansl = median(sl)
  )


#SOLE GLUTAMIC ACID
library(predictmeans)
normalityglu<-glm(glu~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalityglu, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole) #glu is normally distributed
glusole<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site)  , family = "gaussian", data= sole_only)
summary(glusole)



predict(glusole, sole_only, allow.new.levels=TRUE)
ndglusole<-sole_only[1,]
ndglusole$year<-"new"
ndglusole_pop<-data.frame(year=sole_only$year, site=NA, sl = 126) 
glusolepredict<-predict(glusole, newdata=ndglusole_pop, se.fit=TRUE)
as.data.frame(glusolepredict)
glusolepredict$year<-sole_only$year
glusolepredict<-as.data.frame(glusolepredict)

glusolefig<- ggplot() + geom_point(data = sole_only, aes(x = year, y = glu), color="#5445b1") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)")  +theme_classic()


#SOLE PHE

normalityphe<-glm(phe~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalityphe, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole)#phe is normally distrbuted

phesole<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= sole_only)
summary(phesole)



predict(phesole, sole_only, allow.new.levels=TRUE)
ndphesole<-sole_only[1,]
ndphesole$year<-"new"
ndphesole_pop<-data.frame(year=sole_only$year, site=NA, sl = 126) 
phesolepredict<-predict(phesole, newdata=ndphesole_pop, se.fit=TRUE)
as.data.frame(phesolepredict)
phesolepredict$year<-sole_only$year
phesolepredict<-as.data.frame(phesolepredict)

phesolefig<- ggplot() + geom_point(data = sole_only, aes(x = year, y = phe), color="#5445b1") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic()

#SOLE TP

normalitytp<-glm(tp~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalitytp, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole) #tp is normally distributed
tpsole<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= sole_only)
summary(tpsole)

predict(tpsole, sole_only, allow.new.levels=TRUE)
ndtpsole<-sole_only[1,]
ndtpsole$year<-"new"
ndtpsole_pop<-data.frame(year=sole_only$year, site=NA, sl = 126) 
tpsolepredict<-predict(tpsole, newdata=ndtpsole_pop, se.fit=TRUE)
as.data.frame(tpsolepredict)
tpsolepredict$year<-sole_only$year
tpsolepredict<-as.data.frame(tpsolepredict)

tpsolefig<- ggplot()  + geom_point(data = sole_only, aes(x = year, y = tp), color="#5445b1") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()





#Analyses for only hake
hake_only<-allspp%>%
  filter(hostsp == "hake")
#find median value first = 162

hake_only%>%
  summarise(
    mediansl = median(sl)
  )

#HAKE GLU
normalityglu<-glm(glu~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalityglu, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal :)
gluhake<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= hake_only)
summary(gluhake)



predict(gluhake, hake_only, allow.new.levels=TRUE)
ndgluhake<-hake_only[1,]
ndgluhake$year<-"new"
ndgluhake_pop<-data.frame(year=hake_only$year, site=NA, sl = 162) 
gluhakepredict<-predict(gluhake, newdata=ndgluhake_pop, se.fit=TRUE)
as.data.frame(gluhakepredict)
gluhakepredict$year<-hake_only$year
gluhakepredict<-as.data.frame(gluhakepredict)

gluhakefig<- ggplot()  + geom_point(data = hake_only, aes(x = year, y = glu), color="#f3c483") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic()

#HAKE PHE

normalityphe<-glm(phe~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalityphe, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal :) 
phehake<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(phehake)


predict(phehake, hake_only, allow.new.levels=TRUE)
ndphehake<-hake_only[1,]
ndphehake$year<-"new"
ndphehake_pop<-data.frame(year=hake_only$year, site=NA, sl = 126) 
phehakepredict<-predict(phehake, newdata=ndphehake_pop, se.fit=TRUE)
as.data.frame(phehakepredict)
phehakepredict$year<-hake_only$year
phehakepredict<-as.data.frame(phehakepredict)

phehakefig<- ggplot()  + geom_point(data = hake_only, aes(x = year, y = phe), color="#f3c483") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic()

#HAKE TP
normalitytp<-glm(tp~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalitytp, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal
tphake<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(tphake)


predict(tphake, hake_only, allow.new.levels=TRUE)
ndtphake<-hake_only[1,]
ndtphake$year<-"new"
ndtphake_pop<-data.frame(year=hake_only$year, site=NA, sl = 162) 
tphakepredict<-predict(tphake, newdata=ndtphake_pop, se.fit=TRUE)
as.data.frame(tphakepredict)
tphakepredict$year<-hake_only$year
tphakepredict<-as.data.frame(tphakepredict)

tphakefig<- ggplot() + geom_point(data = hake_only, aes(x = year, y = tp), color="#f3c483") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()

#Analyses for only pollock
pollock_only<-allspp%>%
  filter(hostsp == "pollock")

#find median value first = 108
pollock_only %>%
  summarise(
    mediansl = median(sl, na.rm = TRUE)
  )

#POLLOCK GLU
normalityglu<-glm(glu~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalityglu, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock) #normal
glupollock<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(glupollock)


predict(glupollock, pollock_only, allow.new.levels=TRUE)
ndglupollock<-pollock_only[1,]
ndglupollock$year<-"new"
ndglupollock_pop<-data.frame(year=pollock_only$year, site=NA, sl = 162) 
glupollockpredict<-predict(glupollock, newdata=ndglupollock_pop, se.fit=TRUE)
as.data.frame(glupollockpredict)
glupollockpredict$year<-pollock_only$year
glupollockpredict<-as.data.frame(glupollockpredict)

glupollockfig<- ggplot()  + geom_point(data = pollock_only, aes(x = year, y = glu), color="#5c1a33") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic()

#POLLOCK PHE

normalityphe<-glm(phe~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalityphe, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock)#normal
phepollock<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(phepollock)


predict(phepollock, pollock_only, allow.new.levels=TRUE)
ndphepollock<-pollock_only[1,]
ndphepollock$year<-"new"
ndphepollock_pop<-data.frame(year=pollock_only$year, site=NA, sl = 126) 
phepollockpredict<-predict(phepollock, newdata=ndphepollock_pop, se.fit=TRUE)
as.data.frame(phepollockpredict)
phepollockpredict$year<-pollock_only$year
phepollockpredict<-as.data.frame(phepollockpredict)

phepollockfig<- ggplot() +  geom_point(data = pollock_only, aes(x = year, y = phe), color="#5c1a33") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic()


#Polock TP

normalitytp<-glm(tp~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalitytp, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock)
tppollock<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(tppollock)



predict(tppollock, pollock_only, allow.new.levels=TRUE)
ndtppollock<-pollock_only[1,]
ndtppollock$year<-"new"
ndtppollock_pop<-data.frame(year=pollock_only$year, site=NA, sl = 108) 
tppollockpredict<-predict(tppollock, newdata=ndtppollock_pop, se.fit=TRUE)
as.data.frame(tppollockpredict)
tppollockpredict$year<-pollock_only$year
tppollockpredict<-as.data.frame(tppollockpredict)

tppollockfig<- ggplot() + geom_point(data = pollock_only, aes(x = year, y = tp), color="#5c1a33") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()

#Analyses for only herring
herring_only<-allspp%>%
  filter(hostsp == "herring")

#find median value first = 146
herring_only %>%
  summarise(
    mediansl = median(sl, na.rm = TRUE)
  )

#HERRING GLU
normalityglu<-glm(glu~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalityglu, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
gluherring<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(gluherring)


#I can't seem to get the herring figure to work! Not "Error in eval(substitute(expr), data, enclos = parent.frame()) : 
#NAs's causing problems, so I removed them to make the figures only!

herring_only2<-herring_only %>% 
  drop_na(site, year, sl) 


predict(gluherring, herring_only2, allow.new.levels=TRUE)
ndgluherring<-herring_only[1,]
ndgluherring$year<-"new"
ndgluherring_pop<-data.frame(year=herring_only2$year, site=NA, sl = 146) 
gluherringpredict<-predict(gluherring, newdata=ndgluherring_pop, se.fit=TRUE)
as.data.frame(gluherringpredict)
gluherringpredict$year<-herring_only2$year
gluherringpredict<-as.data.frame(gluherringpredict)

gluherringfig<- ggplot() + geom_point(data = herring_only, aes(x = year, y = glu), color="#cd3341") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic()

#Herring PHE


normalityphe<-glm(phe~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalityphe, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
pheherring<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(pheherring)


predict(pheherring, herring_only2, allow.new.levels=TRUE)
ndpheherring<-herring_only[1,]
ndpheherring$year<-"new"
ndpheherring_pop<-data.frame(year=herring_only2$year, site=NA, sl = 146) 
pheherringpredict<-predict(pheherring, newdata=ndpheherring_pop, se.fit=TRUE)
as.data.frame(pheherringpredict)
pheherringpredict$year<-herring_only2$year
pheherringpredict<-as.data.frame(pheherringpredict)

pheherringfig<- ggplot()  + geom_point(data = herring_only, aes(x = year, y = phe), color="#cd3341") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic()

#HERRING TP

normalitytp<-glm(tp~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalitytp, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
tpherring<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), na.action = na.omit, family = "gaussian", data= herring_only)
summary(tpherring)


library(performance)
check_collinearity(
  pheherring,
  component = c("all", "conditional", "count", "zi", "zero_inflated"),
  verbose = TRUE)
check_collinearity(
  gluherring,
  component = c("all", "conditional", "count", "zi", "zero_inflated"),
  verbose = TRUE)
check_collinearity(
  tpherring,
  component = c("all", "conditional", "count", "zi", "zero_inflated"),
  verbose = TRUE)

#lets make a tp figure for herring;find median value first = 146



tpherring<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= herring_only2)

predict(tpherring, herring_only2, allow.new.levels=TRUE)
ndtp<-herring_only2[1,]
ndtp$year<-"new"
ndtp_pop<-data.frame(year=herring_only2$year, site=NA, sl = 146) #median = sl
tppredict<-predict(tpherring, newdata=ndtp_pop, se.fit=TRUE)
as.data.frame(tppredict)
tppredict$year<-herring_only2$year
tppredict<-as.data.frame(tppredict)

tpherringfig<- ggplot() + geom_point(data = herring_only, aes(x = year, y = tp), color="#cd3341") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()



#Analyses for only rock
rock_only<-allspp%>%
  filter(hostsp == "rock")

normalityglu<-glm(glu~year, data = rock_only)
resrock<-simulateResiduals(fittedModel = normalityglu, n = 250)
resrock$scaledResiduals
value<-testUniformity(resrock) #normal
shapiro.test(rock_only$glu)
glurock<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= rock_only)
summary(glurock)

intglurock<-glmmTMB(glu~ scale(year)*scale(sl) +(1|site) , family = "gaussian", data= rock_only)
summary(intglurock)

glurocknoscale<-glmmTMB(glu~ year + sl +(1|site) , family = "gaussian", data= rock_only)
summary(glurocknoscale)


normalityphe<-glm(phe~year, data = rock_only)
resrock<-simulateResiduals(fittedModel = normalityphe, n = 250)
resrock$scaledResiduals
plot(resrock)
testUniformity(resrock) #normal
shapiro.test(rock_only$phe)
pherock<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= rock_only)
summary(pherock)


library(glmmTMB)
normalitytp<-glm(tp~year, data = rock_only)
resrock<-simulateResiduals(fittedModel = normalitytp, n = 250)
resrock$scaledResiduals
plot(resrock)
testUniformity(resrock) 
tprock<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= rock_only)
summary(tprock)
#plot residuals for year to see if there are non-linear trends; look at GAMS
diffrock<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= rock_only)
summary(diffrock)

#supp interaction analysis
tprockint<-glmmTMB(tp~ scale(year) + scale(sl) +scale(year)*scale(sl) +(1|site), family = "gaussian", data= rock_only)
summary(tprockint)

#make the tpfigure for publication
predict(tprock, rock_only, allow.new.levels=TRUE)
ndtp<-rock_only[1,]
ndtp$year<-"new"
ndtp_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) #median = sl
tppredict<-predict(tprock, newdata=ndtp_pop, se.fit=TRUE)
as.data.frame(tppredict)
tppredict$year<-rock_only$year
tppredict<-as.data.frame(tppredict)

tplot1rock<- ggplot() + geom_line(data =tppredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tppredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#0077b6", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = tp), color="#0077b6") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic() + theme(axis.title.x = element_blank()) +ggtitle("p < 0.001") + theme(plot.title = element_text(vjust = -5, hjust = 0.05))

#make figure for glu
predict(glurock, rock_only, allow.new.levels=TRUE)
ndglu<-rock_only[1,]
ndglu$year<-"new"
ndglu_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) 
glupredict<-predict(glurock, newdata=ndglu_pop, se.fit=TRUE)
as.data.frame(glupredict)
glupredict$year<-rock_only$year
glupredict<-as.data.frame(glupredict)

gluplot1rock<- ggplot() + geom_line(data =glupredict, aes(x = year, y = fit)) +
  geom_ribbon(data = glupredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#0077b6", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = glu), color="#0077b6") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic() + theme(axis.title.x = element_blank()) +ggtitle("p < 0.001") + theme(plot.title = element_text(vjust = -5, hjust = 0.05))

#make figure for phe
predict(pherock, rock_only, allow.new.levels=TRUE)
ndphe<-rock_only[1,]
ndphe$year<-"new"
ndphe_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) 
phepredict<-predict(pherock, newdata=ndphe_pop, se.fit=TRUE)
as.data.frame(phepredict)
phepredict$year<-rock_only$year
phepredict<-as.data.frame(phepredict)

pheplot1rock<- ggplot() + geom_point(data = rock_only, aes(x = year, y = phe), color="#0077b6") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic() +  theme(axis.title.x = element_blank()) +ggtitle("p = 0.405") + theme(plot.title = element_text(vjust = -5, hjust = 0.05))

library(patchwork)
rockfig<-tplot1rock| gluplot1rock / pheplot1rock 
rockfig1<- rockfig +  plot_annotation(title = "                               Copper Rockfish") 

ggsave("fig3.jpg")


library(ggeasy)
tpfig<-ggarrange(tphakefig + ggtitle("Pacific Hake \n p = 0.416") + easy_center_title()
 +theme(axis.title = element_blank() ),
                 tppollockfig +ggtitle("Walleye Pollock \n p = 0.975")+ easy_center_title()+
theme(axis.title = element_blank()),
                 tpsolefig + ggtitle("English Sole \n p = 0.443")+ easy_center_title()
+theme(axis.title = element_blank()),
                 tpherringfig+ggtitle("Pacific Herring \n p = 0.842") + easy_center_title()
+theme(axis.title = element_blank()), 
                 nrow= 2, ncol  = 2)
library(ggpubr)
tpfigcomplete<- annotate_figure(tpfig, top = NULL, bottom = "Year collected", left = "Trophic position")
ggsave("figure2.jpg")                               


glufig<-ggarrange(gluhakefig + ggtitle("Pacific Hake \n p = 0.993") + easy_center_title()
                 +theme(axis.title = element_blank() ),
                 glupollockfig +ggtitle("Walleye Pollock \n p = 0.435")+ easy_center_title()+
                   theme(axis.title = element_blank()),
                glusolefig + ggtitle("English Sole \n p = 0.998")+ easy_center_title()
                 +theme(axis.title = element_blank()),
                gluherringfig+ggtitle("Pacific Herring \n p = 0.285") + easy_center_title()
                 +theme(axis.title = element_blank()), 
                 nrow= 2, ncol  = 2) + xlab("Year collected") + ylab("Glutamic Acid (in ‰)")
glufigsupp<-annotate_figure(
  glufig,
  bottom = "Year Collected",
  left = "Glutamic Acid (in ‰)")
ggsave("glufigsupp.jpg")

phefig<-ggarrange(phehakefig + ggtitle("Pacific Hake \n p = 0.361") + easy_center_title()
                  +theme(axis.title = element_blank() ),
                  phepollockfig +ggtitle("Walleye Pollock \n p = 0.320")+ easy_center_title()+
                    theme(axis.title = element_blank()),
                  phesolefig + ggtitle("English Sole \n p = 0.125")+ easy_center_title()
                  +theme(axis.title = element_blank()),
                  pheherringfig+ggtitle("Pacific Herring \n p = 0.127") + easy_center_title()
                  +theme(axis.title = element_blank()), 
                  nrow= 2, ncol  = 2) + xlab("Year collected") + ylab("Phenylanaine(in ‰)")
phefigsupp<-annotate_figure(
  phefig,
  bottom = "Year Collected",
  left = "Phenylalanine (in ‰)")
ggsave("phefigsupp.jpg")




#make data tables
library(kableExtra)
glmmTable <- function(glmmTMBobject, AAtype, filename){
  rawtab <- summary(glmmTMBobject)$coefficients$cond
  rownames(rawtab) <- c("Intercept","scale(year)", "scale(sl)")
  # Making the table
  rawtab[c(2:3),c(1:4)] %>% 
    kbl(col.names = c("Estimate", "Standard Error", "z-value","p-value"), 
        digits=c(3,3,3,3), 
        align = "c") %>% 
    kable_styling(full_width = FALSE, 
                  html_font = "Times New Roman") %>% 
    save_kable(file = filename, self_contained = T) 
}

glmmTable(glusole,"Glu_Sole", "glu_sole.html")
glmmTable(phesole,"Phe_Sole", "phe_sole.html")
glmmTable(tpsole,"TP_Sole", "tp_sole.html")
glmmTable(gluhake,"Glu_Hake", "glu_hake.html")
glmmTable(phehake,"Phe_Hake", "phe_hake.html")
glmmTable(tphake,"TP_Hake", "tp_hake.html")
glmmTable(glupollock,"Glu_pollock", "glu_pollock.html")
glmmTable(phepollock,"Phe_pollock", "phe_pollock.html")
glmmTable(tppollock,"TP_pollock", "tp_pollock.html")
glmmTable(gluherring,"Glu_herring", "glu_herring.html")
glmmTable(pheherring,"Phe_herring", "phe_herring.html")
glmmTable(tpherring,"TP_herring", "tp_herring.html")
glmmTable(glurock,"Glu_rock", "glu_rock.html")
glmmTable(pherock,"Phe_rock", "phe_rock.html")
glmmTable(tprock,"TP_rock", "tp_rock.html")

save(glusole, phesole, tpsole, gluhake, phehake, tphake, glupollock, phepollock, tppollock, gluherring, pheherring, tpherring, glurock, pherock, tprock,  file="models.Rdata")


library(mgcv) #a package to run GAMS
library(mgcViz) #a package to plot GAMS
allspp$site<-as.factor(allspp$site)
rock_only$site<-as.factor(rock_only$site)


phemodgam<-gam(phe~ s(year,k=4)+ sl + s(site,bs="re"), data=rock_only, method="REML")
summary(phemodgam)
glumodgam<-gam(glu~ s(year,k=5)+ sl  + s(site,bs="re"), data=rock_only, method="REML")
summary(glumodgam)

phemodgamsole<-gam(phe~ s(year,k=5)+ sl + s(site,bs="re"), data=sole_only, method="REML")
summary(phemodgamsole)
glumodgamsole<-gam(glu~ s(year,k=5)+ sl  + s(site,bs="re"), data=sole_only, method="REML")
summary(glumodgamsole)

phemodgamhake<-gam(phe~ s(year,k=5)+ sl + s(site,bs="re"), data=hake_only, method="REML")
summary(phemodgamhake)
glumodgamhake<-gam(glu~ s(year,k=5)+ sl  + s(site,bs="re"), data=hake_only, method="REML")
summary(glumodgamhake)

phemodgamherring<-gam(phe~ s(year,k=5)+ sl + s(site,bs="re"), data=herring_only, method="REML")
summary(phemodgamherring)
glumodgamherring<-gam(glu~ s(year,k=5)+ sl  + s(site,bs="re"), data=herring_only, method="REML")
summary(glumodgamherring)


phemodgampollock<-gam(phe~ s(year,k=5)+ sl + s(site,bs="re"), data=pollock_only, method="REML")
summary(phemodgampollock)
glumodgampollock<-gam(glu~ s(year,k=5)+ sl  + s(site,bs="re"), data=pollock_only, method="REML")
summary(glumodgampollock)



fishstatsdecade<-allspp%>%
  group_by(hostsp,decade)%>%
  summarize(
    meantp = mean(tp, na.rm = TRUE),
    tpsd = sd(tp, na.rm = TRUE),
    avgglu = mean(glu, na.rm = TRUE), 
    sdglu = sd (glu, na.rm = TRUE),
    avgphe = mean(phe, na.rm = TRUE),
    sdphe = sd(phe, na.rm = TRUE))
write.csv(fishstats, "descriptiveresultstablebydecade.csv")

fishstatsoverall<-allspp%>%
  group_by(hostsp)%>%
  summarize(
    meantp = mean(tp, na.rm = TRUE),
    tpsd = sd(tp, na.rm = TRUE),
    avgglu = mean(glu, na.rm = TRUE), 
    sdglu = sd (glu, na.rm = TRUE),
    avgphe = mean(phe, na.rm = TRUE),
    sdphe = sd(phe, na.rm = TRUE))
write.csv(fishstats, "descriptiveresultstableoverall.csv")


solen<-sole_only%>%
  group_by(decade)%>%
  summarize(
    fishcount = n_distinct(hostid, na.rm = TRUE))

rockn<-rock_only%>%
  group_by(decade)%>%
  summarize(
    fishcount = n_distinct(hostid, na.rm = TRUE))

pollockn<-pollock_only%>%
  group_by(decade)%>%
  summarize(
    fishcount = n_distinct(hostid, na.rm = TRUE))


haken<-hake_only%>%
  group_by(decade)%>%
  summarize(
    fishcount = n_distinct(hostid, na.rm = TRUE))

herringn<-herring_only%>%
  group_by(decade)%>%
  drop_na(decade)%>%
  summarize(
    fishcount = n_distinct(hostid, na.rm = TRUE))




write.csv(fishstats2, "decadesamplesize.csv")






library(ggplot2)
samplesizefig<- ggplot(fishstats2, aes( x= decade, y = fishcount, fill = hostsp, color=hostsp)) + geom_bar()

#rockfish glu time series
library(strucchange)
ocus.glu<-efp(glurock, type = "Rec-CUSUM", data = rock_only)
bound.ocus.glu<-boundary(ocus.glu, alpha = 0.1)
plot(ocus.glu)
sctest(ocus.glu)
