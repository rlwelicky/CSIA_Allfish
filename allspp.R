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


#lets check out how bad our distribution are...
normalityglu<-glm(glu~year, data = allspp)
resglu<-simulateResiduals(fittedModel = normalityglu, n = 250)
plot(resglu)
resglu$scaledResiduals
testUniformity(resglu)
#whoa, yea for normal glu! onto phe...
normalityphe<-glm(phe~year, data = allspp)
resphe<-simulateResiduals(fittedModel = normalityphe, n = 250)
resphe$scaledResiduals
plot(resphe)
testUniformity(resphe) 
#wow, lucky me--normal!...something will be off soon, haha. now tp...
normalitytp<-glm(tp~year, data = allspp)
restp<-simulateResiduals(fittedModel = normalitytp, n = 250)
restp$scaledResiduals
plot(restp)
testUniformity(restp)  #well, normal : normal, is normal, expected. good.

#curious about normality of the diff
normalitydiff<-glm(diff~year, data = allspp)
resdiff<-simulateResiduals(fittedModel = normalitydiff, n = 250)
resdiff$scaledResiduals
plot(resdiff)
testUniformity(resdiff)

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

#Analyses with all species in one model, no lag time of temp
#does time influence glu?
library(lmerTest)
gluint<-glmmTMB(glu~ year  + hostsp +  (1|site)  , family = "gaussian", data= allspp)
summary(gluint) #AIC 781.5; 749
glu<-glmmTMB(glu~ year + sl + (1|hostsp)  + (1|site) , family = "gaussian", data= allspp)
summary(glu) #AIC 748
glulmer<-lmer(glu~ year  + hostsp +  (1|site), data= allspp)
summary(glulmer)
#does time influence phe?
pheint<-glmmTMB(phe~ year + sl + hostsp  + (1|site) + lag3_temp, family = "gaussian", data= allspp)
summary(pheint) #AIC 811;779
phe<-glmmTMB(phe~ year  + (1|hostid)  + hostsp , family = "gaussian", data= allspp)
summary(phe) #806;776


#does time influence diff?
diffint<-glmmTMB(diff~ year  + sl + (1|site) , family = "gaussian", data= allspp)
summary(diffint) #AIC 857
diff<-glmmTMB(diff~ year + sl + hostsp  + (1|site), family = "gaussian", data= allspp)
summary(diff) #853

#does time influence tp?
tpint<-glmmTMB(tp~ year + hostsp + sl   +(1|site) , family = "gaussian", data= allspp)
summary(tpint) #AIC  18
tp<-glmmTMB(diff~ year + sl + hostsp  + (1|site), family = "gaussian", data= allspp)
summary(tp) 

ab<-glht(tpint, linfct = mcp(hostsp = "Tukey"))
summary(ab)


#Analyses for only Sole

sole_only<-allspp%>%
  filter(hostsp == "sole")
library(predictmeans)
normalityglu<-glm(glu~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalityglu, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole) #glu is normally distributed
glusole<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site)  , family = "gaussian", data= sole_only)
summary(glusole)

glusolelmer<-lmer(glu~ year + sl +(1|site), data= sole_only)
summary(glusolelmer)
library(visreg)
a<-visreg(glusolelmer, "year")
view (a)
normalityphe<-glm(phe~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalityphe, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole)#phe is normally distrbuted

phesole<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= sole_only)
summary(phesole)


normalitytp<-glm(tp~year, data = sole_only)
ressole<-simulateResiduals(fittedModel = normalitytp, n = 250)
ressole$scaledResiduals
plot(ressole)
testUniformity(ressole) #tp is normally distributed
tpsole<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= sole_only)
summary(tpsole)

diffsole<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= sole_only)
summary(diffsole)


#make tp figure for sole
#what's the median sl needed?
sole_only%>%
  summarise(
    mediansl = median(sl)
  )

predict(tpsole, sole_only, allow.new.levels=TRUE)
ndtpsole<-sole_only[1,]
ndtpsole$year<-"new"
ndtpsole_pop<-data.frame(year=sole_only$year, site=NA, sl = 126) #random effects are set to NA, other effects I chose median of climate and mean for sl
tpsolepredict<-predict(tpsole, newdata=ndtpsole_pop, se.fit=TRUE)
as.data.frame(tpsolepredict)
tpsolepredict$year<-sole_only$year
tpsolepredict<-as.data.frame(tpsolepredict)

tpsolefig<- ggplot() + geom_line(data =tpsolepredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tpsolepredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="darkorchid4", alpha=0.3) + geom_point(data = sole_only, aes(x = year, y = tp), color="darkorchid4") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()



#Analyses for only hake
hake_only<-allspp%>%
  filter(hostsp == "hake")

normalityglu<-glm(glu~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalityglu, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal :)
gluhake<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= hake_only)
summary(gluhake)


normalityphe<-glm(phe~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalityphe, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal :) 
phehake<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(phehake)


normalitytp<-glm(tp~year, data = hake_only)
reshake<-simulateResiduals(fittedModel = normalitytp, n = 250)
reshake$scaledResiduals
plot(reshake)
testUniformity(reshake) #normal
tphake<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(tphake)

diffhake<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= hake_only)
summary(diffhake)


#lets make a tp figure for hake;find median value first = 162

hake_only%>%
  summarise(
    mediansl = median(sl)
  )

predict(tphake, hake_only, allow.new.levels=TRUE)
ndtphake<-hake_only[1,]
ndtphake$year<-"new"
ndtphake_pop<-data.frame(year=hake_only$year, site=NA, sl = 162) #random effects are set to NA, other effects I chose median of climate and mean for sl
tphakepredict<-predict(tphake, newdata=ndtphake_pop, se.fit=TRUE)
as.data.frame(tphakepredict)
tphakepredict$year<-hake_only$year
tphakepredict<-as.data.frame(tphakepredict)

tphakefig<- ggplot() + geom_line(data =tphakepredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tphakepredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="darkorchid4", alpha=0.3) + geom_point(data = hake_only, aes(x = year, y = tp), color="darkorchid4") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()

#Analyses for only pollock
pollock_only<-allspp%>%
  filter(hostsp == "pollock")

normalityglu<-glm(glu~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalityglu, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock) #normal
glupollock<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(glupollock)


normalityphe<-glm(phe~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalityphe, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock)#normal
phepollock<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(phepollock)


normalitytp<-glm(tp~year, data = pollock_only)
respollock<-simulateResiduals(fittedModel = normalitytp, n = 250)
respollock$scaledResiduals
plot(respollock)
testUniformity(respollock)
tppollock<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(tppollock)

diffpollock<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= pollock_only)
summary(diffpollock)

#lets make a tp figure for pollock; find median value first = 108
pollock_only %>%
  summarise(
  mediansl = median(sl, na.rm = TRUE)
  )

predict(tppollock, pollock_only, allow.new.levels=TRUE)
ndtppollock<-pollock_only[1,]
ndtppollock$year<-"new"
ndtppollock_pop<-data.frame(year=pollock_only$year, site=NA, sl = 108) #random effects are set to NA, other effects I chose median of climate and mean for sl
tppollockpredict<-predict(tppollock, newdata=ndtppollock_pop, se.fit=TRUE)
as.data.frame(tppollockpredict)
tppollockpredict$year<-pollock_only$year
tppollockpredict<-as.data.frame(tppollockpredict)

tppollockfig<- ggplot() + geom_line(data =tppollockpredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tppollockpredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="darkorchid4", alpha=0.3) + geom_point(data = pollock_only, aes(x = year, y = tp), color="darkorchid4") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()

#Analyses for only herring
herring_only<-allspp%>%
  filter(hostsp == "herring")


normalityglu<-glm(glu~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalityglu, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
gluherring<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(gluherring)


normalityphe<-glm(phe~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalityphe, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
pheherring<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(pheherring)


normalitytp<-glm(tp~year, data = herring_only)
resherring<-simulateResiduals(fittedModel = normalitytp, n = 250)
resherring$scaledResiduals
plot(resherring)
testUniformity(resherring) #normal
tpherring<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), na.action = na.omit, family = "gaussian", data= herring_only)
summary(tpherring)


diffherring<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(diffherring)
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

#lets make a tp figure for herring;find median value first = 142

#I can't seem to get the herring figure to work! Not "Error in eval(substitute(expr), data, enclos = parent.frame()) : 
#NAs's causing problems, so I removed them to make the figures only!
herring_only2<-herring_only %>% 
  drop_na(site, year, sl) 

tpherring<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= herring_only2)

predict(tpherring, herring_only2, allow.new.levels=TRUE)
ndtp<-herring_only2[1,]
ndtp$year<-"new"
ndtp_pop<-data.frame(year=herring_only2$year, site=NA, sl = 142) #median = sl
tppredict<-predict(tpherring, newdata=ndtp_pop, se.fit=TRUE)
as.data.frame(tppredict)
tppredict$year<-herring_only2$year
tppredict<-as.data.frame(tppredict)

tpherringfig<- ggplot() + geom_line(data =tppredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tppredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="darkorchid4", alpha=0.3) + geom_point(data = herring_only, aes(x = year, y = tp), color="darkorchid4") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic()



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
  geom_ribbon(data = tppredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="darkorchid4", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = tp), color="darkorchid4") +  xlab("Year collected") + ylab("Trophic Position") +theme_classic() + theme(axis.title.x = element_blank())

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
  geom_ribbon(data = glupredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="olivedrab", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = glu), color="olivedrab") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic() + theme(axis.title.x = element_blank())

#make figure for phe
predict(pherock, rock_only, allow.new.levels=TRUE)
ndphe<-rock_only[1,]
ndphe$year<-"new"
ndphe_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) 
phepredict<-predict(pherock, newdata=ndphe_pop, se.fit=TRUE)
as.data.frame(phepredict)
phepredict$year<-rock_only$year
phepredict<-as.data.frame(phepredict)

pheplot1rock<- ggplot() + geom_line(data =phepredict, aes(x = year, y = fit)) +
  geom_ribbon(data = phepredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="dodgerblue3", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = phe), color="dodgerblue3") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic() +  theme(axis.title.x = element_blank())

library(patchwork)
rockfig<-tplot1rock| gluplot1rock / pheplot1rock 
ggsave("fig3.jpg")


library(ggeasy)
tpfig<-ggarrange(tphakefig + ggtitle("Pacific Hake") + easy_center_title()
 +theme(axis.title = element_blank() ),
                 tppollockfig +ggtitle("Walleye Pollock")+ easy_center_title()+
theme(axis.title = element_blank()),
                 tpsolefig + ggtitle("English Sole")+ easy_center_title()
+theme(axis.title = element_blank()),
                 tpherringfig+ggtitle("Pacific Herring") + easy_center_title()
+theme(axis.title = element_blank()), 
                 nrow= 2, ncol  = 2)
library(ggpubr)
tpfigcomplete<- annotate_figure(tpfig, top = NULL, bottom = "Year collected", left = "Trophic position")
ggsave("figure2.jpg")                               

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
