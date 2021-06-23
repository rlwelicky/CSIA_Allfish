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
#now lets combine the sp. files into one datatable
speciescombined<-bind_rows(hake,herring, pollock, sole,rock,reruns)
#I need to make a column of final values for each AA for indiv analysis. # I included this values_fn statement because if there are duplicate samples, we don't know which best, so we will take the average of the duplicate.
allspp<-speciescombined%>%
  filter(!is.na(c(finalvalue)))%>%
  pivot_wider(id_cols = c(hostid, hostsp, sl, tl, lat, long, year, decade, lag0_temp, lag1_temp, lag2_temp, lag3_temp, lag12_temp), names_from = AA, values_from = finalvalue)
duplicated(allspp$hostid)#there are a few duplicated but its ok, bc one of the duplicate doensn't have glu and phe so it won't be used anyways
#I need to make o)ne variable called site using lat and long
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
            mediansl = median(sl, na.rm = TRUE),
            fishcount = n_distinct(hostid, na.rm = TRUE),
            averagetp = mean (tp, na.rm = TRUE))

library(DHARMa)
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
library(spdep)
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
library(lmtest)
library(car)
allspp$timegroup<-allspp$year
temporal.glu<-glm(glu~year, data= allspp)
temporaltest.glu<-dwtest(temporal.glu, order.by = NULL, alternative = "two.sided", exact = FALSE, tol = 1e-10) 
temporaltest.glu #yay, no temporal corr
dwtest(temporal.glu)

#Analyses with all species in one model, no lag time of temp
#does time influence glu?
library(glmmTMB)
library(multcomp)
gluint<-glmmTMB(glu~ year  + hostsp +  (1|site)  , family = "gaussian", data= allspp)
summary(gluint) #AIC 781.5; 749
glu<-glmmTMB(glu~ year + sl + (1|hostsp)  + (1|site) , family = "gaussian", data= allspp)
summary(glu) #AIC 748

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
glusole<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site)  , family = "gaussian", data= sole_only)
summary(glusole)
phesole<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= sole_only)
summary(phesole)
tpsole<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= sole_only)
summary(tpsole)



#Analyses for only hake
hake_only<-allspp%>%
  filter(hostsp == "hake")

gluhake<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= hake_only)
summary(gluhake)
phehake<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(phehake)
tphake<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= hake_only)
summary(tphake)


#Analyses for only pollock
pollock_only<-allspp%>%
  filter(hostsp == "pollock")

glupollock<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(glupollock)
phepollock<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(phepollock)
tppollock<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= pollock_only)
summary(tppollock)

#Analyses for only herring
herring_only<-allspp%>%
  filter(hostsp == "herring")

gluherring<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(gluherring)
herringresid<-resid(gluherring)
plot(herringresid, herring_only$year)

pheherring<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(pheherring)
tpherring<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= herring_only)
summary(tpherring)

#Analyses for only rock
rock_only<-allspp%>%
  filter(hostsp == "rock")

glurock<-glmmTMB(glu~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= rock_only)
summary(glurock)
pherock<-glmmTMB(phe~ scale(year) + scale(sl) +(1|site) , family = "gaussian", data= rock_only)
summary(pherock)
diffrock<-glmmTMB(diff~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= rock_only)
summary(diffrock)
tprock<-glmmTMB(tp~ scale(year) + scale(sl) +(1|site), family = "gaussian", data= rock_only)
summary(tprock)
#plot residuals for year to see if there are non-linear trends; look at GAMS

#make the tpfigure for publication

predict(tprock, rock_only, allow.new.levels=TRUE)
ndtp<-rock_only[1,]
ndtp$year<-"new"
ndtp_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) #random effects are set to NA, other effects I chose median of climate and mean for sl
tppredict<-predict(tprock, newdata=ndtp_pop, se.fit=TRUE)
as.data.frame(tppredict)
tppredict$year<-rock_only$year
tppredict<-as.data.frame(tppredict)

tplot1<- ggplot() + geom_line(data =tppredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tppredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = tp), color="#00529c") +  xlab("Year collected") + ylab("tp") +theme_classic()

#make figure for glu
predict(glurock, rock_only, allow.new.levels=TRUE)
ndglu<-rock_only[1,]
ndglu$year<-"new"
ndglu_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) #random effects are set to NA, other effects I chose median of climate and mean for sl
glupredict<-predict(glurock, newdata=ndglu_pop, se.fit=TRUE)
as.data.frame(glupredict)
glupredict$year<-rock_only$year
glupredict<-as.data.frame(glupredict)

gluplot1<- ggplot() + geom_line(data =glupredict, aes(x = year, y = fit)) +
  geom_ribbon(data = glupredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = glu), color="#00529c") +  xlab("Year collected") + ylab("glu") +theme_classic()

#make figure for phe
predict(pherock, rock_only, allow.new.levels=TRUE)
ndphe<-rock_only[1,]
ndphe$year<-"new"
ndphe_pop<-data.frame(year=rock_only$year, site=NA, sl = 180) #random effects are set to NA, other effects I chose median of climate and mean for sl
phepredict<-predict(pherock, newdata=ndphe_pop, se.fit=TRUE)
as.data.frame(phepredict)
phepredict$year<-rock_only$year
phepredict<-as.data.frame(phepredict)

pheplot1<- ggplot() + geom_line(data =phepredict, aes(x = year, y = fit)) +
  geom_ribbon(data = phepredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data = rock_only, aes(x = year, y = phe), color="#00529c") +  xlab("Year collected") + ylab("phe") +theme_classic()

