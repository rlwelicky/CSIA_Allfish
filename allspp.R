#lets input each spp file
sole<-read.csv("Sole_CSIA_withmetadata.csv", header=TRUE)
hake<-read.csv('Hake_CSIA_withmetadata.csv', header=TRUE) #one fish was manually removed because it was from the coast.
library(tidyverse)
library(dplyr)
library(broom)
#now lets combine the sp. files into one datatable
speciescombined<-bind_rows(sole,hake)
#I need to make a column of final values for each AA for indiv analysis. # I included this values_fn statement because if there are duplicate samples, we don't know which best, so we will take the average of the duplicate.
allspp<-speciescombined%>%
  filter(!is.na(c(finalvalue)))%>%
  pivot_wider(id_cols = c(hostid, hostsp, sl, tl, lat, long, year, decade, lag0_temp, lag1_temp, lag2_temp, lag3_temp, lag12_temp), names_from = AA, values_from = finalvalue)
#I need to make o)ne variable called site using lat and long
allspp$site=paste(allspp$lat,allspp$long)
#cI need to reate a column of trophic position calculations
allspp$tp=(((allspp$glu-allspp$phe-3.4)/7.6)+1)
allspp$diff=(allspp$glu-allspp$phe)

allspp$gluadjust= (allspp$glu/allspp$sl)
allspp$diffadj=(allspp$gluadjust-allspp$phe)
allspp$tpadjust=(((allspp$gluadjust-allspp$phe-3.4)/7.6)+1)
summary(allspp$tp)
summary(allspp$diff)
library(DHARMa)
#lets check out how bad our distribution are...
normalityglu<-glm(glu~year, data = allspp)
resglu<-simulateResiduals(fittedModel = normalityglu, n = 250)
resglu$scaledResiduals
plot(resglu)
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
spatialcheck<-allspp%>%
  filter(!is.na(c(lat)))%>%
  filter(!is.na(c(long)))
spatialcheck$latjitt<-jitter(spatialcheck$lat, factor=0.1, amount=NULL)
spatialcheck$longjitt<-jitter(spatialcheck$long, factor=0.1, amount=0)
spatial.allspp<-glm(glu~year, data=spatialcheck)
simspatial.allspp<-simulateResiduals(fittedModel = spatial.allspp)
spatialtest.allspp<-testSpatialAutocorrelation(simulationOutput = simspatial.allspp,  x = spatialcheck$longjitt, y = spatialcheck$latjitt)
spatialtest.allspp #Yay, no spatial corr

#temporal checks for glutamic acid, his would be the same as phe or tp, since all these values are from 1 fish from the same site
library(lmtest)
timegroup<-allspp$year
temporal.glu<-glm(glu~ year, data= allspp)
temporaltest.glu<-dwtest(temporal.glu, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) 
temporaltest.glu #yay, no temporal corr


#Analyses with all species in one model
#does time influence glu?
library(glmmTMB)
glu0<-glmmTMB(glu~ year + sl + hostsp + lag0_temp + (1|site), family = "gaussian", data= allspp)
summary(glu0)
glu1<-glmmTMB(glu~ year + sl + hostsp + lag1_temp + (1|site) , family = "gaussian", data= allspp)
summary(glu1)
glu2<-glmmTMB(glu~ year + sl + hostsp + lag2_temp + (1|site), family = "gaussian", data= allspp)
summary(glu2)
glu3<-glmmTMB(glu~ year + sl + hostsp + lag3_temp + (1|site) , family = "gaussian", data= allspp)
summary(glu3)
glu12<-glmmTMB(glu~ year + sl + hostsp + lag12_temp + (1|site) , family = "gaussian", data= allspp)
summary(glu12)

#does time influence phe?
phe0<-glmmTMB(phe~ year + sl + hostsp + lag0_temp + (1|site),  family = "gaussian", data= allspp)
summary(phe0)
phe1<-glmmTMB(phe~ year + sl + hostsp + lag1_temp + (1|site) , family = "gaussian", data= allspp)
summary(phe1)
phe2<-glmmTMB(phe~ year + sl + hostsp + lag2_temp + (1|site) , family = "gaussian", data= allspp)
summary(phe2)
phe3<-glmmTMB(phe~ year + sl + hostsp + lag3_temp + (1|site) , family = "gaussian", data= allspp)
summary(phe3)
phe12<-glmmTMB(phe~ year + sl + hostsp + lag12_temp + (1|site) , family = "gaussian", data= allspp)
summary(phe12)
#does time influence tp?
tp0<-glmmTMB(tp~ year + sl + hostsp + lag0_temp + (1|site) , family = "gaussian", data= allspp)
summary(tp0)
tp1<-glmmTMB(tp~ year + sl + hostsp + lag1_temp + (1|site) , family = "gaussian", data= allspp)
summary(tp1)
tp2<-glmmTMB(tp~ year + sl + hostsp + lag2_temp + (1|site)  , family = "gaussian", data= allspp)
summary(tp2)
tp3<-glmmTMB(tp~ year + sl + hostsp + lag3_temp + (1|site)  , family = "gaussian", data= allspp)
summary(tp3)
tp12<-glmmTMB(tp~ year + sl + hostsp + lag12_temp + (1|site) , family = "gaussian", data= allspp)
summary(tp12)

#does time influence diff?
diff0<-glmmTMB(diff~ year + sl + hostsp + lag0_temp + (1|site) , family = "gaussian", data= allspp)
summary(diff0)
diff1<-glmmTMB(diff~ year + sl + hostsp + lag1_temp + (1|site) , family = "gaussian", data= allspp)
summary(diff1)
diff2<-glmmTMB(diff~ year + sl + hostsp + lag2_temp + (1|site)  , family = "gaussian", data= allspp)
summary(diff2)
diff3<-glmmTMB(diff~ year + sl + hostsp + lag3_temp + (1|site)  , family = "gaussian", data= allspp)
summary(diff3)
diff12<-glmmTMB(diff~ year + sl + hostsp + lag12_temp + (1|site) , family = "gaussian", data= allspp)
summary(diff12)
  

#Analyses for only Sole
sole_only<-allspp%>%
  filter(hostsp == "sole")
glu2sole<-glmmTMB(glu~ year + sl + lag2_temp + (1|site) , family = "gaussian", data= sole_only)
summary(glu2sole)
glu12sole<-glmmTMB(glu~ year + sl + lag12_temp + (1|site) , family = "gaussian", data= sole_only)
summary(glu12sole)
phe2sole<-glmmTMB(phe~ year + sl + lag2_temp + (1|site) , family = "gaussian", data= sole_only)
summary(phe2sole)
phe12sole<-glmmTMB(phe~ year + sl + lag12_temp + (1|site), family = "gaussian", data= sole_only)
summary(phe12sole)
tp2sole<-glmmTMB(tp~ year + sl + lag2_temp + (1|site) , family = "gaussian", data= sole_only)
summary(tp2sole)
tp12sole<-glmmTMB(tp~ year + sl + lag12_temp + (1|site) , family = "gaussian", data= sole_only)
summary(tp12sole)

diff1<-glmmTMB(diffadj~ year +  (1|site) , family = "gaussian", data= sole_only)
summary(diff1)

#Analyses for only hake
hake_only<-allspp%>%
  filter(hostsp == "Merluccius_productus")
hist(hake_only$phe)
glu2hake<-glmmTMB(glu~ scale(year) + sl + lag2_temp + (1|site) , family = "gaussian", data= hake_only)
summary(glu2hake)
glu12hake<-glmmTMB(glu~ scale(year) + sl + lag12_temp + (1|site) , family = "gaussian", data= hake_only)
summary(glu12hake)
phe2hake<-glmmTMB(phe~ scale(year) + sl + lag2_temp + (1|site) , family = "gaussian", data= hake_only)
summary(phe2hake)
phe2hake<-glmmTMB(phe~ scale(year) + sl + lag2_temp + (1|site) , family = "Gamma", data= hake_only)
summary(phe2hake)
phe12hake<-glmmTMB(phe~ scale(year) + sl + lag12_temp + (1|site)  , family = "gaussian", data= hake_only)
summary(phe12hake)
tp2hake<-glmmTMB(tp~ scale(year) + sl + lag2_temp + (1|site), family = "gaussian", data= hake_only)
summary(tp2hake)
tp12hake<-glmmTMB(tp~ scale(year) + sl + lag12_temp + (1|site)  , family = "gaussian", data= hake_only)
summary(tp12hake)

phemod<-glmmTMB(phe~ year + (1|site) , family = "gaussian", data= sole_only)
summary(phemod)
glumod<-glmmTMB(glu/sl~ year + (1|site) , family = "gaussian", data= sole_only)
summary(glumod)

phemodhake<-glmmTMB(phe~ scale(year) + (1|site) , family = "gaussian", data= hake_only)
summary(phemodhake)
glumodhake<-glmmTMB(glu/sl~ year + (1|site) , family = "gaussian", data= hake_only)
summary(glumodhake)


pheall<-glmmTMB(log(phe)~ scale(year) + (1|hostsp) + (1|site) , family = "gaussian", data= allspp)
summary(pheall)
library(lmerTest)
tpall<-lmer(phe~ year +(1|site), data=hake_only)
summary(tpall)


