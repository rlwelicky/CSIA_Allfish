csiadata<-read.csv("CompiledData.csv", header=TRUE)
indivparasitedata<-read.csv("IndivParasiteforCSIA.csv", header=TRUE)
climatedata<-read.csv("climatedata.csv", header=TRUE)
#groupparasitedata<-read.csv("GroupParasiteforCSIA.csv", header=TRUE)
library(tidyverse)
library(dplyr)
library(lsmeans)
library(lmtest)
library(DHARMa)
library(lme4)
library(geoR)
library(spacetime)
library(gstat)
library(sf)
library(sp)
library(nlme)
library(MASS)
library(ggplot2)
library(patchwork)
library(gdm)
library(vegan)
library(glmmTMB)
library(psych)
library(lme4)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(performance)
library(ggplot2)
library(kableExtra)
as.character(csiadata$hostid)
as.factor(indivparasitedata$decade)
#compiled<-left_join(csiadata,indivparasitedata, by = "hostid", copy = FALSE, suffix = c(".x", ".y"))
#datawide<-dcast(compiled, Fish.ID + hostid  +Year + Trophic.position + SEAA + diff + tl + sl + long + lat ~ AA, value.var = "MeanAA")

csiadatawide<-as_tibble(csiadata)%>%
  pivot_wider(id_cols = hostid, names_from = AA, values_from = MeanAA)
compiled1<-left_join(csiadatawide,indivparasitedata, by = "hostid", copy = FALSE, suffix = c(".x", ".y"))
compiled<-left_join(compiled1,climatedata, by = "hostid", copy = FALSE, suffix = c(".x", ".y"))

#compiled <- compiled %>% rename( MEI = MEIhttps...psl.noaa.gov.enso.mei.ext.table.ext.html)

compiled$sites=paste(compiled$lat,compiled$long)
compiled$logsl= log(compiled$sl)
compiled$logtl= log(compiled$tl)
as.factor(compiled$sites)
as.factor(compiled$PDO_nolag)
as.factor(compiled$PDO_1lag)
as.factor(compiled$PDO_2lag)
as.factor(compiled$PDO_3lag)
as.factor(compiled$MEI)
as.factor(compiled$NPGO)

yearsc<-scale(compiled$year, center = TRUE, scale = TRUE)
PDOsc<-scale(compiled$PDO, center = TRUE, scale = TRUE)
NPGOsc<-scale(compiled$NPGO, center = TRUE, scale = TRUE)
MEIsc<-scale(compiled$MEI, center = TRUE, scale = TRUE)
print(yearsc)
#create table function to store my results
glmmTable <- function(glmmTMBobject, PDOtype, filename){
  rawtab <- summary(glmmTMBobject)$coefficients$cond
  rownames(rawtab) <- c("Intercept","NPGO",
                        PDOtype,"Year", "SL")
  # Making the table
  rawtab[c(2:5),c(1:4)] %>% 
    kbl(col.names = c("Estimate", "Standard Error", "z-value","p-value"), 
        digits=c(1,1,2,3), 
        align = "c") %>% 
    kable_styling(full_width = FALSE, 
                  html_font = "Times New Roman") %>% 
    save_kable(file = filename, self_contained = T) 
}

#create a column of trophic position calculations
compiled$tp=(((compiled$GlutamicAcid-compiled$Phenylalanine-3.4)/7.6)+1)
summary(compiled$tp)
compiled$diff=((compiled$GlutamicAcid-compiled$Phenylalanine))
summary(compiled$diff)

histogram(compiled$GlutamicAcid)
shapiro.test(compiled$GlutamicAcid)
normality<-glm(GlutamicAcid~year, data = compiled)
res<-simulateResiduals(fittedModel = normality, n = 250)
res$scaledResiduals
plot(res)
testUniformity(res)

#spatial checks for glutamic acid
compiled$latjitt<-jitter(compiled$lat, factor=0.1, amount=NULL)
compiled$longjitt<-jitter(compiled$long, factor=0.1, amount=0)
spatial.compiled<-glm(GlutamicAcid~year, data=compiled)
simspatial.compiled<-simulateResiduals(fittedModel = spatial.compiled)
spatialtest.compiled<-testSpatialAutocorrelation(simulationOutput = simspatial.compiled,  x = compiled$longjitt, y = compiled$latjitt)
spatialtest.compiled

#temporal checks for glutamic acid
timegroup<-compiled$year
temporal.glu<-glm(GlutamicAcid~ year, data= compiled)
temporaltest.glu<-dwtest(temporal.glu, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) 
temporaltest.glu

#which envi indices should be included?
cor.test(compiled$NPGO, compiled$MEI) # correlated, use one
cor.test(compiled$NPGO, compiled$PDO_3lag) #not correlated, use both; same for all PDOs
cor.test(compiled$MEI, compiled$PDO_3lag) #correlated use, one; asame for all PDOs

#test for collinearity
library(olsrr)
a<-lm(GlutamicAcid~ NPGO+  PDO_2lag  + sl, family = "gaussian", data= compiled)
ols_coll_diag(a)
#does time influence glu?

glu1<-glmmTMB(GlutamicAcid~ NPGO+  PDO_1lag + scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(glu1)
glu2<-glmmTMB(GlutamicAcid~ NPGO+  PDO_2lag + scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(glu2)
glu0<-glmmTMB(GlutamicAcid~ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(glu0)
glu3<-glmmTMB(GlutamicAcid~ NPGO+  PDO_3lag + scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(glu3)
glmmTable(glu1,"PDO", "glu1.html")
glmmTable(glu2,"PDO", "glu2.html")
glmmTable(glu3,"PDO", "glu3.html")
library(lme4)
predict(glu2, compiled, allow.new.levels=TRUE)
nd_pop<-data.frame(year=compiled$year, sites=NA, NPGO = -0.011230659, PDO_2lag = -0.291501502, sl = 127.33)
predictglu<-predict(glu2, newdata=nd_pop, se.fit=TRUE)
as.data.frame(predictglu)
predictglu$year<-compiled$year
predictglu<-as.data.frame(predictglu)

#make the glu figure for publication
library(ggplot2)

gluplot1<- ggplot() + geom_line(data =predictglu, aes(x = year, y = fit)) +
  geom_ribbon(data = predictglu, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data = compiled, aes(x = year, y = GlutamicAcid), color="#00529c") +  xlab("Year collected") + ylab("Glutamic acid (in ‰)") +theme_classic()
ggsave("gluplot1.jpeg")


library(dplyr)
#decade glu infor
decadeinfoglu<-ddply(compiled, c("decade"), summarise,
                    N= sum(!is.na((GlutamicAcid))),
                    meanglu= mean(((GlutamicAcid)), na.rm = TRUE),
                    sd = sd((GlutamicAcid), na.rm=TRUE),
                    se = sd/sqrt(N),
                    var = var((GlutamicAcid), na.rm = TRUE))
as.numeric(decadeinfoGlutamicAcid$decade)


#phe does phe change over time?
compiled$logphe= log(compiled$Phenylalanine)
histogram(compiled$logphe)
shapiro.test(compiled$Phenylalanine)
shapiro.test(compiled$logphe)

normality<-glm(logphe~year, data = compiled)
res<-simulateResiduals(fittedModel = normality, n = 250)
res$scaledResiduals
plot(res)
testUniformity(res)

propphe<-(compiled$logphe/compiled$logsl)

#spatial checks for phe
spatial.compiled.phe<-glm(logphe~year, data=compiled)
simspatial.compiled.phe<-simulateResiduals(fittedModel = spatial.compiled.phe)
spatialtest.compiled.phe<-testSpatialAutocorrelation(simulationOutput = simspatial.compiled.phe,  x = compiled$longjitt, y = compiled$latjitt)
spatialtest.compiled.phe

#temporal checks for phe
temporal.phe<-glm(logphe~ year, data= compiled)
temporaltest.phe<-dwtest(temporal.phe, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) 
temporaltest.phe

#test for collinearity

b<-lm(Phenylalanine~ NPGO+  PDO_2lag, data= compiled)
ols_coll_diag(b)

phe1<-glmmTMB(log(Phenylalanine)~ NPGO + PDO_1lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(phe1)
phe2<-glmmTMB(log(Phenylalanine)~ NPGO + PDO_2lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(phe2)
phe3<-glmmTMB(log(Phenylalanine)~ NPGO + PDO_3lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(phe3)
phe0<-glmmTMB(log(Phenylalanine)~  scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(phe0)

glmmTable(phe1,"PDO", "phe1.html")
glmmTable(phe2,"PDO", "phe2.html")
glmmTable(phe3,"PDO", "phe3.html")
phe2nolog<-glmmTMB(Phenylalanine~ NPGO + PDO_2lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
#make the phe figure for publication

predict(phe2, compiled, allow.new.levels=TRUE)
ndphe<-compiled[1,]
ndphe$year<-"new"
ndphe_pop<-data.frame(year=compiled$year, sites=NA, NPGO = .2925, PDO_2lag = -0.06, sl = 127.33)
predictphe<-predict(phe2nolog, newdata=ndphe_pop, se.fit=TRUE)
as.data.frame(predictphe)
predictphe$year<-compiled$year
predictphe<-as.data.frame(predictphe)

#make the phefigure for publication


pheplot1<- ggplot() + geom_line(data =predictphe, aes(x = year, y = fit)) +
  geom_ribbon(data = predictphe, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data = compiled, aes(x = year, y = Phenylalanine), color="#00529c") +  xlab("Year collected") + ylab("Phenylalanine (in ‰)") +theme_classic()
ggsave("pheplot1.jpeg")



#phe decade info
decadeinfophe<-ddply(compiled, c("decade"), summarise,
                    N= sum(!is.na((Phenylalanine))),
                    meanglu= mean(((Phenylalanine)), na.rm = TRUE),
                    sd = sd((Phenylalanine), na.rm=TRUE),
                    se = sd/sqrt(N),
                    var = var((Phenylalanine), na.rm = TRUE))
as.numeric(decadeinfophe$decade)

#does time influence tp? yes, as years increases, tp is increasing

histogram(compiled$tp)
shapiro.test(compiled$tp)
normality<-glm(tp~year, data = compiled)
res<-simulateResiduals(fittedModel = normality, n = 250)
res$scaledResiduals
plot(res)
testUniformity(res)
#test for collinearity
c<-lm(tp~ NPGO+  PDO_2lag, data= compiled)
ols_coll_diag(c)

tp1<-glmmTMB(tp~ NPGO + PDO_1lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(tp1)
tp2<-glmmTMB(tp~ NPGO + PDO_2lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(tp2)
tp3<-glmmTMB(tp~ NPGO + PDO_3lag+ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(tp3)
tp0<-glmmTMB(tp~ scale(year) +(1|sites)  + sl, family = "gaussian", data= compiled)
summary(tp0)

glmmTable(tp1,"PDO", "tp1.html")
glmmTable(tp2,"PDO", "tp2.html")
glmmTable(tp3,"PDO", "tp3.html")
tp2inter<-glmmTMB(tp~ NPGO + PDO_2lag+ scale(year)*sl +(1|sites)  + sl, family = "gaussian", data= compiled)
#psthoc tp host size analyses 
cor.test(compiled$sl, compiled$year)

#make the tp figure for publication

predict(tp2, compiled, allow.new.levels=TRUE)
ndtp<-compiled[1,]
ndtp$year<-"new"
ndtp_pop<-data.frame(year=compiled$year, sites=NA, NPGO = .2925, PDO_2lag = -0.06, sl = 127.33) #random effects are set to NA, other effects I chose median of climate and mean for sl
tppredict<-predict(tp2, newdata=ndtp_pop, se.fit=TRUE)
as.data.frame(tppredict)
tppredict$year<-compiled$year
tppredict<-as.data.frame(tppredict)

#make the phefigure for publication


library(PupillometryR)

tplot1<- ggplot() + geom_line(data =tppredict, aes(x = year, y = fit)) +
  geom_ribbon(data = tppredict, aes(x = year, ymin = fit-se.fit, ymax = fit+se.fit), fill="#00529c", alpha=0.3) + geom_point(data= compiled, aes(x = year, y = tp), color="#00529c") +  xlab("") + ylab("Trophic position") + theme_classic() 

ggsave("tpplot1.jpeg")
library(PupillometryR)


samplesize<-  ggplot(compiled, aes (x = year, y = hostsp)) + geom_jitter(height = 0.15, width = 0.1, color = "#31a354")  + theme(panel.background = element_blank()) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position = "none")+ xlab("Year collected") +ylab("Replicates") 

samplegrob<-ggplotGrob(samplesize)

t


library (patchwork)
tplot1 / samplesize+   plot_layout(heights = c(3, 2)) + plot_annotation(tag_levels = 'A')
ggsave("tpplotsamplesize.jpg", height = 6, width = 5.5 )


par(mar=c(0,0,0,0))
tplot1 / samplesize +  plot_layout(heights = c(3, 2)) + plot_annotation(tag_levels = 'A')

#density<-ggplot(vio2, aes(x = year, y  = fishsp)) + 
 # geom_density_ridges(fill= NA, rel_min_height = 0.02, scale = 30) +  xlim(1930,2020) + theme_ridges() + theme_blank() + ylab("") +xlab("") + theme(axis.title.x=element_blank(),
#axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank())

#https://www.r-bloggers.com/2019/02/plots-within-plots-with-ggplot2-and-ggmap/


library(pwr)



mod<-lm(tp~PDO_2lag, data = compiled)
summary(mod)

#decade info
#decade info
library(dplyr)
library(tidyverse)
library(plyr)
library(reshape2)
decadeinfotp<-ddply(compiled, c("decade"), summarise,
                    N= sum(!is.na((tp))),
                    meantp= mean(((tp)), na.rm = TRUE),
                    sd = sd((tp), na.rm=TRUE),
                    se = sd/sqrt(N),
                    var = var((tp), na.rm = TRUE), 
                    median=(median(tp)),
                    min=(min(tp)),
                    max = (max(tp)))
as.numeric(decadeinfotp$decade)

