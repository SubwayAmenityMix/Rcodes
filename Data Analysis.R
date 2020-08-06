##########################################################################
################# DATA ANALYSIS ##########################################
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§#
pkgs <- c("maptools","rgdal","sp", "sf", "jpeg", "textreadr", "GEOmap", "RANN", "NISTunits",  "pracma", "celestial", "openxlsx", "data.table", "purrr",  "geonames",  "foreign", "tidyverse",  
          "ggplot2", "cowplot",  "psych",  "tidyimpute", "sjPlot", "ggpubr", "devtools", "car", "lattice",  "openintro",  
          "raster",  "spatialEco","xlsx", "spData",  "leaflet",  "geojsonio",  "rjson", "RJSONIO", "jsonlite","EconGeo", "finalfit",  
          "rstan", "boot",  "plm", "MatchIt", "MASS", "QuantPsyc", "reghelper", "jpeg", "ggfortify", "gapminder", "leaps", "usdm", "ggplot2", "glmnet", "MASS", "car", "corrplot", "texreg","MPV", "ebal")
sapply(pkgs, require, character.only = T) #load 
rm(pkgs)

######## Loading the data##########

setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona")
Barcecentroid_data <- read.csv("Barcecentroid_data.csv", header = TRUE)
Barcevenuehexagonbuffers_tot <- read.csv("Barcevenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")
Barce_PSM_dataset <- read.csv("Barce_PSM_dataset.csv", header = T)
Barcecentroid_data_metroandcontrol <- read.csv("Barcecentroid_data_metroandcontrol.csv", header = T)
write.csv(Barcecentroid_data, "Barcecentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna")
Viennacentroid_data <- read.csv("Viennacentroid_data.csv", header = TRUE)
Viennavenuehexagonbuffers_tot <- read.csv("Viennavenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")
Vienna_PSM_dataset <- read.csv("Vienna_PSM_dataset.csv", header = T)
Viennacentroid_data_metroandcontrol <- read.csv("Viennacentroid_data_metroandcontrol.csv", header = T)
write.csv(Viennacentroid_data, "Viennacentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome")
Romecentroid_data <- read.csv("Romecentroid_data.csv", header = TRUE)
Romevenuehexagonbuffers_tot <- read.csv("Romevenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
Rome_PSM_dataset <- read.csv("Rome_PSM_dataset.csv", header = T)
Romecentroid_data_metroandcontrol <- read.csv("Romecentroid_data_metroandcontrol.csv", header = T)
write.csv(Romecentroid_data, "Romecentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest")
Budapestcentroid_data <- read.csv("Budapestcentroid_data.csv", header = TRUE)
Budapestvenuehexagonbuffers_tot <- read.csv("Budapestvenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")
Budapest_PSM_dataset <- read.csv("Budapest_PSM_dataset.csv", header = T)
Budapestcentroid_data_metroandcontrol <- read.csv("Budapestcentroid_data_metroandcontrol.csv", header = T)
write.csv(Budapestcentroid_data, "Budapestcentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart")
Stuttgartcentroid_data <- read.csv("Stuttgartcentroid_data.csv", header = TRUE)
Stuttgartvenuehexagonbuffers_tot <- read.csv("Stuttgartvenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")
Stuttgart_PSM_dataset <- read.csv("Stuttgart_PSM_dataset.csv", header = T)
Stuttgartcentroid_data_metroandcontrol <- read.csv("Stuttgartcentroid_data_metroandcontrol.csv", header = T)
write.csv(Stuttgartcentroid_data, "Stuttgartcentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki")
Helsinkicentroid_data <- read.csv("Helsinkicentroid_data.csv", header = TRUE)
Helsinkivenuehexagonbuffers_tot <- read.csv("Helsinkivenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")
Helsinki_PSM_dataset <- read.csv("Helsinki_PSM_dataset.csv", header = T)
Helsinkicentroid_data_metroandcontrol <- read.csv("Helsinkicentroid_data_metroandcontrol.csv", header = T)
write.csv(Helsinkicentroid_data, "Helsinkicentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
Warsawcentroid_data <- read.csv("Warsawcentroid_data.csv", header = TRUE)
Warsawvenuehexagonbuffers_tot <- read.csv("Warsawvenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")
Warsaw_PSM_dataset <- read.csv("Warsaw_PSM_dataset.csv", header = T)
Warsawcentroid_data_metroandcontrol <- read.csv("Warsawcentroid_data_metroandcontrol.csv", header = T)
write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia")
Sofiacentroid_data <- read.csv("Sofiacentroid_data.csv", header = TRUE)
Sofiavenuehexagonbuffers_tot <- read.csv("Sofiavenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")
Sofia_PSM_dataset <- read.csv("Sofia_PSM_dataset.csv", header = T)
Sofiacentroid_data_metroandcontrol <- read.csv("Sofiacentroid_data_metroandcontrol.csv", header = T)
write.csv(Sofiacentroid_data, "Sofiacentroid_data.csv")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan")
Milancentroid_data <- read.csv("Milancentroid_data.csv", header = TRUE)
Milanvenuehexagonbuffers_tot <- read.csv("Milanvenuehexagonbuffers_tot.csv", header = TRUE)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")
Milan_PSM_dataset <- read.csv("Milan_PSM_dataset.csv", header = T)
Milancentroid_data_metroandcontrol <- read.csv("Milancentroid_data_metroandcontrol.csv", header = T)
write.csv(Milancentroid_data, "Milancentroid_data.csv")

###############################################################################
####################### Propensity score matching #############################
###############################################################################

#Rome
Romecentroidcontrolcases <- Romecentroid_data[which(Romecentroid_data$NEAR_DIST > 1500),]
Romecentroid_data_metroneighbouroods <- Romecentroid_data[which(Romecentroid_data$NEAR_DIST < 1200 & Romecentroid_data$Year_of_op > 2013),]
Romecentroid_data_metroneighbouroods$metroneighbourhood <- 1
Romecentroidcontrolcases$metroneighbourhood <- 0
Romecentroid_data_metroandcontrol <-  rbind(Romecentroid_data_metroneighbouroods, Romecentroidcontrolcases)

Romecentroid_data_metroandcontrol <- Romecentroid_data_metroandcontrol[!is.na(Romecentroid_data_metroandcontrol$streetdistance),]
Romecentroid_data_metroandcontrol <- Romecentroid_data_metroandcontrol[which(!is.na(Romecentroid_data_metroandcontrol$perc_unemployed) & !is.na(Romecentroid_data_metroandcontrol$population_2012) & !is.na(Romecentroid_data_metroandcontrol$centroid_distance_center) & !is.na(Romecentroid_data_metroandcontrol$rowsum_13 & !is.na(Romecentroid_data_metroandcontrol$streetdistance))),]
Romecentroid_data_metroandcontrol$streetdistance <- as.numeric(Romecentroid_data_metroandcontrol$streetdistance)
Romecentroid_data_metroandcontrol <- Romecentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST"  ,"Station_Na" ,  "Year_of_op" ,   "population_2012" , "population_2018",
                                                                          "centroid_distance_center","pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "entropy_13",   "rowsum_13" , "entropy_16", "rowsum_16" , "entropy_diff", 
                                                                          "density_increase" ,  "nr_Food_13_300mbuff" , "nr_ShopsServ_13_300mbuff",  "nr_ArtsEnter_13_300mbuff",  "nr_Nightlife_13_300mbuff"  ,"nr_Proff_13_300mbuff" ,
                                                                          "nr_Outdoors_13_300mbuff","nr_Food_1316_300mbuff", "nr_ShopsServ_1316_300mbuff",  "nr_ArtsEnter_1316_300mbuff",  "nr_Nightlife_1316_300mbuff" ,   
                                                                          "nr_Proff_1316_300mbuff", "nr_Outdoors_1316_300mbuff",  "perc_unemployed"  , "edu_tertiary_diplomas",
                                                                          "edu_highschool"  , "edu_lower_average" ,  "edu_elementary" , "foreigners","pop_dens","subcenter_amen_density","perc_foreign",
                                                                          "traveltime",  "streetdistance",  "metroneighbourhood", "multifunctionality_diff", "rowentropy_pre14")]


#Propensity Score Matching
summary(lm(metroneighbourhood ~  centroid_distance_center+ perc_unemployed + pop_2012_buffer + rowsum_13, Romecentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_distance_center + perc_unemployed+ pop_2012_buffer + rowsum_13, data = Romecentroid_data_metroandcontrol, method="nearest", ratio=1)
Rome_PSM_dataset <- match.data(match.it, distance ="pscore")

setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
sink("Rome PSM summary.txt")
summary(match.it)
sink()
sink("Rome PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~ centroid_distance_center+ perc_unemployed + pop_2012_buffer + rowsum_13, Romecentroid_data_metroandcontrol))
sink()
jpeg("Rome PSM results.jpeg", width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Rome PSM Histogram.jpeg", width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Rome PSM Jitter.jpeg", width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Rome PSM T-test.txt")
t.test(Rome_PSM_dataset$density_increase[Rome_PSM_dataset$metroneighbourhood== 1],Rome_PSM_dataset$density_increase[Rome_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Rome_PSM_dataset, "Rome_PSM_dataset.csv")

#Entropy Balancing Matching
Rome_eb <- ebal::ebalance(Treatment = Romecentroid_data_metroandcontrol$metroneighbourhood, X= Romecentroid_data_metroandcontrol[,c("centroid_distance_center", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "rowsum_13")], max.iterations = 200, print.level = 3)
Rome_eb_trim <- ebalance.trim(Rome_eb, max.weight = NULL,
              min.weight = 0, max.trim.iterations = 200,
              max.weight.increment = 0.92,
              min.weight.increment = 1.08,
              print.level = 0)


#covariate means of treated units
apply(Romecentroid_data_metroandcontrol[Romecentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_distance_center", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "rowsum_13", "perc_foreign")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Romecentroid_data_metroandcontrol[Romecentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_distance_center", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "rowsum_13", "perc_foreign")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Romecentroid_data_metroandcontrol[Romecentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_distance_center", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "rowsum_13", "perc_foreign")],2,weighted.mean, w= Rome_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Romecentroid_data_metroandcontrol[Romecentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_distance_center", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "rowsum_13", "perc_foreign")],2,weighted.mean, w= Rome_eb_trim$w)

Romecentroid_data_metroandcontrol$eb_weights[Romecentroid_data_metroandcontrol$metroneighbourhood==0] <- Rome_eb_trim$w
Romecentroid_data_metroandcontrol$eb_weights[Romecentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Romecentroid_data_metroandcontrol, "Romecentroid_data_metroandcontrol.csv")



##Barcelona
Barcecentroidcontrolcases <- Barcecentroid_data[which(Barcecentroid_data$distance_nearstation > 1200), ]
Barcecentroid_data_metroneighbouroods <- Barcecentroid_data[which(Barcecentroid_data$distance_nearstation < 1200 & Barcecentroid_data$station_Openingyear > 2015),]
Barcecentroid_data_metroneighbouroods$metroneighbourhood <- 1
Barcecentroidcontrolcases$metroneighbourhood <- 0
Barcecentroid_data_metroandcontrol <-  rbind(Barcecentroid_data_metroneighbouroods, Barcecentroidcontrolcases)

Barcecentroid_data_metroandcontrol <- Barcecentroid_data_metroandcontrol[!is.na(Barcecentroid_data_metroandcontrol$streetdistance),]
Barcecentroid_data_metroandcontrol <- Barcecentroid_data_metroandcontrol[!is.na(Barcecentroid_data_metroandcontrol$perc_unemployed_2015),]
Barcecentroid_data_metroandcontrol <- Barcecentroid_data_metroandcontrol[Barcecentroid_data_metroandcontrol$station_na != "Aeroport T2" & Barcecentroid_data_metroandcontrol$station_na !="Aeroport T1",]
Barcecentroid_data_metroandcontrol <- Barcecentroid_data_metroandcontrol[which(!is.na(Barcecentroid_data_metroandcontrol$pop_2012_buffer) & !is.na(Barcecentroid_data_metroandcontrol$centroid_distance_center) & !is.na(Barcecentroid_data_metroandcontrol$rowsum15) &  !is.na(Barcecentroid_data_metroandcontrol$multifunctionality_increase) ),]
Barcecentroid_data_metroandcontrol <- Barcecentroid_data_metroandcontrol[,c("ORIG_FID", "distance_nearstation", "station_na", "station_Openingyear", "population_2012", "population_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", 
                                                                            "perc_unemployed_2015", "centroid_distance_center", "entropy_pre15", "entropy18", "entropy_diff", "rowsum15", "rowsum18", "density_increase", "nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
                                                                            "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff", "nr_ArtsEnter_18_300mbuff", "nr_ArtsEnter_15_300mbuff",
                                                                            "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
                                                                            "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff", "rowsum15_2", "rowsum18_2", "density_increase2",
                                                                            "traveltime",  "streetdistance",  "rowentropy_16" , "rowentropy_19" ,   "multifunctionality_increase", "rowentropy_pre16" , "rowentropy_pre19" ,   "multifunctionality_diff", "metroneighbourhood") ]

summary(Barcecentroid_data$perc_unemployed_2015)
match.it <- matchit(metroneighbourhood ~ centroid_distance_center  + pop_2012_buffer + rowsum15, data = Barcecentroid_data_metroandcontrol, method="nearest",  ratio=1)

## PSM 
summary(lm(metroneighbourhood ~  centroid_distance_center  + pop_2012_buffer + perc_unemployed_2015+ rowsum15, Barcecentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_distance_center +perc_unemployed_2015 + pop_2012_buffer + rowsum15, data = Barcecentroid_data_metroandcontrol, method="nearest",  ratio=1)
Barce_PSM_dataset <- match.data(match.it, distance ="pscore")

summary(Barcecentroid_data_metroandcontrol$multifunctionality_diff)

setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")
sink("Barcelona PSM summary.txt")
summary(match.it)
sink()
sink("Barcelona PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_distance_center + perc_unemployed_2015 + population_2012 + rowsum15, Barcecentroid_data_metroandcontrol))
sink()
jpeg("Barcelona PSM results.jpeg", width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Barcelona PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Barcelona PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Barcelona PSM T-test.txt")
t.test(Barce_PSM_dataset$density_increase[Barce_PSM_dataset$metroneighbourhood== 1],Barce_PSM_dataset$density_increase[Barce_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Barce_PSM_dataset, "Barce_PSM_dataset.csv")

#Entropy Balancing Matching
Barce_eb <- ebal::ebalance(Treatment = Barcecentroid_data_metroandcontrol$metroneighbourhood, X= Barcecentroid_data_metroandcontrol[,c("centroid_distance_center"  , "pop_2012_buffer" , "rowsum15")], max.iterations = 200, print.level = 3)

Barce_eb <- ebal::ebalance(Treatment = Barcecentroid_data_metroandcontrol$metroneighbourhood, X= Barcecentroid_data_metroandcontrol[,c("centroid_distance_center"  , "pop_2012_buffer" , "perc_unemployed_2015", "rowsum15")], max.iterations = 200, print.level = 3)
Barce_eb_trim <- ebalance.trim(Barce_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Barcecentroid_data_metroandcontrol[Barcecentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_distance_center"  , "pop_2012_buffer" , "perc_unemployed_2015", "rowsum15")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Barcecentroid_data_metroandcontrol[Barcecentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_distance_center"  , "pop_2012_buffer" , "perc_unemployed_2015", "rowsum15")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Barcecentroid_data_metroandcontrol[Barcecentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_distance_center"  , "pop_2012_buffer" , "perc_unemployed_2015", "rowsum15")],2,weighted.mean, w= Barce_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Barcecentroid_data_metroandcontrol[Barcecentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_distance_center"  , "pop_2012_buffer" , "perc_unemployed_2015", "rowsum15")],2,weighted.mean, w= Barce_eb_trim$w)

Barcecentroid_data_metroandcontrol$eb_weights[Barcecentroid_data_metroandcontrol$metroneighbourhood==0] <- Barce_eb$w
Barcecentroid_data_metroandcontrol$eb_weights[Barcecentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Barcecentroid_data_metroandcontrol, "Barcecentroid_data_metroandcontrol.csv")



## Vienna
Viennacentroidcontrolcases <- Viennacentroid_data[which(Viennacentroid_data$NEAR_DIST > 1500), ]
Viennacentroid_data_metroneighbouroods <- Viennacentroid_data[which(Viennacentroid_data$NEAR_DIST < 1200 & Viennacentroid_data$station_openingyear > 2016),]
Viennacentroid_data_metroneighbouroods$metroneighbourhood <- 1
Viennacentroidcontrolcases$metroneighbourhood <- 0
Viennacentroid_data_metroandcontrol <-  rbind(Viennacentroid_data_metroneighbouroods, Viennacentroidcontrolcases)

Viennacentroid_data_metroandcontrol <- Viennacentroid_data_metroandcontrol[!is.na(Viennacentroid_data_metroandcontrol$streetdistance),]
Viennacentroid_data_metroandcontrol <- Viennacentroid_data_metroandcontrol[which(!is.na(Viennacentroid_data_metroandcontrol$perc_unemployed) & !is.na(Viennacentroid_data_metroandcontrol$pop_2012) & !is.na(Viennacentroid_data_metroandcontrol$distance_center) & !is.na(Viennacentroid_data_metroandcontrol$rowsum_16)),]
Viennacentroid_data_metroandcontrol <- Viennacentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "station_name", "station_openingyear", "pop_2012", "pop_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "perc_unemployed", "distance_center", "entropy_16", "entropy_19", "entropy_diff", "rowsum_16", "rowsum_19", "density_increase", "nr_Food_16_300mbuff", "nr_ShopsServ_16_300mbuff", "nr_ArtsEnter_16_300mbuff", 
                                                                              "nr_Nightlife_16_300mbuff", "nr_Proff_16_300mbuff", "nr_Outdoors_16_300mbuff",
                                                                              "nr_Food_1619_300mbuff", "nr_ShopsServ_1619_300mbuff", "nr_ArtsEnter_1619_300mbuff",
                                                                              "nr_Nightlife_1619_300mbuff", "nr_Proff_1619_300mbuff", "nr_Outdoors_1619_300mbuff",  "multifunctionality_increase", "rowentropy_16", "rowentropy_19", "multifunctionality_diff", "rowentropy_pre17", "rowentropy_pre20", "traveltime",  "streetdistance", "metroneighbourhood") ]

#PSM
summary(lm(metroneighbourhood ~  distance_center + perc_unemployed + pop_2012_buffer + rowsum_16, Viennacentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ distance_center + perc_unemployed + pop_2012_buffer + rowsum_16, data = Viennacentroid_data_metroandcontrol, method="nearest", discard = "both", ratio=1)
Vienna_PSM_dataset <- match.data(match.it, distance ="pscore")

summary(Viennacentroid_data_metroandcontrol$multifunctionality_increase)

setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")
sink("Vienna PSM summary.txt")
summary(match.it)
sink()
sink("Vienna PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  distance_center + perc_unemployed + pop_2012 + rowsum_16, Viennacentroid_data_metroandcontrol))
sink()
jpeg("Vienna PSM results.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it)
dev.off() 
jpeg("Vienna PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Vienna PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Vienna PSM T-test.txt")
t.test(Vienna_PSM_dataset$density_increase[Vienna_PSM_dataset$metroneighbourhood== 1],Vienna_PSM_dataset$density_increase[Vienna_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Vienna_PSM_dataset, "Vienna_PSM_dataset.csv")


#Entropy Balancing Matching
Vienna_eb <- ebal::ebalance(Treatment = Viennacentroid_data_metroandcontrol$metroneighbourhood, X= Viennacentroid_data_metroandcontrol[,c("distance_center" , "perc_unemployed" , "pop_2012_buffer" , "rowsum_16")], max.iterations = 200, print.level = 3)
Vienna_eb_trim <- ebalance.trim(Vienna_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Viennacentroid_data_metroandcontrol[Viennacentroid_data_metroandcontrol$metroneighbourhood==1, c("distance_center" , "perc_unemployed" , "pop_2012_buffer" , "rowsum_16")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Viennacentroid_data_metroandcontrol[Viennacentroid_data_metroandcontrol$metroneighbourhood==0, c("distance_center" , "perc_unemployed" , "pop_2012_buffer" , "rowsum_16")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Viennacentroid_data_metroandcontrol[Viennacentroid_data_metroandcontrol$metroneighbourhood==0,c("distance_center" , "perc_unemployed" , "pop_2012_buffer" , "rowsum_16")],2,weighted.mean, w= Vienna_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Viennacentroid_data_metroandcontrol[Viennacentroid_data_metroandcontrol$metroneighbourhood==0,c("distance_center" , "perc_unemployed" , "pop_2012_buffer" , "rowsum_16")],2,weighted.mean, w= Vienna_eb_trim$w)

Viennacentroid_data_metroandcontrol$eb_weights[Viennacentroid_data_metroandcontrol$metroneighbourhood==0] <- Vienna_eb_trim$w
Viennacentroid_data_metroandcontrol$eb_weights[Viennacentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Viennacentroid_data_metroandcontrol, "Viennacentroid_data_metroandcontrol.csv")




## Helsinki
Helsinkicentroid_data_metroneighbouroods <- Helsinkicentroid_data[which(Helsinkicentroid_data$NEAR_DIST < 1200 & Helsinkicentroid_data$Opened > 2016),]
Helsinkicentroidcontrolcases <- Helsinkicentroid_data[which(Helsinkicentroid_data$NEAR_DIST > 1500), ]
Helsinkicentroid_data_metroneighbouroods$metroneighbourhood <- 1
Helsinkicentroidcontrolcases$metroneighbourhood <- 0
Helsinkicentroid_data_metroandcontrol <-  rbind(Helsinkicentroid_data_metroneighbouroods, Helsinkicentroidcontrolcases)

Helsinkicentroid_data_metroandcontrol <- Helsinkicentroid_data_metroandcontrol[!is.na(Helsinkicentroid_data_metroandcontrol$streetdistance),]
Helsinkicentroid_data_metroandcontrol <- Helsinkicentroid_data_metroandcontrol[which(!is.na(Helsinkicentroid_data_metroandcontrol$unemployment_perc2012) & !is.na(Helsinkicentroid_data_metroandcontrol$pop_2012) & !is.na(Helsinkicentroid_data_metroandcontrol$centroid_dist_center) & !is.na(Helsinkicentroid_data_metroandcontrol$rowsum_16)),]
Helsinkicentroid_data_metroandcontrol <- Helsinkicentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "Station", "Opened", "pop_2012", "pop_2018","pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", 
                                                                                  "unemployment_perc2012", "centroid_dist_center", "entropy_16", "entropy_19", "entropy_diff", "rowsum_16", "rowsum_19", "density_increase", "nr_Food_16_300mbuff", "nr_ShopsServ_16_300mbuff", "nr_ArtsEnter_16_300mbuff", 
                                                                                  "nr_Nightlife_16_300mbuff", "nr_Proff_16_300mbuff", "nr_Outdoors_16_300mbuff",
                                                                                  "nr_Food_1619_300mbuff", "nr_ShopsServ_1619_300mbuff", "nr_ArtsEnter_1619_300mbuff",
                                                                                  "nr_Nightlife_1619_300mbuff", "nr_Proff_1619_300mbuff", "nr_Outdoors_1619_300mbuff",
                                                                                  "traveltime",  "streetdistance",  "rowentropy_pre17" , "rowentropy_pre20" ,   "multifunctionality_diff", "metroneighbourhood") ]


#PSM
summary(lm(metroneighbourhood ~  centroid_dist_center + unemployment_perc2012 + pop_2012_buffer + rowsum_16, Helsinkicentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_dist_center + unemployment_perc2012 +pop_2012_buffer + rowsum_16, data = Helsinkicentroid_data_metroandcontrol, method="nearest",  ratio=1)
Helsinki_PSM_dataset <- match.data(match.it)

summary(Helsinkicentroid_data_metroandcontrol$multifunctionality_diff)

setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")
sink("Helsinki PSM summary.txt")
summary(match.it)
sink()
sink("Helsinki PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_dist_center + unemployment_perc2012 + pop_2012 + rowsum_16, Helsinkicentroid_data_metroandcontrol))
sink()
jpeg("Helsinki PSM results.jpeg",  width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Helsinki PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Helsinki PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Helsinki PSM T-test.txt")
t.test(Helsinki_PSM_dataset$density_increase[Helsinki_PSM_dataset$metroneighbourhood== 1],Helsinki_PSM_dataset$density_increase[Helsinki_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Helsinki_PSM_dataset, "Helsinki_PSM_dataset.csv")

#Entropy Balancing Matching
Helsinki_eb <- ebal::ebalance(Treatment = Helsinkicentroid_data_metroandcontrol$metroneighbourhood, X= Helsinkicentroid_data_metroandcontrol[,c("centroid_dist_center" , "unemployment_perc2012" , "pop_2012_buffer" , "rowsum_16")], max.iterations = 200, print.level = 3)
Helsinki_eb_trim <- ebalance.trim(Helsinki_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Helsinkicentroid_data_metroandcontrol[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_dist_center" , "unemployment_perc2012" , "pop_2012_buffer" , "rowsum_16")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Helsinkicentroid_data_metroandcontrol[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_dist_center" , "unemployment_perc2012" , "pop_2012_buffer" , "rowsum_16")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Helsinkicentroid_data_metroandcontrol[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "unemployment_perc2012" , "pop_2012_buffer" , "rowsum_16")],2,weighted.mean, w= Helsinki_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Helsinkicentroid_data_metroandcontrol[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "unemployment_perc2012" , "pop_2012_buffer" , "rowsum_16")],2,weighted.mean, w= Helsinki_eb_trim$w)

Helsinkicentroid_data_metroandcontrol$eb_weights[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==0] <- Helsinki_eb_trim$w
Helsinkicentroid_data_metroandcontrol$eb_weights[Helsinkicentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Helsinkicentroid_data_metroandcontrol, "Helsinkicentroid_data_metroandcontrol.csv")



## Budapest
Budapestcentroidcontrolcases <- Budapestcentroid_data[which(Budapestcentroid_data$NEAR_DIST > 1500), ]
Budapestcentroid_data_metroneighbouroods <- Budapestcentroid_data[which(Budapestcentroid_data$NEAR_DIST < 1200 & Budapestcentroid_data$Year_of_op > 2013),]
Budapestcentroid_data_metroneighbouroods$metroneighbourhood <- 1
Budapestcentroidcontrolcases$metroneighbourhood <- 0
Budapestcentroid_data_metroandcontrol <-  rbind(Budapestcentroid_data_metroneighbouroods, Budapestcentroidcontrolcases)

Budapestcentroid_data_metroandcontrol <- Budapestcentroid_data_metroandcontrol[!is.na(Budapestcentroid_data_metroandcontrol$streetdistance),]
Budapestcentroid_data_metroandcontrol <- Budapestcentroid_data_metroandcontrol[which(!is.na(Budapestcentroid_data_metroandcontrol$percentage_unemployed) & !is.na(Budapestcentroid_data_metroandcontrol$pop_2012) & !is.na(Budapestcentroid_data_metroandcontrol$distance_center) & !is.na(Budapestcentroid_data_metroandcontrol$rowsum14)),]
Budapestcentroid_data_metroandcontrol <- Budapestcentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "Station", "Year_of_op", "pop_2012", "pop_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "percentage_unemployed", "distance_center", "entropy_pre14", "entropy_pre17", "entropy_diff", "rowsum14", "rowsum17", "density_increase", "nr_ArtsEnter_1316_300mbuff", "nr_Nightlife_1316_300mbuff", "nr_Proff_1316_300mbuff" , 
                                                                                  "nr_Outdoors_1316_300mbuff", "nr_ShopsServ_1316_300mbuff" , "nr_Food_1316_300mbuff" ,
                                                                                  "nr_Proff_13_300mbuff", "nr_Outdoors_13_300mbuff",  "nr_ShopsServ_13_300mbuff", 
                                                                                  "nr_ArtsEnter_13_300mbuff",  "nr_Food_13_300mbuff", "nr_Nightlife_13_300mbuff",   "multifunctionality_diff", "rowentropy_pre14", "rowentropy_pre17", "traveltime",  "streetdistance", "metroneighbourhood") ]

#PSM
summary(lm(metroneighbourhood ~  distance_center + percentage_unemployed + pop_2012_buffer + rowsum14, Budapestcentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ distance_center + percentage_unemployed + pop_2012_buffer + rowsum14, data = Budapestcentroid_data_metroandcontrol, method="nearest", discard = "both", ratio=1)
Budapest_PSM_dataset <- match.data(match.it, distance ="pscore")


setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")
sink("Budapest PSM summary.txt")
summary(match.it)
sink()
sink("Budapest PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  distance_center + percentage_unemployed + pop_2012 + rowsum14, Budapestcentroid_data_metroandcontrol))
sink()
jpeg("Budapest PSM results.jpeg", width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Budapest PSM Histogram.jpeg", width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Budapest PSM Jitter.jpeg", width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Budapest PSM T-test.txt")
t.test(Budapest_PSM_dataset$density_increase[Budapest_PSM_dataset$metroneighbourhood== 1],Budapest_PSM_dataset$density_increase[Budapest_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Budapest_PSM_dataset, "Budapest_PSM_dataset.csv")

#Entropy Balancing Matching
Budapest_eb <- ebal::ebalance(Treatment = Budapestcentroid_data_metroandcontrol$metroneighbourhood, X= Budapestcentroid_data_metroandcontrol[,c("distance_center" , "percentage_unemployed" , "pop_2012_buffer" , "rowsum14")], max.iterations = 200, print.level = 3)
Budapest_eb_trim <- ebalance.trim(Budapest_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Budapestcentroid_data_metroandcontrol[Budapestcentroid_data_metroandcontrol$metroneighbourhood==1, c("distance_center" , "percentage_unemployed" , "pop_2012_buffer" , "rowsum14")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Budapestcentroid_data_metroandcontrol[Budapestcentroid_data_metroandcontrol$metroneighbourhood==0, c("distance_center" , "percentage_unemployed" , "pop_2012_buffer" , "rowsum14")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Budapestcentroid_data_metroandcontrol[Budapestcentroid_data_metroandcontrol$metroneighbourhood==0,c("distance_center" , "percentage_unemployed" , "pop_2012_buffer" , "rowsum14")],2,weighted.mean, w= Budapest_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Budapestcentroid_data_metroandcontrol[Budapestcentroid_data_metroandcontrol$metroneighbourhood==0,c("distance_center" , "percentage_unemployed" , "pop_2012_buffer" , "rowsum14")],2,weighted.mean, w= Budapest_eb_trim$w)

Budapestcentroid_data_metroandcontrol$eb_weights[Budapestcentroid_data_metroandcontrol$metroneighbourhood==0] <- Budapest_eb_trim$w
Budapestcentroid_data_metroandcontrol$eb_weights[Budapestcentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Budapestcentroid_data_metroandcontrol, "Budapestcentroid_data_metroandcontrol.csv")



## Stuttgart
Stuttgartcentroidcontrolcases <- Stuttgartcentroid_data[which(Stuttgartcentroid_data$NEAR_DIST > 1500), ]
Stuttgartcentroid_data_metroneighbouroods <- Stuttgartcentroid_data[which(Stuttgartcentroid_data$NEAR_DIST < 1200 & Stuttgartcentroid_data$openingyea > 2013),]
Stuttgartcentroid_data_metroneighbouroods$metroneighbourhood <- 1
Stuttgartcentroidcontrolcases$metroneighbourhood <- 0
Stuttgartcentroid_data_metroandcontrol <-  rbind(Stuttgartcentroid_data_metroneighbouroods, Stuttgartcentroidcontrolcases)

Stuttgartcentroid_data_metroandcontrol <- Stuttgartcentroid_data_metroandcontrol[!is.na(Stuttgartcentroid_data_metroandcontrol$streetdistance),]
Stuttgartcentroid_data_metroandcontrol <- Stuttgartcentroid_data_metroandcontrol[which(!is.na(Stuttgartcentroid_data_metroandcontrol$perc_unemployed12) & !is.na(Stuttgartcentroid_data_metroandcontrol$pop_2012) & !is.na(Stuttgartcentroid_data_metroandcontrol$centroid_dist_center) & !is.na(Stuttgartcentroid_data_metroandcontrol$rowsum_15)),]
Stuttgartcentroid_data_metroandcontrol <- Stuttgartcentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "Name", "openingyea", "pop_2012", "pop_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "perc_unemployed12", "centroid_dist_center", "entropy_15", "entropy_18", "entropy_diff", "rowsum_15", "rowsum_18", "density_increase", "nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
                                                                                    "centroid_lat", "centroid_lon", "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
                                                                                    "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
                                                                                    "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff", "multifunctionality_diff", "rowentropy_pre16", "rowentropy_pre19", "traveltime",  "streetdistance", "metroneighbourhood") ]


#PSM
summary(lm(metroneighbourhood ~ centroid_dist_center + perc_unemployed12 + pop_2012_buffer + rowsum_15, Stuttgartcentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_dist_center + perc_unemployed12 + pop_2012_buffer + rowsum_15, data = Stuttgartcentroid_data_metroandcontrol, method="nearest", discard = "both", ratio=1)
Stuttgart_PSM_dataset <- match.data(match.it, distance ="pscore")

summary(Stuttgartcentroid_data_metroandcontrol$multifunctionality_increase)

setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")
sink("Stuttgart PSM summary.txt")
summary(match.it)
sink()
sink("Stuttgart PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_dist_center + perc_unemployed12 + pop_2012 + rowsum_15, Stuttgartcentroid_data_metroandcontrol))
sink()
jpeg("Stuttgart PSM results.jpeg", width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Stuttgart PSM Histogram.jpeg", width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Stuttgart PSM Jitter.jpeg", width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Stuttgart PSM T-test.txt")
t.test(Stuttgart_PSM_dataset$density_increase[Stuttgart_PSM_dataset$metroneighbourhood== 1],Stuttgart_PSM_dataset$density_increase[Stuttgart_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Stuttgart_PSM_dataset, "Stuttgart_PSM_dataset.csv")

#Entropy Balancing Matching
Stuttgart_eb <- ebal::ebalance(Treatment = Stuttgartcentroid_data_metroandcontrol$metroneighbourhood, X= Stuttgartcentroid_data_metroandcontrol[,c("centroid_dist_center" , "perc_unemployed12" , "pop_2012_buffer" , "rowsum_15")], max.iterations = 200, print.level = 3)
Stuttgart_eb_trim <- ebalance.trim(Stuttgart_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Stuttgartcentroid_data_metroandcontrol[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_dist_center" , "perc_unemployed12" , "pop_2012_buffer" , "rowsum_15")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Stuttgartcentroid_data_metroandcontrol[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_dist_center" , "perc_unemployed12" , "pop_2012_buffer" , "rowsum_15")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Stuttgartcentroid_data_metroandcontrol[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "perc_unemployed12" , "pop_2012_buffer" , "rowsum_15")],2,weighted.mean, w= Stuttgart_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Stuttgartcentroid_data_metroandcontrol[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "perc_unemployed12" , "pop_2012_buffer" , "rowsum_15")],2,weighted.mean, w= Stuttgart_eb_trim$w)

Stuttgartcentroid_data_metroandcontrol$eb_weights[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==0] <- Stuttgart_eb_trim$w
Stuttgartcentroid_data_metroandcontrol$eb_weights[Stuttgartcentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Stuttgartcentroid_data_metroandcontrol, "Stuttgartcentroid_data_metroandcontrol.csv")


## Warsaw
Warsawcentroid_data_metroneighbouroods <- Warsawcentroid_data[which(Warsawcentroid_data$NEAR_DIST < 1200 & Warsawcentroid_data$station_openingyear > 2014 & Warsawcentroid_data$station_openingyear < 2016),]
Warsawcentroidcontrolcases <- Warsawcentroid_data[which(Warsawcentroid_data$NEAR_DIST > 1500), ]
Warsawcentroid_data_metroneighbouroods$metroneighbourhood <- 1
Warsawcentroidcontrolcases$metroneighbourhood <- 0
Warsawcentroid_data_metroandcontrol <-  rbind(Warsawcentroid_data_metroneighbouroods, Warsawcentroidcontrolcases)

Warsawcentroid_data_metroandcontrol <- Warsawcentroid_data_metroandcontrol[!is.na(Warsawcentroid_data_metroandcontrol$streetdistance),]
Warsawcentroid_data_metroandcontrol <- Warsawcentroid_data_metroandcontrol[which(!is.na(Warsawcentroid_data_metroandcontrol$unemployed_2010) & !is.na(Warsawcentroid_data_metroandcontrol$pop_2012) & !is.na(Warsawcentroid_data_metroandcontrol$centroid_dist_center) & !is.na(Warsawcentroid_data_metroandcontrol$rowsum_14)),]
Warsawcentroid_data_metroandcontrol <- Warsawcentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "stationid", "station_openingyear", "pop_2012", "pop_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "unemployed_2010", "centroid_dist_center", "entropy_14", "entropy18", "entropy_diff", "rowsum_14", "rowsum18", "density_increase", "nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
                                                                                  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
                                                                                  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
                                                                                  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff",  "multifunctionality_increase", "rowentropy_15", "rowentropy_18", "multifunctionality_diff", "rowentropy_pre15", "rowentropy_pre19","traveltime",  "streetdistance", "metroneighbourhood") ]

#PSM
summary(lm(metroneighbourhood ~  centroid_dist_center + unemployed_2010 + pop_2012_buffer + rowsum_14, Warsawcentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_dist_center  + unemployed_2010 + pop_2012_buffer + rowsum_14, data = Warsawcentroid_data_metroandcontrol, method="nearest",  ratio=1)
Warsaw_PSM_dataset <- match.data(match.it)

summary(Warsawcentroid_data_metroandcontrol$multifunctionality_diff)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")
sink("Warsaw PSM summary.txt")
summary(match.it)
sink()
sink("Warsaw PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_dist_center + unemployed_2010 + pop_2012 + rowsum_14, Warsawcentroid_data_metroandcontrol))
sink()
jpeg("Warsaw PSM results.jpeg",  width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Warsaw PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Warsaw PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Warsaw PSM T-test.txt")
t.test(Warsaw_PSM_dataset$density_increase[Warsaw_PSM_dataset$metroneighbourhood== 1],Warsaw_PSM_dataset$density_increase[Warsaw_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Warsaw_PSM_dataset, "Warsaw_PSM_dataset.csv")

#Entropy Balancing Matching
Warsaw_eb <- ebal::ebalance(Treatment = Warsawcentroid_data_metroandcontrol$metroneighbourhood, X= Warsawcentroid_data_metroandcontrol[,c("centroid_dist_center" , "unemployed_2010" , "pop_2012_buffer" , "rowsum_14")], max.iterations = 200, print.level = 3)
Warsaw_eb_trim <- ebalance.trim(Warsaw_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Warsawcentroid_data_metroandcontrol[Warsawcentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_dist_center" , "unemployed_2010" , "pop_2012_buffer" , "rowsum_14")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Warsawcentroid_data_metroandcontrol[Warsawcentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_dist_center" , "unemployed_2010" , "pop_2012_buffer" , "rowsum_14")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Warsawcentroid_data_metroandcontrol[Warsawcentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "unemployed_2010" , "pop_2012_buffer" , "rowsum_14")],2,weighted.mean, w= Warsaw_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Warsawcentroid_data_metroandcontrol[Warsawcentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "unemployed_2010" , "pop_2012_buffer" , "rowsum_14")],2,weighted.mean, w= Warsaw_eb_trim$w)

Warsawcentroid_data_metroandcontrol$eb_weights[Warsawcentroid_data_metroandcontrol$metroneighbourhood==0] <- Warsaw_eb$w
Warsawcentroid_data_metroandcontrol$eb_weights[Warsawcentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Warsawcentroid_data_metroandcontrol, "Warsawcentroid_data_metroandcontrol.csv")



## Sofia
Sofiacentroid_data_metroneighbouroods <- Sofiacentroid_data[which(Sofiacentroid_data$NEAR_DIST < 1200 & Sofiacentroid_data$Openingyear> 2014),]
Sofiacentroidcontrolcases <- Sofiacentroid_data[which(Sofiacentroid_data$NEAR_DIST > 1500), ]
Sofiacentroid_data_metroneighbouroods$metroneighbourhood <- 1
Sofiacentroidcontrolcases$metroneighbourhood <- 0
Sofiacentroid_data_metroandcontrol <-  rbind(Sofiacentroid_data_metroneighbouroods, Sofiacentroidcontrolcases)

Sofiacentroid_data_metroandcontrol <- Sofiacentroid_data_metroandcontrol[!is.na(Sofiacentroid_data_metroandcontrol$streetdistance),]
Sofiacentroid_data_metroandcontrol <- Sofiacentroid_data_metroandcontrol[which(!is.na(Sofiacentroid_data_metroandcontrol$pop_2012) & !is.na(Sofiacentroid_data_metroandcontrol$centroid_dist_center) & !is.na(Sofiacentroid_data_metroandcontrol$rowsum_14)),]
Sofiacentroid_data_metroandcontrol <- Sofiacentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "Station_latin", "Openingyear", "pop_2012", "pop_2018", "pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer","centroid_dist_center", "entropy_14", "entropy18", "entropy_diff", "rowsum_14", "rowsum18", "density_increase", "nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
                                                                                  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
                                                                                  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
                                                                                  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff", "multifunctionality_diff", "rowentropy_pre15", "rowentropy_pre18", "traveltime",  "streetdistance", "metroneighbourhood") ]

#PSM
summary(lm(metroneighbourhood ~  centroid_dist_center + pop_2012_buffer + rowsum_14, Sofiacentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_dist_center + pop_2012_buffer  + rowsum_14, data = Sofiacentroid_data_metroandcontrol, method="nearest",  ratio=1)
Sofia_PSM_dataset <- match.data(match.it)

summary(Sofiacentroid_data_metroandcontrol$multifunctionality_diff)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")
sink("Sofia PSM summary.txt")
summary(match.it)
sink()
sink("Sofia PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_dist_center + pop_2012 + rowsum_14, Sofiacentroid_data_metroandcontrol))
sink()
jpeg("Sofia PSM results.jpeg",  width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Sofia PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Sofia PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Sofia PSM T-test.txt")
t.test(Sofia_PSM_dataset$density_increase[Sofia_PSM_dataset$metroneighbourhood== 1],Sofia_PSM_dataset$density_increase[Sofia_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Sofia_PSM_dataset, "Sofia_PSM_dataset.csv")

#Entropy Balancing Matching
Sofia_eb <- ebal::ebalance(Treatment = Sofiacentroid_data_metroandcontrol$metroneighbourhood, X= Sofiacentroid_data_metroandcontrol[,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14")], max.iterations = 200, print.level = 3)
Sofia_eb_trim <- ebalance.trim(Sofia_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Sofiacentroid_data_metroandcontrol[Sofiacentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Sofiacentroid_data_metroandcontrol[Sofiacentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Sofiacentroid_data_metroandcontrol[Sofiacentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14")],2,weighted.mean, w= Sofia_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Sofiacentroid_data_metroandcontrol[Sofiacentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14")],2,weighted.mean, w= Sofia_eb_trim$w)

Sofiacentroid_data_metroandcontrol$eb_weights[Sofiacentroid_data_metroandcontrol$metroneighbourhood==0] <- Sofia_eb_trim$w
Sofiacentroid_data_metroandcontrol$eb_weights[Sofiacentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Sofiacentroid_data_metroandcontrol, "Sofiacentroid_data_metroandcontrol.csv")



## Milan
Milancentroid_data_metroneighbouroods <- Milancentroid_data[which(Milancentroid_data$NEAR_DIST < 1200 & Milancentroid_data$Openingyear > 2013),]
Milancentroidcontrolcases <- Milancentroid_data[which(Milancentroid_data$NEAR_DIST > 1500), ]
Milancentroid_data_metroneighbouroods$metroneighbourhood <- 1
Milancentroidcontrolcases$metroneighbourhood <- 0
Milancentroid_data_metroandcontrol <-  rbind(Milancentroid_data_metroneighbouroods, Milancentroidcontrolcases)

Milancentroid_data_metroandcontrol <- Milancentroid_data_metroandcontrol[!is.na(Milancentroid_data_metroandcontrol$streetdistance),]
Milancentroid_data_metroandcontrol$streetdistance <- as.numeric(Milancentroid_data_metroandcontrol$streetdistance)
Milancentroid_data_metroandcontrol <- Milancentroid_data_metroandcontrol[which( !is.na(Milancentroid_data_metroandcontrol$perc_unemployed2011)),]
Milancentroid_data_metroandcontrol <- Milancentroid_data_metroandcontrol[which(!is.na(Milancentroid_data_metroandcontrol$pop_2012_buffer) & !is.na(Milancentroid_data_metroandcontrol$centroid_dist_center) & !is.na(Milancentroid_data_metroandcontrol$rowsum_14)),]
Milancentroid_data_metroandcontrol <- Milancentroid_data_metroandcontrol[,c("ORIG_FID", "NEAR_DIST", "MetroStationname", "Openingyear", "pop_2012", "pop_2018","pop_2012_buffer", "pop_2018_buffer", "popchange", "popchange_buffer", "centroid_dist_center", "entropy_14", "entropy17", "entropy_diff", "rowsum_14", "rowsum17", "density_increase", "nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
                                                                            "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff", "perc_unemployed2011", 
                                                                            "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
                                                                            "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff",  "multifunctionality_diff", "rowentropy_pre15", "rowentropy_pre18", "traveltime",  "streetdistance", "metroneighbourhood") ]


#PSM
summary(lm(metroneighbourhood ~  centroid_dist_center + pop_2012_buffer + rowsum_14 + perc_unemployed2011, Milancentroid_data_metroandcontrol))
match.it <- matchit(metroneighbourhood ~ centroid_dist_center + perc_unemployed2011 + pop_2012_buffer + rowsum_14, data = Milancentroid_data_metroandcontrol, method="nearest",  ratio=1)
Milan_PSM_dataset <- match.data(match.it)


summary(Milancentroid_data_metroandcontrol$multifunctionality_increase)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")
sink("Milan PSM summary.txt")
summary(match.it)
sink()
sink("Milan PSM treatment prediction.txt")
summary(lm(metroneighbourhood ~  centroid_dist_center + pop_2012 + rowsum_14, Milancentroid_data_metroandcontrol))
sink()
jpeg("Milan PSM results.jpeg",  width = 8, height = 8, units = "in", res = 1000)
g <- plot(match.it)
dev.off() 
jpeg("Milan PSM Histogram.jpeg",  width = 6, height = 6, units = "in", res = 1000)
plot(match.it, type = "hist")
dev.off() 
jpeg("Milan PSM Jitter.jpeg",  width = 8, height = 8, units = "in", res = 1000)
plot(match.it, type = "jitter", interactive = F)
dev.off() 
sink("Milan PSM T-test.txt")
t.test(Milan_PSM_dataset$density_increase[Milan_PSM_dataset$metroneighbourhood== 1],Milan_PSM_dataset$density_increase[Milan_PSM_dataset$metroneighbourhood==0],paired=TRUE)
sink()
write.csv(Milan_PSM_dataset, "Milan_PSM_dataset.csv")

#Entropy Balancing Matching
Milan_eb <- ebal::ebalance(Treatment = Milancentroid_data_metroandcontrol$metroneighbourhood, X= Milancentroid_data_metroandcontrol[,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14" , "perc_unemployed2011")], max.iterations = 200, print.level = 3)
Milan_eb_trim <- ebalance.trim(Milan_eb, max.weight = NULL,
                              min.weight = 0, max.trim.iterations = 200,
                              max.weight.increment = 0.92,
                              min.weight.increment = 1.08,
                              print.level = 0)


#covariate means of treated units
apply(Milancentroid_data_metroandcontrol[Milancentroid_data_metroandcontrol$metroneighbourhood==1, c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14" , "perc_unemployed2011")],2,mean)
#covariate means of eligible untreated units before reweighting
apply(Milancentroid_data_metroandcontrol[Milancentroid_data_metroandcontrol$metroneighbourhood==0, c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14" , "perc_unemployed2011")],2,mean)
#covariate means of eligible untreated units after entropy balancing
apply(Milancentroid_data_metroandcontrol[Milancentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14" , "perc_unemployed2011")],2,weighted.mean, w= Milan_eb$w)
#covariate means of eligible untreated units after entropy balancing and trimming
apply(Milancentroid_data_metroandcontrol[Milancentroid_data_metroandcontrol$metroneighbourhood==0,c("centroid_dist_center" , "pop_2012_buffer" , "rowsum_14" , "perc_unemployed2011")],2,weighted.mean, w= Milan_eb_trim$w)

Milancentroid_data_metroandcontrol$eb_weights[Milancentroid_data_metroandcontrol$metroneighbourhood==0] <- Milan_eb_trim$w
Milancentroid_data_metroandcontrol$eb_weights[Milancentroid_data_metroandcontrol$metroneighbourhood==1] <- 1
write.csv(Milancentroid_data_metroandcontrol, "Milancentroid_data_metroandcontrol.csv")




#################################################################################################
####################### Create Dummies for Spatial Bufferrings arround stations ##########################
################################################################################################

## Rome
Romecentroid_data_onlymetroneighbour <- Romecentroid_data[which(Romecentroid_data$Year_of_op > 2014 & Romecentroid_data$NEAR_DIST < 1201),]
Romecentroid_data_onlymetroneighbour$buffer400m <- 0
Romecentroid_data_onlymetroneighbour$buffer800m <- 0
Romecentroid_data_onlymetroneighbour$buffer1200m <- 0
Romecentroid_data_onlymetroneighbour$buffer400m[Romecentroid_data_onlymetroneighbour$NEAR_DIST <= 400] <- 1
Romecentroid_data_onlymetroneighbour$buffer800m[Romecentroid_data_onlymetroneighbour$NEAR_DIST > 400 & Romecentroid_data_onlymetroneighbour$NEAR_DIST <= 800] <- 1
Romecentroid_data_onlymetroneighbour$buffer1200m[Romecentroid_data_onlymetroneighbour$NEAR_DIST > 800 & Romecentroid_data_onlymetroneighbour$NEAR_DIST <= 1200] <- 1


## Barcelona
Barcecentroid_data_onlymetroneighbour <- Barcecentroid_data[which(Barcecentroid_data$station_Openingyear > 2015 & Barcecentroid_data$distance_nearstation < 1201),]
Barcecentroid_data_onlymetroneighbour$buffer400m <- 0
Barcecentroid_data_onlymetroneighbour$buffer800m <- 0
Barcecentroid_data_onlymetroneighbour$buffer1200m <- 0
Barcecentroid_data_onlymetroneighbour$buffer400m[Barcecentroid_data_onlymetroneighbour$distance_nearstation <= 400] <- 1
Barcecentroid_data_onlymetroneighbour$buffer800m[Barcecentroid_data_onlymetroneighbour$distance_nearstation > 400 & Barcecentroid_data_onlymetroneighbour$distance_nearstation <= 800] <- 1
Barcecentroid_data_onlymetroneighbour$buffer1200m[Barcecentroid_data_onlymetroneighbour$distance_nearstation > 800 & Barcecentroid_data_onlymetroneighbour$distance_nearstation <= 1200] <- 1


## Helsinki
Helsinkicentroid_data_onlymetroneighbour <- Helsinkicentroid_data[which(Helsinkicentroid_data$Opened > 2016 & Helsinkicentroid_data$NEAR_DIST < 1201),]
Helsinkicentroid_data_onlymetroneighbour$buffer400m <- 0
Helsinkicentroid_data_onlymetroneighbour$buffer800m <- 0
Helsinkicentroid_data_onlymetroneighbour$buffer1200m <- 0
Helsinkicentroid_data_onlymetroneighbour$buffer400m[Helsinkicentroid_data_onlymetroneighbour$NEAR_DIST <= 400] <- 1
Helsinkicentroid_data_onlymetroneighbour$buffer800m[Helsinkicentroid_data_onlymetroneighbour$NEAR_DIST > 400 & Helsinkicentroid_data_onlymetroneighbour$NEAR_DIST <= 800] <- 1
Helsinkicentroid_data_onlymetroneighbour$buffer1200m[Helsinkicentroid_data_onlymetroneighbour$NEAR_DIST > 800 & Helsinkicentroid_data_onlymetroneighbour$NEAR_DIST <= 1200] <- 1


## Budapest
Budapestcentroid_data_onlymetroneighbour <- Budapestcentroid_data[which(Budapestcentroid_data$Year_of_op > 2013 & Budapestcentroid_data$NEAR_DIST < 1201),]
Budapestcentroid_data_onlymetroneighbour$buffer400m <- 0
Budapestcentroid_data_onlymetroneighbour$buffer800m <- 0
Budapestcentroid_data_onlymetroneighbour$buffer1200m <- 0
Budapestcentroid_data_onlymetroneighbour$buffer400m[Budapestcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400] <- 1
Budapestcentroid_data_onlymetroneighbour$buffer800m[Budapestcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Budapestcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800] <- 1
Budapestcentroid_data_onlymetroneighbour$buffer1200m[Budapestcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Budapestcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200] <- 1


## Vienna
Viennacentroid_data_onlymetroneighbour <- Viennacentroid_data[which(Viennacentroid_data$station_openingyear > 2016 & Viennacentroid_data$NEAR_DIST < 1201),]

Viennacentroid_data_onlymetroneighbour$buffer400m <- 0
Viennacentroid_data_onlymetroneighbour$buffer800m <- 0
Viennacentroid_data_onlymetroneighbour$buffer1200m <- 0

for (i in 1:nrow(Viennacentroid_data_onlymetroneighbour)){
  if(Viennacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400){
    Viennacentroid_data_onlymetroneighbour$buffer400m[i] <- 1
  }
}
for (i in 1:nrow(Viennacentroid_data_onlymetroneighbour)){
  if(Viennacentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Viennacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800){
    Viennacentroid_data_onlymetroneighbour$buffer800m[i] <- 1
  }
}
for (i in 1:nrow(Viennacentroid_data_onlymetroneighbour)){
  if(Viennacentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Viennacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200){
    Viennacentroid_data_onlymetroneighbour$buffer1200m[i] <- 1
  }
}

## Stuttgart
Stuttgartcentroid_data_onlymetroneighbour <- Stuttgartcentroid_data[which(Stuttgartcentroid_data$openingyea > 2015 & Stuttgartcentroid_data$NEAR_DIST < 1201),]

Stuttgartcentroid_data_onlymetroneighbour$buffer400m <- 0
Stuttgartcentroid_data_onlymetroneighbour$buffer800m <- 0
Stuttgartcentroid_data_onlymetroneighbour$buffer1200m <- 0

for (i in 1:nrow(Stuttgartcentroid_data_onlymetroneighbour)){
  if(Stuttgartcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400){
    Stuttgartcentroid_data_onlymetroneighbour$buffer400m[i] <- 1
  }
}
for (i in 1:nrow(Stuttgartcentroid_data_onlymetroneighbour)){
  if(Stuttgartcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Stuttgartcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800){
    Stuttgartcentroid_data_onlymetroneighbour$buffer800m[i] <- 1
  }
}
for (i in 1:nrow(Stuttgartcentroid_data_onlymetroneighbour)){
  if(Stuttgartcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Stuttgartcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200){
    Stuttgartcentroid_data_onlymetroneighbour$buffer1200m[i] <- 1
  }
}

## Warsaw
Warsawcentroid_data_onlymetroneighbour <- Warsawcentroid_data[which(Warsawcentroid_data$station_openingyear == 2015 & Warsawcentroid_data$NEAR_DIST < 1201),]

Warsawcentroid_data_onlymetroneighbour$buffer400m <- 0
Warsawcentroid_data_onlymetroneighbour$buffer800m <- 0
Warsawcentroid_data_onlymetroneighbour$buffer1200m <- 0

for (i in 1:nrow(Warsawcentroid_data_onlymetroneighbour)){
  if(Warsawcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400){
    Warsawcentroid_data_onlymetroneighbour$buffer400m[i] <- 1
  }
}
for (i in 1:nrow(Warsawcentroid_data_onlymetroneighbour)){
  if(Warsawcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Warsawcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800){
    Warsawcentroid_data_onlymetroneighbour$buffer800m[i] <- 1
  }
}
for (i in 1:nrow(Warsawcentroid_data_onlymetroneighbour)){
  if(Warsawcentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Warsawcentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200){
    Warsawcentroid_data_onlymetroneighbour$buffer1200m[i] <- 1
  }
}


## Sofia
Sofiacentroid_data_onlymetroneighbour <- Sofiacentroid_data[which(Sofiacentroid_data$Openingyear > 2014 & Sofiacentroid_data$NEAR_DIST < 1201),]

Sofiacentroid_data_onlymetroneighbour$buffer400m <- 0
Sofiacentroid_data_onlymetroneighbour$buffer800m <- 0
Sofiacentroid_data_onlymetroneighbour$buffer1200m <- 0

for (i in 1:nrow(Sofiacentroid_data_onlymetroneighbour)){
  if(Sofiacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400){
    Sofiacentroid_data_onlymetroneighbour$buffer400m[i] <- 1
  }
}
for (i in 1:nrow(Sofiacentroid_data_onlymetroneighbour)){
  if(Sofiacentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Sofiacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800){
    Sofiacentroid_data_onlymetroneighbour$buffer800m[i] <- 1
  }
}
for (i in 1:nrow(Sofiacentroid_data_onlymetroneighbour)){
  if(Sofiacentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Sofiacentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200){
    Sofiacentroid_data_onlymetroneighbour$buffer1200m[i] <- 1
  }
}



## Milan
Milancentroid_data_onlymetroneighbour <- Milancentroid_data[which(Milancentroid_data$Openingyear > 2013 & Milancentroid_data$NEAR_DIST < 1201),]
Milancentroid_data_onlymetroneighbour <- Milancentroid_data_onlymetroneighbour[which(!is.na(Milancentroid_data_onlymetroneighbour$streetdistance)),]
Milancentroid_data_onlymetroneighbour$streetdistance <- as.numeric(Milancentroid_data_onlymetroneighbour$streetdistance)

Milancentroid_data_onlymetroneighbour$buffer400m <- 0
Milancentroid_data_onlymetroneighbour$buffer800m <- 0
Milancentroid_data_onlymetroneighbour$buffer1200m <- 0

for (i in 1:nrow(Milancentroid_data_onlymetroneighbour)){
  if(Milancentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 400){
    Milancentroid_data_onlymetroneighbour$buffer400m[i] <- 1
  }
}
for (i in 1:nrow(Milancentroid_data_onlymetroneighbour)){
  if(Milancentroid_data_onlymetroneighbour$NEAR_DIST[i] > 400 & Milancentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 800){
    Milancentroid_data_onlymetroneighbour$buffer800m[i] <- 1
  }
}
for (i in 1:nrow(Milancentroid_data_onlymetroneighbour)){
  if(Milancentroid_data_onlymetroneighbour$NEAR_DIST[i] > 800 & Milancentroid_data_onlymetroneighbour$NEAR_DIST[i] <= 1200){
    Milancentroid_data_onlymetroneighbour$buffer1200m[i] <- 1
  }
}

#########################################################################################
############################## Regression models ########################################
########################################################################################

save_table <- function(lm, tablename){
  x <- as.data.frame(lm$coefficients)
  x$significance <- ''
  x$significance[x$`Pr(>|t|)` < 0.05] <- "*"
  x$significance[x$`Pr(>|t|)` < 0.1] <- "."
  x$significance[x$`Pr(>|t|)` < 0.01] <- "**"
  x$significance[x$`Pr(>|t|)` < 0.001] <- "***"
  x$adjustedrsq <- lm$adj.r.squared
  x$fstatistics_value <- lm$fstatistic[1]
  x$fstatistics_Nrdf <- lm$fstatistic[2]
  x$fstatistics_dendf <- lm$fstatistic[3]
  x$call <- lm$call$formula[1]
  x$type <- as.character(tablename)
  return(x)
}

## Rome
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + entropy_13 +streetdistance , Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + entropy_13 + streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + entropy_13 + streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed+ centroid_distance_center + rowsum_13  + streetdistance,  Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed+ centroid_distance_center  + streetdistance,  Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowsum_13  + streetdistance , Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowsum_13  + streetdistance , Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowsum_13  + streetdistance + metroneighbourhood:rowsum_13 , Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowsum_13  + streetdistance + metroneighbourhood:rowsum_13 , Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowentropy_pre14 +streetdistance , Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +streetdistance , Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowentropy_pre14 + streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowentropy_pre14 + streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Rome_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff, density_diff2, density_PSM, density_ebal, density_interaction, density_interaction2,  density_interaction2,  multifunctionality_diff, multifunctionality_diff2, multifunctionality_PSM, multifunctionality_ebal)
Rome_reg_results$city <- "Rome"


lm <- lm(rowentropy_13 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + streetdistance, Rome_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_13 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

plot(lm, which = 1)

par(mfrow=c(2,2))
plot(lm)


vif.1 <- usdm::vif(Rome_PSM_dataset[,c("metroneighbourhood", "population_2012", "popchange", "perc_unemployed", "centroid_distance_center", "entropy_13", "streetdistance")])
vif.1
vif.1 <- usdm::vif(Rome_PSM_dataset[,c("metroneighbourhood", "population_2012", "popchange", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")])
vif.1
install.packages("corpcor")
library(corpcor)
cor2pcor(cov(Rome_PSM_dataset[,c("metroneighbourhood", "pop_2012_buffer", "popchange_buffer", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")]))
cor2pcor(cov(Romecentroid_data_onlymetroneighbour[!is.na(Romecentroid_data_onlymetroneighbour[,c("buffer400m", "buffer800m","pop_2012_buffer", "popchange_buffer", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")]), c("buffer400m", "buffer800m","pop_2012_buffer", "popchange_buffer", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")]))
cor2pcor(cov(x))

lm <- lm(buffer400m ~ rowsum_13, Romecentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)


x <-Romecentroid_data_onlymetroneighbour[!is.na(Romecentroid_data_onlymetroneighbour[,c("buffer400m", "buffer800m","pop_2012_buffer", "popchange_buffer", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")]), c("buffer400m", "buffer800m","pop_2012_buffer", "popchange_buffer", "perc_unemployed", "centroid_distance_center", "rowsum_13", "streetdistance")]


Rome <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + rowsum_13  + streetdistance, Rome_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
jpeg("Rome Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Rome)
dev.off()

Rome <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + entropy_13  + streetdistance, Rome_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
jpeg("Rome Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Rome)
dev.off()


## Barcelona
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + entropy_pre15 + streetdistance, Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + entropy_pre15 + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + entropy_pre15 + streetdistance,Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase2 ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance , Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase2 ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance , Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(density_increase ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance+ metroneighbourhood:rowsum15 , Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance + metroneighbourhood:rowsum15, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + rowentropy_pre16 + streetdistance, Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center  + streetdistance, Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + rowentropy_pre16 + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + rowentropy_pre16 + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")


Barce_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff, density_diff2, density_PSM, density_ebal, density_interaction, density_interaction2,  density_interaction2,  multifunctionality_diff, multifunctionality_diff2, multifunctionality_PSM, multifunctionality_ebal)
Barce_reg_results$city <- "Barcelona"


lm <- lm(rowentropy_16 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + streetdistance, Barce_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_16 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum15, Barcecentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)

Barce <- lm(density_increase ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + rowsum15 + streetdistance, Barce_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")
jpeg("Barce Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Barce)
dev.off()

Barce <-lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed_2015 + centroid_distance_center + entropy_pre16 + streetdistance, Barce_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")
jpeg("Barce Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Barce)
dev.off()


## Helsinki
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + entropy_16 + streetdistance, Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + entropy_16 + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + entropy_16 + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance + metroneighbourhood:rowsum_16, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance + metroneighbourhood:rowsum_16, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowentropy_pre17 + streetdistance, Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + streetdistance, Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowentropy_pre17  + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowentropy_pre17  + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")


Helsinki_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff, density_diff2, density_PSM, density_interaction, density_interaction2,  multifunctionality_diff, multifunctionality_diff2, multifunctionality_PSM, density_ebal, multifunctionality_ebal)
Helsinki_reg_results$city <- "Helsinki"

write.csv(Helsinki_reg_results, "Helsinki_reg_results.csv")

lm <- lm(rowentropy_17 ~ metroneighbourhood + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center  + streetdistance, Helsinki_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_17 ~ metroneighbourhood + pop_2012_buffer  + popchange_buffer + unemployment_perc2012 + centroid_dist_center  + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum_16, Helsinkicentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)



Helsinki <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + rowsum_16 + streetdistance, Helsinki_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")
jpeg("Helsinki Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Helsinki)
dev.off()

Helsinki <-lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + entropy_16 + streetdistance, Helsinki_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")
jpeg("Helsinki Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Helsinki)
dev.off()


## Vienna
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + entropy_16 + streetdistance, Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + entropy_16 + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + entropy_16 + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance , Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance + metroneighbourhood:rowsum_16, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance + metroneighbourhood:rowsum_16, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowentropy_pre17 + streetdistance, Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center  + streetdistance, Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowentropy_pre17 + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance , Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowentropy_pre17 + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Vienna_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff, density_PSM, density_diff2, density_interaction, density_interaction2,  multifunctionality_diff2, multifunctionality_diff, multifunctionality_PSM, density_ebal, multifunctionality_ebal)
Vienna_reg_results$city <- "Vienna"


lm <- lm(rowentropy_16 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + streetdistance, Vienna_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_16 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))


lm <- lm(buffer400m ~ rowsum_16, Viennacentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)


Vienna <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + rowsum_16 + streetdistance , Vienna_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")
jpeg("Vienna Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Vienna)
dev.off()

Vienna <-lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + entropy_16 + streetdistance, Vienna_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")
jpeg("Vienna Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Vienna)
dev.off()


## Budapest
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + entropy_pre14 + streetdistance , Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + entropy_pre14 + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + entropy_pre14 + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance,  Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + streetdistance,  Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(density_increase  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance + metroneighbourhood:rowsum14, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance + metroneighbourhood:rowsum14, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowentropy_pre14 + streetdistance , Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + streetdistance , Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowentropy_pre14 + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(multifunctionality_diff  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowentropy_pre14 + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Budapest_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff2, density_interaction, density_interaction2,  multifunctionality_diff2, density_diff, density_PSM, multifunctionality_diff, multifunctionality_PSM, density_ebal, multifunctionality_ebal)
Budapest_reg_results$city <- "Budapest"

lm <- lm(rowentropy_14  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + streetdistance, Budapest_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_14  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))


lm <- lm(buffer400m ~ rowsum14, Budapestcentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)


Budapest <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + rowsum14 + streetdistance, Budapest_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")
jpeg("Budapest Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Budapest)
dev.off()

Budapest <-lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + entropy_pre14 + streetdistance , Budapest_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")
jpeg("Budapest Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Budapest)
dev.off()



## Stuttgart
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer +  popchange_buffer + perc_unemployed12 + centroid_dist_center + entropy_15 + streetdistance, Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + entropy_15 + centroid_dist_center+ streetdistance , Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + entropy_15 + centroid_dist_center+ streetdistance , Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer +  popchange_buffer + perc_unemployed12 + centroid_dist_center + rowsum_15 + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer +  popchange_buffer + perc_unemployed12 + centroid_dist_center  + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center  + rowsum_15 + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center  + rowsum_15 + streetdistance + metroneighbourhood:rowsum_15, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center  + rowsum_15 + streetdistance + metroneighbourhood:rowsum_15,  Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + rowentropy_pre16 + streetdistance, Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center  + streetdistance, Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + rowentropy_pre16 + streetdistance , Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center  + rowsum_15 + streetdistance,  Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + rowentropy_pre16 + streetdistance , Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Stuttgart_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff2, density_interaction, density_interaction2,  multifunctionality_diff2, density_diff, density_PSM, multifunctionality_diff, multifunctionality_PSM, density_ebal, multifunctionality_ebal)
Stuttgart_reg_results$city <- "Stuttgart"


lm <- lm(rowentropy_15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + streetdistance , Stuttgart_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum_15, Stuttgartcentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)


Stuttgart <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + rowsum_15 + streetdistance, Stuttgart_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")
jpeg("Stuttgart Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Stuttgart)
dev.off()

Stuttgart <-lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + entropy_15 + streetdistance, Stuttgart_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")
jpeg("Stuttgart Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Stuttgart)
dev.off()

## Warsaw
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + entropy_14 + streetdistance, Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + entropy_14 + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + entropy_14 + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_PSM")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14 + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14+ streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14+ streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowentropy_pre15 + streetdistance, Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center  + streetdistance, Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowentropy_pre15 + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowentropy_pre15 + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Warsaw_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff2, density_interaction, density_interaction2,  multifunctionality_diff2, density_diff, density_PSM, density_ebal, multifunctionality_diff, multifunctionality_PSM, multifunctionality_ebal)
Warsaw_reg_results$city <- "Warsaw"

lm <- lm(rowentropy_pre15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + streetdistance, Warsaw_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_pre15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum_14, Warsawcentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)


Warsaw <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + rowsum_14 + streetdistance, Warsaw_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")
jpeg("Warsaw Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Warsaw)
dev.off()

Warsaw <-lm(multifunctionality_increase ~  metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + entropy_14 + streetdistance , Warsaw_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")
jpeg("Warsaw Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Warsaw)
dev.off()



## Sofia
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + entropy_14 + streetdistance, Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + entropy_14 + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + entropy_14 + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowsum_14 + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center + rowsum_14+ streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center + rowsum_14+ streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Sofia_PSM_dataset)
lm <-  beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowentropy_pre15 + streetdistance, Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center  + streetdistance, Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowentropy_pre15 + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowentropy_pre15 + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Sofia_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff2, density_interaction, density_interaction2,  multifunctionality_diff2, density_diff, density_PSM, density_ebal, multifunctionality_diff, multifunctionality_PSM, multifunctionality_ebal)
Sofia_reg_results$city <- "Sofia"


lm <- lm(rowentropy_pre15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + streetdistance, Sofia_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_pre15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum_14, Sofiacentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)

lm <- lm(centroid_dist_cente ~ streetdistance, Sofia_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE)

Sofia <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowsum_14 + streetdistance, Sofia_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")
jpeg("Sofia Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Sofia)
dev.off()

Sofia <-lm(entropy_diff ~  metroneighbourhood + pop_2012   + centroid_dist_center + entropy_14 + streetdistance , Sofia_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")
jpeg("Sofia Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Sofia)
dev.off()

## Milan
#$$$$$$$$$$$$$#
setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")

lm <- lm(entropy_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011 +centroid_dist_center + entropy_14 + streetdistance, Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
entropy_diff <- save_table(lm, "entropy_diff")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + perc_unemployed2011 + centroid_dist_center + entropy_14 + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_PSM <- save_table(lm, "entropy_PSM")

lm <- lm(entropy_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + perc_unemployed2011 + centroid_dist_center + entropy_14 + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
entropy_ebal <- save_table(lm, "entropy_ebal")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011+ centroid_dist_center + rowsum_14 + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff <- save_table(lm, "density_diff")

lm <- lm(density_increase ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011+ centroid_dist_center + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
density_diff2 <- save_table(lm, "density_diff2")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowsum_14+ streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_PSM <- save_table(lm, "density_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction <- save_table(lm, "density_interaction")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_interaction2 <- save_table(lm, "density_interaction2")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011+ centroid_dist_center + rowentropy_pre15 + streetdistance, Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff <- save_table(lm, "multifunctionality_diff")

lm <- lm(multifunctionality_diff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011+ centroid_dist_center + streetdistance, Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
multifunctionality_diff2 <- save_table(lm, "multifunctionality_diff2")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowentropy_pre15 + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_PSM <- save_table(lm, "multifunctionality_PSM")

lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowsum_14+ streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
density_ebal <- save_table(lm, "density_ebal")

lm <- lm(multifunctionality_diff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowentropy_pre15 + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
multifunctionality_ebal <- save_table(lm, "multifunctionality_ebal")

Milan_reg_results <- rbind(entropy_diff, entropy_PSM, entropy_ebal, density_diff2, density_interaction, density_interaction2,  density_interaction2, multifunctionality_diff2, density_diff, density_PSM, multifunctionality_diff, multifunctionality_PSM, density_ebal, multifunctionality_ebal)
Milan_reg_results$city <- "Milan"


lm <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + rowsum_14+ streetdistance + metroneighbourhood:rowsum_14, Milan_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(buffer400m ~ rowsum_14, Milancentroid_data_onlymetroneighbour)
beta(model = lm,  x = TRUE, y = TRUE)

lm <- lm(rowentropy_pre15 ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center  + streetdistance, Milan_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(rowentropy_pre15  ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center  + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))


Milan <- lm(density_increase ~ metroneighbourhood + pop_2012_buffer + popchange_buffer  + centroid_dist_center + rowsum_14 + streetdistance, Milan_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")
jpeg("Milan Density Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Milan)
dev.off()

Milan <-lm(multifunctionality_increase ~  metroneighbourhood + pop_2012   + centroid_dist_center + rowentropy_14 + streetdistance , Milan_PSM_dataset)
setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")
jpeg("Milan Multifunctionality Model Assumptions Test.jpeg",  width=7, height=7, units = 'in', res = 1500)
autoplot(Milan)
dev.off()

###################### Comparison

allRegression_results <- rbind(Rome_reg_results,
                               Barce_reg_results,
                               Vienna_reg_results,
                               Budapest_reg_results,
                               Helsinki_reg_results,
                               Stuttgart_reg_results,
                               Sofia_reg_results,
                               Warsaw_reg_results,
                               Milan_reg_results)

Milan_reg_results
Rome_reg_results
Barce_reg_results
Helsinki_reg_results
Vienna_reg_results
Sofia_reg_results
Budapest_reg_results
Warsaw_reg_results
Stuttgart_reg_results

setwd("C:/Dokumente/Utrecht/Master_these/results")
density_diff2 <- allRegression_results[allRegression_results$type == "density_diff2", ]
write.csv(density_diff2, "density_diff2.csv")

density_interaction <- allRegression_results[allRegression_results$type == "density_interaction", ]
write.csv(density_interaction, "density_interaction.csv")

density_interaction2 <- allRegression_results[allRegression_results$type == "density_interaction2", ]
write.csv(density_interaction2, "density_interaction2.csv")

multifunctionality_diff2 <- allRegression_results[allRegression_results$type == "multifunctionality_diff2", ]
write.csv(multifunctionality_diff2, "multifunctionality_diff2.csv")

allRegression_results <- allRegression_results[order(allRegression_results$type), ]
allRegression_results <- allRegression_results[allRegression_results$type != "entropy_diff", ]
allRegression_results <- allRegression_results[allRegression_results$type != "entropy_PSM", ]
allRegression_results <- allRegression_results[allRegression_results$type != "entropy_ebal", ]

write.csv(allRegression_results, "allRegression_results.csv")



####### Testing venue type heterogeneity of number of openings ############################

Heterogeneity_output <- function(table){
  output <- matrix(nrow = 10, ncol = 20)
  output <- as.data.frame(output)
  rownames(output) <-  c("(Intercept)",
                         "Metroneighbourhood", "Population 2012", "Pop change 2012-2018", 
                         "Percentage Unemployed", "Distance to Center", "Nr of venue type before station", 
                         "Street-level Distance to nearest Subcenter", "", "Adjusted R-squared:")
  colnames(output) <- c("Food_ebal_coeff","Food_ebal_p", "Food_PSM_coeff", "Food_PSM_p", "Proff_ebal_coeff", 
                        "Proff_ebal_p", "Proff_PSM_coeff", "Proff_PSM_p", "Arts_ebal_coeff", 
                        "Arts_ebal_p", "Arts_PSM_coeff", "Arts_PSM_p", "Shops_ebal_coeff", "Shops_ebal_p", "Shops_PSM_coeff",  
                        "Shops_PSM_p", "Nightlife_ebal_coeff", "Nightlife_ebal_p", "Nightlife_PSM_coeff", "Nightlife_PSM_p")
  
  output[1:8, 1] <- table$Estimate[table$type == "Food_ebal"]
  output[1:8, 2] <- table$significance[table$type == "Food_ebal"]
  output[10, 1] <- table$adjustedrsq[table$type == "Food_ebal"][1]
  
  output[1:8, 3] <- table$Estimate[table$type == "Food_PSM"]
  output[1:8, 4] <- table$significance[table$type == "Food_PSM"]
  output[10, 3] <- table$adjustedrsq[table$type == "Food_PSM"][1]
  
  output[1:8, 5] <- table$Estimate[table$type == "Proff_ebal"]
  output[1:8, 6] <- table$significance[table$type == "Proff_ebal"]
  output[10, 5] <- table$adjustedrsq[table$type == "Proff_ebal"][1]
  
  output[1:8, 7] <- table$Estimate[table$type == "Proff_PSM"]
  output[1:8, 8] <- table$significance[table$type == "Proff_PSM"]
  output[10, 7] <- table$adjustedrsq[table$type == "Proff_PSM"][1]
  
  output[1:8, 9] <- table$Estimate[table$type == "Arts_ebal"]
  output[1:8, 10] <- table$significance[table$type == "Arts_ebal"]
  output[10, 9] <- table$adjustedrsq[table$type == "Arts_ebal"][1]
  
  output[1:8, 11] <- table$Estimate[table$type == "Arts_PSM"]
  output[1:8, 12] <- table$significance[table$type == "Arts_PSM"]
  output[10, 11] <- table$adjustedrsq[table$type == "Arts_PSM"][1]
  
  output[1:8, 13] <- table$Estimate[table$type == "Shops_ebal"]
  output[1:8, 14] <- table$significance[table$type == "Shops_ebal"]
  output[10, 13] <- table$adjustedrsq[table$type == "Shops_ebal"][1]
  
  output[1:8, 15] <- table$Estimate[table$type == "Shops_PSM"]
  output[1:8, 16] <- table$significance[table$type == "Shops_PSM"]
  output[10, 15] <- table$adjustedrsq[table$type == "Shops_PSM"][1]
  
  output[1:8, 17] <- table$Estimate[table$type == "Nightlife_ebal"]
  output[1:8, 18] <- table$significance[table$type == "Nightlife_ebal"]
  output[10, 17] <- table$adjustedrsq[table$type == "Nightlife_ebal"][1]
  
  output[1:8, 19] <- table$Estimate[table$type == "Nightlife_PSM"]
  output[1:8, 20] <- table$significance[table$type == "Nightlife_PSM"]
  output[10, 19] <- table$adjustedrsq[table$type == "Nightlife_PSM"][1]
  return(output)
}

setwd("C:/Dokumente/Utrecht/Master_these/results/Heterogeneity")

## Rome
#$$$$$$$$$$$$$#
c("nr_ArtsEnter_1316_300mbuff", "nr_Nightlife_1316_300mbuff", "nr_Proff_1316_300mbuff" , 
  "nr_Outdoors_1316_300mbuff", "nr_ShopsServ_1316_300mbuff" , "nr_Food_1316_300mbuff" ,
  "nr_Proff_13_300mbuff", "nr_Outdoors_13_300mbuff",  "nr_ShopsServ_13_300mbuff", 
  "nr_ArtsEnter_13_300mbuff",  "nr_Food_13_300mbuff", "nr_Nightlife_13_300mbuff")

lm <- lm(nr_Food_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + nr_Food_13_300mbuff+ streetdistance, Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Food_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + nr_Proff_13_300mbuff + streetdistance, Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Proff_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + nr_ArtsEnter_13_300mbuff+ streetdistance, Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + nr_ShopsServ_13_300mbuff+ streetdistance, Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ShopsServ_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center + nr_Nightlife_13_300mbuff+ streetdistance, Romecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Nightlife_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")
lm <- lm(nr_Food_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Food_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Proff_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ShopsServ_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Nightlife_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")
lm <- lm(nr_Outdoors_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Outdoors_13_300mbuff+ streetdistance, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")
lm <- lm(nr_Outdoors_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Outdoors_13_300mbuff+ streetdistance, Rome_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")

Rome_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)

lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Rome_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Nightlife_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Rome_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + centroid_distance_center +  nr_Nightlife_13_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Romecentroid_data_metroandcontrol, weights = Romecentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

#Rome_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Rome_hetero_output <- Heterogeneity_output(Rome_venueheterogeneity)
write.csv(Rome_hetero_output, "Rome_hetero_output.csv")
Rome_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Rome_venueheterogeneity , "Rome_outdoors.csv")

## Barcelona
#$$$$$$$$$$$$$#
c("nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff")

Barcecentroid_data_metroandcontrol$nr_ArtsEnter_1518_300mbuff<- Barcecentroid_data$nr_ArtsEnter_1518_300mbuff[Barcecentroid_data_metroandcontrol$ORIG_FID]
Barcecentroid_data_metroandcontrol$nr_ArtsEnter_15_300mbuff<- Barcecentroid_data$nr_ArtsEnter_15_300mbuff[Barcecentroid_data_metroandcontrol$ORIG_FID]
Barcecentroid_data_metroandcontrol$nr_Nightlife_1518_300mbuff<- Barcecentroid_data$nr_Nightlife_1518_300mbuff[Barcecentroid_data_metroandcontrol$ORIG_FID]
Barcecentroid_data_metroandcontrol$nr_Nightlife_15_300mbuff<- Barcecentroid_data$nr_Nightlife_15_300mbuff[Barcecentroid_data_metroandcontrol$ORIG_FID]
Barce_PSM_dataset$nr_ArtsEnter_1518_300mbuff<- Barcecentroid_data$nr_ArtsEnter_1518_300mbuff[Barce_PSM_dataset$ORIG_FID]
Barce_PSM_dataset$nr_ArtsEnter_15_300mbuff<- Barcecentroid_data$nr_ArtsEnter_15_300mbuff[Barce_PSM_dataset$ORIG_FID]
Barce_PSM_dataset$nr_Nightlife_1518_300mbuff<- Barcecentroid_data$nr_Nightlife_1518_300mbuff[Barce_PSM_dataset$ORIG_FID]
Barce_PSM_dataset$nr_Nightlife_15_300mbuff<- Barcecentroid_data$nr_Nightlife_15_300mbuff[Barce_PSM_dataset$ORIG_FID]


lm <- lm(nr_Food_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Food_15_300mbuff + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Food_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Proff_15_300mbuff + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Proff_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ArtsEnter_15_300mbuff + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ArtsEnter_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ShopsServ_15_300mbuff + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ShopsServ_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Nightlife_15_300mbuff + streetdistance,  Barcecentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Nightlife_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")

lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Food_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Proff_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ShopsServ_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_ArtsEnter_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Nightlife_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Outdoors_15_300mbuff + streetdistance, Barce_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")

lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed_2015+ centroid_distance_center + nr_Outdoors_15_300mbuff + streetdistance, Barcecentroid_data_metroandcontrol, weights = Barcecentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")



Barce_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)
#Barce_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Shops_diff, Shops_PSM)
Barce_hetero_output <- Heterogeneity_output(Barce_venueheterogeneity)
write.csv(Barce_hetero_output, "Barce_hetero_output.csv")

write.csv(Barce_venueheterogeneity, "Barce_venueheterogeneity.csv")

Barce_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM) 
write.csv(Barce_venueheterogeneity, "Barce_outdoors.csv")

## Vienna
#$$$$$$$$$$$$$#
c("nr_Food_16_300mbuff", "nr_ShopsServ_16_300mbuff", "nr_ArtsEnter_16_300mbuff", 
  "nr_Nightlife_16_300mbuff", "nr_Proff_16_300mbuff", "nr_Outdoors_16_300mbuff",
  "nr_Food_1619_300mbuff", "nr_ShopsServ_1619_300mbuff", "nr_ArtsEnter_1619_300mbuff",
  "nr_Nightlife_1619_300mbuff", "nr_Proff_1619_300mbuff", "nr_Outdoors_1619_300mbuff")

lm <- lm(nr_Food_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Food_16_300mbuff + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Food_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Proff_16_300mbuff + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Proff_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ArtsEnter_16_300mbuff + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ArtsEnter_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ShopsServ_16_300mbuff + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ShopsServ_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Nightlife_16_300mbuff + streetdistance,  Viennacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Nightlife_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")

lm <- lm(nr_Food_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Food_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Proff_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ArtsEnter_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_ShopsServ_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Nightlife_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")


lm <- lm(nr_Outdoors_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Outdoors_16_300mbuff + streetdistance, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")
lm <- lm(nr_Outdoors_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center + nr_Outdoors_16_300mbuff + streetdistance, Vienna_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")


lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center +  nr_ArtsEnter_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Vienna_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center +  nr_ArtsEnter_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center +  nr_Nightlife_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Vienna_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed + distance_center +  nr_Nightlife_16_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Viennacentroid_data_metroandcontrol, weights = Viennacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))


Vienna_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)
#Vienna_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Vienna_hetero_output <- Heterogeneity_output(Vienna_venueheterogeneity)
write.csv(Vienna_hetero_output, "Vienna_hetero_output.csv")

Vienna_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Vienna_venueheterogeneity, "Vienna_outdoors.csv")

## Budapest
#$$$$$$$$$$$$$#
c("nr_ArtsEnter_1316_300mbuff", "nr_Nightlife_1316_300mbuff", "nr_Proff_1316_300mbuff" , 
  "nr_Outdoors_1316_300mbuff", "nr_ShopsServ_1316_300mbuff" , "nr_Food_1316_300mbuff" ,
  "nr_Proff_13_300mbuff", "nr_Outdoors_13_300mbuff",  "nr_ShopsServ_13_300mbuff", 
  "nr_ArtsEnter_13_300mbuff",  "nr_Food_13_300mbuff", "nr_Nightlife_13_300mbuff")

lm <- lm(nr_Food_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + nr_Food_13_300mbuff + streetdistance, Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Food_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + nr_Proff_13_300mbuff + streetdistance, Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Proff_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + nr_ArtsEnter_13_300mbuff + streetdistance, Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ArtsEnter_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + nr_ShopsServ_13_300mbuff + streetdistance, Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ShopsServ_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1316_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center + nr_Nightlife_13_300mbuff + streetdistance, Budapestcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Nightlife_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")

lm <- lm(nr_Food_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Food_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Proff_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ArtsEnter_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ShopsServ_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Nightlife_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Outdoors_13_300mbuff + streetdistance, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")
lm <- lm(nr_Outdoors_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Outdoors_13_300mbuff + streetdistance, Budapest_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")


lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Budapest_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_ArtsEnter_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Nightlife_13_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Budapest_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1316_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + percentage_unemployed + distance_center +  nr_Nightlife_13_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Budapestcentroid_data_metroandcontrol, weights = Budapestcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Budapest_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)

#Budapest_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Budapest_hetero_output <- Heterogeneity_output(Budapest_venueheterogeneity)
write.csv(Budapest_hetero_output, "Budapest_hetero_output.csv")

Budapest_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Budapest_venueheterogeneity, "Budapest_Outdoors.csv")

## Stuttgart
#$$$$$$$$$$$$$#
c("nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff")

lm <- lm(nr_Food_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance , Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance,  Stuttgartcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")

lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance , Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")
lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Stuttgart_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")


lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Stuttgart_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Stuttgart_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed12 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Stuttgartcentroid_data_metroandcontrol, weights = Stuttgartcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Stuttgart_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)

#Stuttgart_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Stuttgart_hetero_output <- Heterogeneity_output(Stuttgart_venueheterogeneity)
write.csv(Stuttgart_hetero_output, "Stuttgart_hetero_output.csv")

Stuttgart_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM) 
write.csv(Stuttgart_venueheterogeneity, "Stuttgart_Outdoors.csv")


## Helsinki
#$$$$$$$$$$$$$#
c("nr_Food_16_300mbuff", "nr_ShopsServ_16_300mbuff", "nr_ArtsEnter_16_300mbuff", 
  "nr_Nightlife_16_300mbuff", "nr_Proff_16_300mbuff", "nr_Outdoors_16_300mbuff",
  "nr_Food_1619_300mbuff", "nr_ShopsServ_1619_300mbuff", "nr_ArtsEnter_1619_300mbuff",
  "nr_Nightlife_1619_300mbuff", "nr_Proff_1619_300mbuff", "nr_Outdoors_1619_300mbuff")

lm <- lm(nr_Food_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Food_16_300mbuff + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Food_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Proff_16_300mbuff + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Proff_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ArtsEnter_16_300mbuff + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ArtsEnter_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ShopsServ_16_300mbuff + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ShopsServ_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1619_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Nightlife_16_300mbuff + streetdistance,  Helsinkicentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Nightlife_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")
lm <- lm(nr_Food_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Food_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Proff_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ArtsEnter_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_ShopsServ_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Nightlife_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Outdoors_16_300mbuff + streetdistance, Helsinki_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")
lm <- lm(nr_Outdoors_1619_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center + nr_Outdoors_16_300mbuff + streetdistance, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")

lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center +  nr_ArtsEnter_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Helsinki_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center +  nr_ArtsEnter_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center +  nr_Nightlife_16_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Helsinki_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1619_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployment_perc2012 + centroid_dist_center +  nr_Nightlife_16_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Helsinkicentroid_data_metroandcontrol, weights = Helsinkicentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Helsinki_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)
#Helsinki_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Helsinki_hetero_output <- Heterogeneity_output(Helsinki_venueheterogeneity)
write.csv(Helsinki_hetero_output, "Helsinki_hetero_output.csv")

Helsinki_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Helsinki_venueheterogeneity, "Helsinki_Outdoors.csv")


## Warsaw
#$$$$$$$$$$$$$#
c("nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff")

lm <- lm(nr_Food_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance,  Warsawcentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Warsaw_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")
lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")


lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Warsaw_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Warsaw_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + unemployed_2010 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Warsawcentroid_data_metroandcontrol, weights = Warsawcentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Warsaw_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)
#Warsaw_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Warsaw_hetero_output <- Heterogeneity_output(Warsaw_venueheterogeneity)
write.csv(Warsaw_hetero_output, "Warsaw_hetero_output.csv")

Warsaw_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Warsaw_venueheterogeneity, "Warsaw_outdoors.csv")

## Sofia
#$$$$$$$$$$$$$#
c("nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff")

lm <- lm(nr_Food_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_Food_15_300mbuff + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_Proff_15_300mbuff + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance,  Sofiacentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer  + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")


lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Sofia_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")
lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")

lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Sofia_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Sofia_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Sofiacentroid_data_metroandcontrol, weights = Sofiacentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Heterogeneity_output3 <- function(table){
  output <- matrix(nrow = 10, ncol = 20)
  output <- as.data.frame(output)
  rownames(output) <-  c("(Intercept)", 
                         "Metroneighbourhood", "Population 2012", "Pop change 2012-2018", 
                         "Percentage Unemployed", "Distance to Center", "Nr of venue type before station", 
                         "Street-level Distance to nearest Subcenter", "", "Adjusted R-squared:")
  colnames(output) <- c("Food_ebal_coeff","Food_ebal_p", "Food_PSM_coeff", "Food_PSM_p", "Proff_ebal_coeff", 
                        "Proff_ebal_p", "Proff_PSM_coeff", "Proff_PSM_p", "Arts_ebal_coeff", 
                        "Arts_ebal_p", "Arts_PSM_coeff", "Arts_PSM_p", "Shops_ebal_coeff", "Shops_ebal_p", "Shops_PSM_coeff",  
                        "Shops_PSM_p", "Nightlife_ebal_coeff", "Nightlife_ebal_p", "Nightlife_PSM_coeff", "Nightlife_PSM_p")
  
  output[c(1,2,3,5,6,7,8), 1] <- table$Estimate[table$type == "Food_ebal"]
  output[c(1,2,3,5,6,7,8), 2] <- table$significance[table$type == "Food_ebal"]
  output[10, 1] <- table$adjustedrsq[table$type == "Food_ebal"][1]
  
  output[c(1,2,3,5,6,7,8), 3] <- table$Estimate[table$type == "Food_PSM"]
  output[c(1,2,3,5,6,7,8), 4] <- table$significance[table$type == "Food_PSM"]
  output[10, 3] <- table$adjustedrsq[table$type == "Food_PSM"][1]
  
  output[c(1,2,3,5,6,7,8), 5] <- table$Estimate[table$type == "Proff_ebal"]
  output[c(1,2,3,5,6,7,8), 6] <- table$significance[table$type == "Proff_ebal"]
  output[10, 5] <- table$adjustedrsq[table$type == "Proff_ebal"][1]
  
  output[c(1,2,3,5,6,7,8), 7] <- table$Estimate[table$type == "Proff_PSM"]
  output[c(1,2,3,5,6,7,8), 8] <- table$significance[table$type == "Proff_PSM"]
  output[10, 7] <- table$adjustedrsq[table$type == "Proff_PSM"][1]
  
  output[c(1,2,3,5,6,7,8), 9] <- table$Estimate[table$type == "Arts_ebal"]
  output[c(1,2,3,5,6,7,8), 10] <- table$significance[table$type == "Arts_ebal"]
  output[10, 9] <- table$adjustedrsq[table$type == "Arts_ebal"][1]
  
  output[c(1,2,3,5,6,7,8), 11] <- table$Estimate[table$type == "Arts_PSM"]
  output[c(1,2,3,5,6,7,8), 12] <- table$significance[table$type == "Arts_PSM"]
  output[10, 11] <- table$adjustedrsq[table$type == "Arts_PSM"][1]
  
  output[c(1,2,3,5,6,7,8), 13] <- table$Estimate[table$type == "Shops_ebal"]
  output[c(1,2,3,5,6,7,8), 14] <- table$significance[table$type == "Shops_ebal"]
  output[10, 13] <- table$adjustedrsq[table$type == "Shops_ebal"][1]
  
  output[c(1,2,3,5,6,7,8), 15] <- table$Estimate[table$type == "Shops_PSM"]
  output[c(1,2,3,5,6,7,8), 16] <- table$significance[table$type == "Shops_PSM"]
  output[10, 15] <- table$adjustedrsq[table$type == "Shops_PSM"][1]
  
  output[c(1,2,3,5,6,7,8), 17] <- table$Estimate[table$type == "Nightlife_ebal"]
  output[c(1,2,3,5,6,7,8), 18] <- table$significance[table$type == "Nightlife_ebal"]
  output[10, 17] <- table$adjustedrsq[table$type == "Nightlife_ebal"][1]
  
  output[c(1,2,3,5,6,7,8), 19] <- table$Estimate[table$type == "Nightlife_PSM"]
  output[c(1,2,3,5,6,7,8), 20] <- table$significance[table$type == "Nightlife_PSM"]
  output[10, 19] <- table$adjustedrsq[table$type == "Nightlife_PSM"][1]
  return(output)
}



Sofia_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)

#Sofia_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Sofia_hetero_output <- Heterogeneity_output3(Sofia_venueheterogeneity)
write.csv(Sofia_hetero_output, "Sofia_hetero_output.csv")

Sofia_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Sofia_venueheterogeneity , "Sofia_outdoors.csv")


## Milan
#$$$$$$$$$$$$$#
c("nr_Food_15_300mbuff", "nr_ShopsServ_15_300mbuff", "nr_ArtsEnter_15_300mbuff", 
  "nr_Nightlife_15_300mbuff", "nr_Proff_15_300mbuff", "nr_Outdoors_15_300mbuff",
  "nr_Food_1518_300mbuff", "nr_ShopsServ_1518_300mbuff", "nr_ArtsEnter_1518_300mbuff",
  "nr_Nightlife_1518_300mbuff", "nr_Proff_1518_300mbuff", "nr_Outdoors_1518_300mbuff")

lm <- lm(nr_Food_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer  + perc_unemployed2011+ centroid_dist_center + nr_Food_15_300mbuff + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Food_diff <- save_table(lm, tablename = "Food_diff")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011+ centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_PSM <- save_table(lm, tablename = "Food_PSM")
lm <- lm(nr_Proff_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed2011+ centroid_dist_center + nr_Proff_15_300mbuff + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Proff_diff <- save_table(lm, tablename = "Proff_diff")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011+ centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_PSM <- save_table(lm, tablename = "Proff_PSM")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Arts_diff <- save_table(lm, tablename = "Arts_diff")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_PSM <- save_table(lm, tablename = "Arts_PSM")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Shops_diff <- save_table(lm, tablename = "Shops_diff")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_PSM <- save_table(lm, tablename = "Shops_PSM")
lm <- lm(nr_Nightlife_1518_300mbuff ~ buffer400m + buffer800m + pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance,  Milancentroid_data_onlymetroneighbour)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("buffer400m", "buffer800m"))
Nightlife_diff <- save_table(lm, tablename = "Nightlife_diff")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_PSM <- save_table(lm, tablename = "Nightlife_PSM")
lm <- lm(nr_Food_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011+ centroid_dist_center + nr_Food_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Food_ebal <- save_table(lm, tablename = "Food_ebal")
lm <- lm(nr_Proff_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011+ centroid_dist_center + nr_Proff_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Proff_ebal <- save_table(lm, tablename = "Proff_ebal")
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Arts_ebal <- save_table(lm, tablename = "Arts_ebal")
lm <- lm(nr_ShopsServ_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center + nr_ShopsServ_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Shops_ebal <- save_table(lm, tablename = "Shops_ebal")
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_Nightlife_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Nightlife_ebal <- save_table(lm, tablename = "Nightlife_ebal")

lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Milan_PSM_dataset)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_PSM <- save_table(lm, tablename = "Outdoors_PSM")
lm <- lm(nr_Outdoors_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_Outdoors_15_300mbuff + streetdistance, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
lm <- beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
Outdoors_ebal <- save_table(lm, tablename = "Outdoors_ebal")


lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Milan_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center +  nr_ArtsEnter_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance + metroneighbourhood:popchange_buffer, Milan_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))
lm <- lm(nr_Nightlife_1518_300mbuff ~ metroneighbourhood + pop_2012_buffer + popchange_buffer + perc_unemployed2011 + centroid_dist_center +  nr_Nightlife_15_300mbuff+ streetdistance +metroneighbourhood:popchange_buffer, Milancentroid_data_metroandcontrol, weights = Milancentroid_data_metroandcontrol$eb_weights)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

Milan_venueheterogeneity <- rbind(Food_ebal, Food_PSM, Proff_ebal, Proff_PSM, Arts_ebal, Arts_PSM, Shops_ebal, Shops_PSM, Nightlife_ebal, Nightlife_PSM)

#Milan_venueheterogeneity <- rbind(Food_diff, Food_PSM, Proff_diff, Proff_PSM, Arts_diff, Arts_PSM, Shops_diff, Shops_PSM, Nightlife_diff, Nightlife_PSM)
Milan_hetero_output <- Heterogeneity_output(Milan_venueheterogeneity)
write.csv(Milan_hetero_output, "Milan_hetero_output.csv")

Milan_venueheterogeneity <- rbind(Outdoors_ebal, Outdoors_PSM)
write.csv(Milan_venueheterogeneity, "Milan_outdoors.csv")

lm <- lm(nr_ArtsEnter_1518_300mbuff ~ metroneighbourhood +  pop_2012_buffer + popchange_buffer +  perc_unemployed2011 + centroid_dist_center + nr_ArtsEnter_15_300mbuff + streetdistance + metroneighbourhood:popchange_buffer, Milan_PSM_dataset)
beta(model = lm,  x = TRUE, y = TRUE, skip = c("metroneighbourhood"))

setwd("C:/Dokumente/Utrecht/Master_these/results/Heterogeneity")


Venueheterogeneity_all <- rbind(
  Rome_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Warsaw_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Sofia_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Vienna_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Helsinki_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Stuttgart_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Budapest_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ],
  Milan_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ]
)

Venueheterogeneity_Barce <- Barce_hetero_output[c("Metroneighbourhood", "Nr of venue type before station", 9, "Adjusted R-squared:"), ]
write.csv(Venueheterogeneity_all, "Venueheterogeneity_all.csv")
write.csv(Venueheterogeneity_Barce, "Venueheterogeneity_Barce.csv")


Venueheterogeneity_complete <- rbind(
  Rome_hetero_output,
  Warsaw_hetero_output,
  Sofia_hetero_output,
  Vienna_hetero_output,
  Helsinki_hetero_output,
  Stuttgart_hetero_output,
  Budapest_hetero_output,
  Milan_hetero_output
)

Venueheterogeneity_complete[81:90, c(1:8,13:16)] <- Barce_hetero_output


write.csv(Venueheterogeneity_complete, "Venueheterogeneity_complete.csv")

Barce_hetero_output


rownames(Milan_hetero_output)

##Coefficients and Significance comparison
coeff_compare <- as.data.frame(matrix(nrow = 9, ncol = 10))
colnames(coeff_compare) <- c("Food_PSM_coeff", "Food_PSM_p",  "Proff_PSM_coeff", "Proff_PSM_p",  "Arts_PSM_coeff", "Arts_PSM_p", 
                             "Shops_PSM_coeff", "Shops_PSM_p",  "Nightlife_PSM_coeff", "Nightlife_PSM_p")
rownames(coeff_compare) <- c("Rome", "Barcelona", "Milan", "Vienna", "Budapest", "Sofia", "Warsaw", "Helsinki", "Stuttgart")

types <- c("Food_PSM", "Proff_PSM", "Arts_PSM", "Shops_PSM", "Nightlife_PSM")

datasets <- c(Rome_venueheterogeneity, Barce_venueheterogeneity, Milan_venueheterogeneity, Vienna_venueheterogeneity )

coeff_compare[1, 1:2] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[2, 1:2] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[3, 1:2] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[4, 1:2] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[5, 1:2] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[6, 1:2] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[7, 1:2] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[8, 1:2] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]
coeff_compare[9, 1:2] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Food_PSM", c("Estimate","significance")][2,]


coeff_compare[1, 3:4] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[2, 3:4] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[3, 3:4] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[4, 3:4] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[5, 3:4] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[6, 3:4] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[7, 3:4] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[8, 3:4] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]
coeff_compare[9, 3:4] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Proff_PSM", c("Estimate","significance")][2,]


coeff_compare[1, 5:6] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[2, 5:6] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[3, 5:6] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[4, 5:6] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[5, 5:6] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[6, 5:6] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[7, 5:6] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[8, 5:6] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]
coeff_compare[9, 5:6] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Arts_PSM", c("Estimate","significance")][2,]

coeff_compare[1, 7:8] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[2, 7:8] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[3, 7:8] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[4, 7:8] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[5, 7:8] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[6, 7:8] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[7, 7:8] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[8, 7:8] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]
coeff_compare[9, 7:8] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Shops_PSM", c("Estimate","significance")][2,]

coeff_compare[1, 9:10] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[2, 9:10] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[3, 9:10] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[4, 9:10] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[5, 9:10] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[6, 9:10] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[7, 9:10] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[8, 9:10] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]
coeff_compare[9, 9:10] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Nightlife_PSM", c("Estimate","significance")][2,]

 write.csv(coeff_compare, "coeff_compare_PSM.csv")

 coeff_compare <- as.data.frame(matrix(nrow = 9, ncol = 10))
 colnames(coeff_compare) <- c("Food_ebal_coeff", "Food_ebal_p",  "Proff_ebal_coeff", "Proff_ebal_p",  "Arts_ebal_coeff", "Arts_ebal_p", 
                              "Shops_ebal_coeff", "Shops_ebal_p",  "Nightlife_ebal_coeff", "Nightlife_ebal_p")
 rownames(coeff_compare) <- c("Rome", "Barcelona", "Milan", "Vienna", "Budapest", "Sofia", "Warsaw", "Helsinki", "Stuttgart")
 
 types <- c("Food_ebal", "Proff_ebal", "Arts_ebal", "Shops_ebal", "Nightlife_ebal")
 
 datasets <- c(Rome_venueheterogeneity, Barce_venueheterogeneity, Milan_venueheterogeneity, Vienna_venueheterogeneity )
 
 coeff_compare[1, 1:2] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[2, 1:2] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[3, 1:2] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[4, 1:2] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[5, 1:2] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[6, 1:2] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[7, 1:2] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[8, 1:2] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 coeff_compare[9, 1:2] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Food_ebal", c("Estimate","significance")][2,]
 
 
 coeff_compare[1, 3:4] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[2, 3:4] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[3, 3:4] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[4, 3:4] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[5, 3:4] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[6, 3:4] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[7, 3:4] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[8, 3:4] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 coeff_compare[9, 3:4] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Proff_ebal", c("Estimate","significance")][2,]
 
 
 coeff_compare[1, 5:6] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[2, 5:6] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[3, 5:6] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[4, 5:6] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[5, 5:6] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[6, 5:6] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[7, 5:6] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[8, 5:6] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 coeff_compare[9, 5:6] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Arts_ebal", c("Estimate","significance")][2,]
 
 coeff_compare[1, 7:8] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[2, 7:8] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[3, 7:8] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[4, 7:8] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[5, 7:8] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[6, 7:8] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[7, 7:8] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[8, 7:8] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 coeff_compare[9, 7:8] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Shops_ebal", c("Estimate","significance")][2,]
 
 coeff_compare[1, 9:10] <- Rome_venueheterogeneity[Rome_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[2, 9:10] <- Barce_venueheterogeneity[Barce_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[3, 9:10] <- Milan_venueheterogeneity[Milan_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[4, 9:10] <- Vienna_venueheterogeneity[Vienna_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[5, 9:10] <- Budapest_venueheterogeneity[Budapest_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[6, 9:10] <- Sofia_venueheterogeneity[Sofia_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[7, 9:10] <- Warsaw_venueheterogeneity[Warsaw_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[8, 9:10] <- Helsinki_venueheterogeneity[Helsinki_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 coeff_compare[9, 9:10] <- Stuttgart_venueheterogeneity[Stuttgart_venueheterogeneity$type == "Nightlife_ebal", c("Estimate","significance")][2,]
 
 write.csv(coeff_compare, "coeff_compare_ebal.csv")
 

###############################

Decentralization <- as.data.frame(matrix(nrow = 9, ncol = 2))
rownames(Decentralization) <- c( "Rome",
                                 "Barcelona",
                                 "Vienna",
                                 "Budapest",
                                 "Helsinki",
                                 "Stuttgart",
                                 "Sofia",
                                 "Warsaw",
                                 "Milan")
colnames(Decentralization) <- c("Before", "After")

#Rome
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Romecentroid_data)))
x[1, ] <- Romecentroid_data$rowsum_13
Decentralization$Before[1]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Romecentroid_data)))
x[1, ] <- Romecentroid_data$rowsum_16
Decentralization$After[1]  <-  entropy(x)

#Barcelona
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Barcecentroid_data)))
x[1, ] <- Barcecentroid_data$rowsum15
Decentralization$Before[2]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Barcecentroid_data)))
x[1, ] <- Barcecentroid_data$rowsum18
Decentralization$After[2]  <-  entropy(x)

#Vienna
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Viennacentroid_data)))
x[1, ] <- Viennacentroid_data$rowsum_16
Decentralization$Before[3]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Viennacentroid_data)))
x[1, ] <- Viennacentroid_data$rowsum_19
Decentralization$After[3]  <-  entropy(x) 

#Budapest
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Budapestcentroid_data)))
x[1, ] <- Budapestcentroid_data$rowsum14
Decentralization$Before[4]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Budapestcentroid_data)))
x[1, ] <- Budapestcentroid_data$rowsum17
Decentralization$After[4]  <-  entropy(x) 

#Helsinki
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Helsinkicentroid_data)))
x[1, ] <- Helsinkicentroid_data$rowsum_16
Decentralization$Before[5]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Helsinkicentroid_data)))
x[1, ] <- Helsinkicentroid_data$rowsum_19
Decentralization$After[5]  <-  entropy(x) 

#Stuttgart
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Stuttgartcentroid_data)))
x[1, ] <- Stuttgartcentroid_data$rowsum_15
Decentralization$Before[6]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Stuttgartcentroid_data)))
x[1, ] <- Stuttgartcentroid_data$rowsum_18
Decentralization$After[6]  <-  entropy(x) 

#Sofia
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Sofiacentroid_data)))
x[1, ] <- Sofiacentroid_data$rowsum_14
Decentralization$Before[7]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Sofiacentroid_data)))
x[1, ] <- Sofiacentroid_data$rowsum18
Decentralization$After[7]  <-  entropy(x) 

#Warsaw
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Warsawcentroid_data)))
x[1, ] <- Warsawcentroid_data$rowsum_14
Decentralization$Before[8]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Warsawcentroid_data)))
x[1, ] <- Warsawcentroid_data$rowsum18
Decentralization$After[8]  <-  entropy(x) 

#Milan
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Milancentroid_data)))
x[1, ] <- Milancentroid_data$rowsum_14
Decentralization$Before[9]  <-  entropy(x)
x <- as.data.frame(matrix(nrow = 1, ncol = nrow(Milancentroid_data)))
x[1, ] <- Milancentroid_data$rowsum17
Decentralization$After[9]  <-  entropy(x) 

colnames(Decentralization) <- c("Before subway expansion", "After subway expansion")

setwd("C:/Dokumente/Utrecht/Master_these/results")
write.csv(Decentralization, "Decentralization2.csv")


############################# plotting the openings of venue types over time differentiating over metroneighbourhod or not ####
#Rome
Romevenuehexagonbuffers_clean <- subset(Romevenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Romevenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "Target_FID")
Rome_metroamenities <- merge(Romevenuehexagonbuffers_clean, Rome_PSM_dataset, by= "Target_FID", all = F)
Rome_metroamenities <- subset(Rome_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, Year_of_op, metroneighbourhood))
Rome_metroamenities <- distinct(Rome_metroamenities)
Rome_metroamenities <- Rome_metroamenities[Rome_metroamenities$venue_openingyear > 2013 & Rome_metroamenities$venue_openingyear < 2020,]
for(i in 1:nrow(Rome_metroamenities)){
  if(as.numeric(Rome_metroamenities$venue_openingmonth[i]) < 10){
    Rome_metroamenities$venue_openingmonth[i] <-paste("0", as.character(Rome_metroamenities$venue_openingmonth[i]), sep = "")
  }
}
Rome_metroamenities$venueopeningdate <- paste(as.character(Rome_metroamenities$venue_openingmonth), "-01-", as.character(Rome_metroamenities$venue_openingyear, sep = ""))
Rome_metroamenities$venueopeningdate <- gsub(" ", "", Rome_metroamenities$venueopeningdate)
Rome_metroamenities$venueopeningdate <- as.Date.character(as.character(Rome_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
jpeg("Distribution of Amenity Openings in Rome.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Nightlife" & Rome_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE))+
  geom_density(n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Nightlife" & Rome_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE)+
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Food" & Rome_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Food" & Rome_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Proffesional" & Rome_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Proffesional" & Rome_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Shops & Services" & Rome_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Shops & Services" & Rome_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Arts & Entertainment" & Rome_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 72) +
  geom_density(data= Rome_metroamenities[which(Rome_metroamenities$Metacatego == "Arts & Entertainment" & Rome_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 72) +
  labs(title = "Density of total openings per venue types over time in Rome", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Helsinki
Helsinkivenuehexagonbuffers_clean <- subset(Helsinkivenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Helsinkivenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Helsinki_metroamenities <- merge(Helsinkivenuehexagonbuffers_clean, Helsinki_PSM_dataset, by= "ORIG_FID", all = F)
Helsinki_metroamenities <- subset(Helsinki_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, Opened, metroneighbourhood))
Helsinki_metroamenities <- distinct(Helsinki_metroamenities)
Helsinki_metroamenities <- Helsinki_metroamenities[Helsinki_metroamenities$venue_openingyear > 2014 & Helsinki_metroamenities$venue_openingyear < 2020,]
for(i in 1:nrow(Helsinki_metroamenities)){
  if(as.numeric(Helsinki_metroamenities$venue_openingmonth[i]) < 10){
    Helsinki_metroamenities$venue_openingmonth[i] <-paste("0", as.character(Helsinki_metroamenities$venue_openingmonth[i]), sep = "")
  }
}
Helsinki_metroamenities$venueopeningdate <- paste(as.character(Helsinki_metroamenities$venue_openingmonth), "-01-", as.character(Helsinki_metroamenities$venue_openingyear, sep = ""))
Helsinki_metroamenities$venueopeningdate <- gsub(" ", "", Helsinki_metroamenities$venueopeningdate)
Helsinki_metroamenities$venueopeningdate <- as.Date.character(as.character(Helsinki_metroamenities$venueopeningdate), format = '%m-%d-%Y')

setwd("C:/Dokumente/Utrecht/Master_these/Data/Helsinki/PSMresults")
jpeg("Distribution of Amenity Openings in Helsinki.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Nightlife" & Helsinki_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Nightlife" & Helsinki_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Food" & Helsinki_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Food" & Helsinki_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Professional" & Helsinki_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Professional" & Helsinki_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Shops & Services" & Helsinki_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Shops & Services" & Helsinki_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Art & Entertainment" & Helsinki_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Helsinki_metroamenities[which(Helsinki_metroamenities$Metacatego == "Art & Entertainment" & Helsinki_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Helsinki", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Vienna
Viennavenuehexagonbuffers_clean <- subset(Viennavenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Viennavenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Vienna_metroamenities <- merge(Viennavenuehexagonbuffers_clean, Vienna_PSM_dataset, by= "ORIG_FID", all = F)
Vienna_metroamenities <- subset(Vienna_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, station_openingyear, metroneighbourhood))
Vienna_metroamenities <- distinct(Vienna_metroamenities)
Vienna_metroamenities <- Vienna_metroamenities[Vienna_metroamenities$venue_openingyear > 2014 & Vienna_metroamenities$venue_openingyear < 2020,]
for(i in 1:nrow(Vienna_metroamenities)){
  if(as.numeric(Vienna_metroamenities$venue_openingmonth[i]) < 10){
    Vienna_metroamenities$venue_openingmonth[i] <-paste("0", as.character(Vienna_metroamenities$venue_openingmonth[i]), sep = "")
  }
}
Vienna_metroamenities$venueopeningdate <- paste(as.character(Vienna_metroamenities$venue_openingmonth), "-01-", as.character(Vienna_metroamenities$venue_openingyear, sep = ""))
Vienna_metroamenities$venueopeningdate <- gsub(" ", "", Vienna_metroamenities$venueopeningdate)
Vienna_metroamenities$venueopeningdate <- as.Date.character(as.character(Vienna_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Vienna/PSMresults")
jpeg("Distribution of Amenity Openings in Vienna.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Nightlife" & Vienna_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Nightlife" & Vienna_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Food" & Vienna_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Food" & Vienna_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Professional" & Vienna_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Professional" & Vienna_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Shops & Services" & Vienna_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Shops & Services" & Vienna_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Art & Entertainment" & Vienna_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Vienna_metroamenities[which(Vienna_metroamenities$Metacatego == "Art & Entertainment" & Vienna_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Vienna", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Stuttgart
Stuttgartvenuehexagonbuffers_clean <- subset(Stuttgartvenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Stuttgartvenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Stuttgart_metroamenities <- merge(Stuttgartvenuehexagonbuffers_clean, Stuttgart_PSM_dataset, by= "ORIG_FID", all = F)
Stuttgart_metroamenities <- subset(Stuttgart_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, openingyea, metroneighbourhood))
Stuttgart_metroamenities <- distinct(Stuttgart_metroamenities)
Stuttgart_metroamenities <- Stuttgart_metroamenities[Stuttgart_metroamenities$venue_openingyear > 2014 & Stuttgart_metroamenities$venue_openingyear < 2020,]
for(i in 1:nrow(Stuttgart_metroamenities)){
  if(as.numeric(Stuttgart_metroamenities$venue_openingmonth[i]) < 10){
    Stuttgart_metroamenities$venue_openingmonth[i] <-paste("0", as.character(Stuttgart_metroamenities$venue_openingmonth[i]), sep = "")
  }
}
Stuttgart_metroamenities$venueopeningdate <- paste(as.character(Stuttgart_metroamenities$venue_openingmonth), "-01-", as.character(Stuttgart_metroamenities$venue_openingyear, sep = ""))
Stuttgart_metroamenities$venueopeningdate <- gsub(" ", "", Stuttgart_metroamenities$venueopeningdate)
Stuttgart_metroamenities$venueopeningdate <- as.Date.character(as.character(Stuttgart_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Stuttgart/PSMresults")
jpeg("Distribution of Amenity Openings in Stuttgart.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Nightlife" & Stuttgart_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Nightlife" & Stuttgart_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Food" & Stuttgart_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Food" & Stuttgart_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Professional" & Stuttgart_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Professional" & Stuttgart_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Shops & Services" & Stuttgart_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Shops & Services" & Stuttgart_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Art & Entertainment" & Stuttgart_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Stuttgart_metroamenities[which(Stuttgart_metroamenities$Metacatego == "Art & Entertainment" & Stuttgart_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Stuttgart", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Barcelona
Barcevenuehexagonbuffers_clean <- subset(Barcevenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Barcevenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Barce_metroamenities <- merge(Barcevenuehexagonbuffers_clean, Barce_PSM_dataset, by= "ORIG_FID", all = F)
Barce_metroamenities <- subset(Barce_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, station_Openingyear, metroneighbourhood))
Barce_metroamenities <- distinct(Barce_metroamenities)
Barce_metroamenities <- Barce_metroamenities[Barce_metroamenities$venue_openingyear > 2014 & Barce_metroamenities$venue_openingyear < 2020,]

Barce_metroamenities$venueopeningdate <- paste(as.character(Barce_metroamenities$venue_openingmonth), "-01-", as.character(Barce_metroamenities$venue_openingyear, sep = ""))
Barce_metroamenities$venueopeningdate <- gsub(" ", "", Barce_metroamenities$venueopeningdate)
Barce_metroamenities$venueopeningdate <- as.Date.character(as.character(Barce_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Barcelona/PSMresults")
jpeg("Distribution of Amenity Openings in Barcelona.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Food" & Barce_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Food" & Barce_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Proffesional" & Barce_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Proffesional" & Barce_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Shops & Services" & Barce_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Barce_metroamenities[which(Barce_metroamenities$Metacatego == "Shops & Services" & Barce_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Barcelona", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Budapest
Budapestvenuehexagonbuffers_clean <- subset(Budapestvenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, ORIG_FID))
colnames(Budapestvenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Budapest_metroamenities <- merge(Budapestvenuehexagonbuffers_clean, Budapest_PSM_dataset, by= "ORIG_FID", all = F)
Budapest_metroamenities <- subset(Budapest_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, Year_of_op, metroneighbourhood))
Budapest_metroamenities <- distinct(Budapest_metroamenities)
Budapest_metroamenities <- Budapest_metroamenities[Budapest_metroamenities$venue_openingyear > 2012 & Budapest_metroamenities$venue_openingyear < 2020,]

Budapest_metroamenities$venueopeningdate <- paste(as.character(Budapest_metroamenities$venue_openingmonth), "-01-", as.character(Budapest_metroamenities$venue_openingyear, sep = ""))
Budapest_metroamenities$venueopeningdate <- gsub(" ", "", Budapest_metroamenities$venueopeningdate)
Budapest_metroamenities$venueopeningdate <- as.Date.character(as.character(Budapest_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Budapest/PSMresults")
jpeg("Distribution of Amenity Openings in Budapest.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Nightlife" & Budapest_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Nightlife" & Budapest_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Food" & Budapest_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Food" & Budapest_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Professional" & Budapest_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Professional" & Budapest_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Shops & Services" & Budapest_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Shops & Services" & Budapest_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Art & Entertainment" & Budapest_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Budapest_metroamenities[which(Budapest_metroamenities$Metacatego == "Art & Entertainment" & Budapest_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Budapest", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Warsaw
Warsawvenuehexagonbuffers_clean <- subset(Warsawvenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryna, Metacatego, creationye, creationmo, ORIG_FID))
colnames(Warsawvenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Warsaw_metroamenities <- merge(Warsawvenuehexagonbuffers_clean, Warsaw_PSM_dataset, by= "ORIG_FID", all = F)
Warsaw_metroamenities <- subset(Warsaw_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, station_openingyear, metroneighbourhood))
Warsaw_metroamenities <- distinct(Warsaw_metroamenities)
Warsaw_metroamenities <- Warsaw_metroamenities[Warsaw_metroamenities$venue_openingyear > 2013 & Warsaw_metroamenities$venue_openingyear < 2020,]
Warsaw_metroamenities$venueopeningdate <- paste(as.character(Warsaw_metroamenities$venue_openingmonth), "-01-", as.character(Warsaw_metroamenities$venue_openingyear, sep = ""))
Warsaw_metroamenities$venueopeningdate <- gsub(" ", "", Warsaw_metroamenities$venueopeningdate)
Warsaw_metroamenities$venueopeningdate <- as.Date.character(as.character(Warsaw_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/PSMresults")
jpeg("Distribution of Amenity Openings in Warsaw.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Nightlife" & Warsaw_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Nightlife" & Warsaw_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Food" & Warsaw_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Food" & Warsaw_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Professional" & Warsaw_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Professional" & Warsaw_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Shops & Services" & Warsaw_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Shops & Services" & Warsaw_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Arts & Entertainment" & Warsaw_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Warsaw_metroamenities[which(Warsaw_metroamenities$Metacatego == "Arts & Entertainment" & Warsaw_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Warsaw", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


#Sofia
Sofiavenuehexagonbuffers_clean <- subset(Sofiavenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryname, Metacategory, creationyear, creationmonths, ORIG_FID))
colnames(Sofiavenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Sofia_metroamenities <- merge(Sofiavenuehexagonbuffers_clean, Sofia_PSM_dataset, by= "ORIG_FID", all = F)
Sofia_metroamenities <- subset(Sofia_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, Openingyear, metroneighbourhood))
Sofia_metroamenities <- distinct(Sofia_metroamenities)
Sofia_metroamenities <- Sofia_metroamenities[Sofia_metroamenities$venue_openingyear > 2013 & Sofia_metroamenities$venue_openingyear < 2020,]
Sofia_metroamenities$venueopeningdate <- paste(as.character(Sofia_metroamenities$venue_openingmonth), "-01-", as.character(Sofia_metroamenities$venue_openingyear, sep = ""))
Sofia_metroamenities$venueopeningdate <- gsub(" ", "", Sofia_metroamenities$venueopeningdate)
Sofia_metroamenities$venueopeningdate <- as.Date.character(as.character(Sofia_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Sofia/PSMresults")
jpeg("Distribution of Amenity Openings in Sofia.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Nightlife" & Sofia_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Nightlife" & Sofia_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Food" & Sofia_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Food" & Sofia_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Professional" & Sofia_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Professional" & Sofia_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Shops & Services" & Sofia_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Shops & Services" & Sofia_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Arts & Entertainment" & Sofia_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Sofia_metroamenities[which(Sofia_metroamenities$Metacatego == "Arts & Entertainment" & Sofia_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Sofia", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()

#Milan
Milanvenuehexagonbuffers_clean <- subset(Milanvenuehexagonbuffers_tot, select = c(venueid, venuename, categoryid, categoryname, Metacategory, creationyear, creationmonths, ORIG_FID))
colnames(Milanvenuehexagonbuffers_clean) <- c("venueid", "venuename", "categoryid", "categoryna", "Metacatego", "venue_openingyear", "venue_openingmonth", "ORIG_FID")
Milan_metroamenities <- merge(Milanvenuehexagonbuffers_clean, Milan_PSM_dataset, by= "ORIG_FID", all = F)
Milan_metroamenities <- subset(Milan_metroamenities, select= c(venueid, venuename, categoryid, categoryna, Metacatego, venue_openingyear, venue_openingmonth, Openingyear, metroneighbourhood))
Milan_metroamenities <- distinct(Milan_metroamenities)
Milan_metroamenities <- Milan_metroamenities[Milan_metroamenities$venue_openingyear > 2013 & Milan_metroamenities$venue_openingyear < 2020,]
Milan_metroamenities$venueopeningdate <- paste(as.character(Milan_metroamenities$venue_openingmonth), "-01-", as.character(Milan_metroamenities$venue_openingyear, sep = ""))
Milan_metroamenities$venueopeningdate <- gsub(" ", "", Milan_metroamenities$venueopeningdate)
Milan_metroamenities$venueopeningdate <- as.Date.character(as.character(Milan_metroamenities$venueopeningdate), format = '%m-%d-%Y')


setwd("C:/Dokumente/Utrecht/Master_these/Data/Milan/PSMresults")
jpeg("Distribution of Amenity Openings in Milan.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Nightlife" & Milan_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Nightlife", show.legend = TRUE, trim = TRUE, n =60))+
  geom_density(n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Nightlife" & Milan_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Nightlife"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n=60)+
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Food" & Milan_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Food"), show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Food" & Milan_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Food"), linetype = "dotted", show.legend = TRUE, trim = TRUE, n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Professional" & Milan_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Professional"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Professional" & Milan_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Professional"),linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Shops & Services" & Milan_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Shops and Services"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Shops & Services" & Milan_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Shops and Services"), linetype = "dotted", show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Art & Entertainment" & Milan_metroamenities$metroneighbourhood == 1),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), show.legend = TRUE, trim = TRUE,n = 60) +
  geom_density(data= Milan_metroamenities[which(Milan_metroamenities$Metacatego == "Art & Entertainment" & Milan_metroamenities$metroneighbourhood == 0),], aes(x=venueopeningdate, colour ="Arts and Entertainment"), linetype = "dotted",show.legend = TRUE, trim = TRUE,n = 60) +
  labs(title = "Density of total openings per venue types over time in Milan", subtitle = "(Solid Line = Metroneighbourhoods, Dotted Line = PSM Control Areas)")+
  labs(y="Density")+
  labs(x="Openingdate") + 
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(name = "Venue Types:",
                     breaks = c("Food", "Shops and Services", "Professional", "Nightlife","Arts and Entertainment"),
                     values = c("Food" = "cyan",  "Shops and Services" = "firebrick1", "Professional" = "blue", "Nightlife" = "green2", "Arts and Entertainment" = "purple"),
                     aesthetics = c("colour", "fill"))
dev.off()


################################################
############ Dissimilarity Index ###############
################################################

# Rome
setwd("C:/Dokumente/Utrecht/Master_these/Data/Rome/PSMresults")
jpeg("Rome Decentralization before.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Romecentroid_data, aes(x= centroid_distance_center, y = rowsum_13))+
  geom_point( alpha = 0.1 , colour = "black")+
  geom_smooth(method='lm', colour = "red")+
  labs(title = "Rome before the Subway Expansion")+
  labs(y="Density of Amenities")+
  labs(x="Distance to CBD")
dev.off()

jpeg("Rome Decentralization after.jpeg",  width=7, height=6, units = 'in', res = 2000 )
ggplot(data =  Romecentroid_data, aes(x= centroid_distance_center, y = rowsum_16))+
  geom_point( alpha = 0.1 , colour = "black")+
  geom_smooth(method='lm', colour = "red")+
  labs(title = "Rome 3 years after the Subway Expansion")+
  labs(y="Density of Amenities")+
  labs(x="Distance to CBD")
dev.off()

Decentralization_results <- matrix(ncol = 4, nrow = 9)
colnames(Decentralization_results) <- c("City", "beta_densdist_before", "beta_densdist_after", "beta_diff")
Decentralization_results <- as.data.frame(Decentralization_results)
Decentralization_results$City <- c("Rome", "Barcelona", "Helsinki", "Vienna", "Milan", "Budapest", "Sofia", "Warsaw", "Stuttgart")
Decentralization_results$beta_densdist_before[1] <- lm(rowsum_13 ~ centroid_distance_center, Romecentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[1] <- lm(rowsum_16 ~ centroid_distance_center, Romecentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[2] <- lm(rowsum15 ~ centroid_distance_center, Barcecentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[2] <- lm(rowsum18 ~ centroid_distance_center, Barcecentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[3] <- lm(rowsum_16 ~ centroid_dist_center, Helsinkicentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[3] <- lm(rowsum_19 ~ centroid_dist_center, Helsinkicentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[4] <- lm(rowsum_16 ~ distance_center, Viennacentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[4] <- lm(rowsum_19 ~ distance_center, Viennacentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[5] <- lm(rowsum_14 ~ centroid_dist_center, Milancentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[5] <- lm(rowsum17~ centroid_dist_center, Milancentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[6] <- lm(rowsum14 ~ distance_center, Budapestcentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[6] <- lm(rowsum17~ distance_center, Budapestcentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[7] <- lm(rowsum_14 ~ centroid_dist_center, Sofiacentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[7] <- lm(rowsum18~ centroid_dist_center, Sofiacentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[8] <- lm(rowsum_14 ~ centroid_dist_center, Warsawcentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[8] <- lm(rowsum18~ centroid_dist_center, Warsawcentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_before[9] <- lm(rowsum_15 ~ centroid_dist_center, Stuttgartcentroid_data)$coefficients[2]
Decentralization_results$beta_densdist_after[9] <- lm(rowsum_18~ centroid_dist_center, Stuttgartcentroid_data)$coefficients[2]

Decentralization_results$beta_diff <- Decentralization_results$beta_densdist_after - Decentralization_results$beta_densdist_before


setwd("C:/Dokumente/Utrecht/Master_these/Data")
write.csv(Decentralization_results, "Decentralization_results.csv")
