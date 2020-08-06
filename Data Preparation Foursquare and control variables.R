pkgs <- c("maptools","rgdal","sp", "sf", "jpeg", "textreadr", "GEOmap", "RANN", "NISTunits",  "pracma", "celestial", "openxlsx", "data.table", "purrr",  "geonames",  "foreign", "tidyverse",  
          "ggplot2", "cowplot",  "psych",  "tidyimpute", "sjPlot", "ggpubr", "devtools", "car", "lattice",  "openintro","dbf" , 
          "raster",  "spatialEco","xlsx", "spData",  "leaflet",  "geojsonio",  "rjson", "RJSONIO", "jsonlite","EconGeo", "finalfit",  
          "rstan", "boot",  "plm", "phylin")
sapply(pkgs, require, character.only = T) #load 
rm(pkgs)

####### Code written for Warsaw, but is equally applied to Rome, Barcelona, Vienna, Budapest, Milan, Helsinki, Stuttgart and Sofia


###########################################################################
########################## Reading the Foursquare data ####################
###########################################################################

setwd("C:/Dokumente/Utrecht/Master_these/Data/Foursquare/Warsaw/inferred_openings")
Warsaw_Fsq =read.csv("Warsaw_all.csv", header = T)
coordinates(Warsaw_Fsq)= ~lon+lat
proj4string(Warsaw_Fsq)=CRS("+proj=longlat +datum=WGS84") 
writeOGR(Warsaw_Fsq, dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw" ,layer="Warsaw_Fsq",driver="ESRI Shapefile")
Warsaw_Fsq <- readOGR(dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw" ,layer="Warsaw_Fsq")
plot(Warsaw_Fsq)



#######################################################################################################
####### getting the Nightlight defined urban area from Nightlight Satelite imagery ###################
#######################################################################################################

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/Nightlight")
Warsaw_Nightlight <- raster("LuoJia1-01_LR201808294726_20180828203503_HDR_0019_gec.tif")

x_coord <- c(21.6, 21.6, 20.0, 20.0)
y_coord <- c(51.5, 53.0, 53.0, 51.5)
pol <- cbind(x_coord, y_coord)
poly <- Polygon(pol)
ps = Polygons(list(poly),1)
Clip_extent = SpatialPolygons(list(ps))
plot(Clip_extent)

Warsaw_Nightlight <- crop(Warsaw_Nightlight, extent(Clip_extent))
summary(Warsaw_Nightlight)
plot(Warsaw_Nightlight)

jpeg("Warsaw_Nightlight_hist.jpeg")
Warsaw_Nightlight_hist <- hist(as.integer(Warsaw_Nightlight),
     main = "Distribution of Nightlight values Warsaw",
     breaks = 500, xlim = c(0,60000), 
     xlab = "Light Gradient", ylab = "Frequency",
     col = "springgreen")

dev.off() 

## determine the cutoff value out of the distribution --> should be the infliction value
cutoff_value <- 1000
Warsaw_NLight_def_urbanarea <- Warsaw_Nightlight > cutoff_value
plot(Warsaw_NLight_def_urbanarea)

Warsaw_urbanarea <- rasterToPolygons(Warsaw_NLight_def_urbanarea,fun=function(x){x == 1}, dissolve= F)
Warsaw_urbanarea  <- aggregate(Warsaw_urbanarea, dissolve = T)
Warsaw_urbanarea <- disaggregate(Warsaw_urbanarea)
Warsaw_urbanarea$area <- areaPolygon(Warsaw_urbanarea)
biggest <- which(Warsaw_urbanarea$area == max(Warsaw_urbanarea$area))
Warsaw_urbanarea <-Warsaw_urbanarea[biggest,]
jpeg("Warsaw_urbanarea.jpeg")
plot(Warsaw_urbanarea)
dev.off() 

plot(Warsaw_urbanarea, pch = 20, col = "purple")
plot(Warsaw_Fsq, add = T)


writeOGR(Warsaw_urbanarea, dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw",layer="Warsaw_urbanarea",driver="ESRI Shapefile")
Warsaw_urbanarea <- readOGR(dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw",layer="Warsaw_urbanarea")


############################################################################################
############################# Creating regular point grid in the urban area #########################
############################################################################################

polygoncoord <- as.data.frame(Warsaw_urbanarea@polygons[[1]]@Polygons[[1]]@coords)
names(polygoncoord) <- c("x", "y")
max_coord_x <- max(polygoncoord$x)
min_coord_x <- min(polygoncoord$x)
max_coord_y <- max(polygoncoord$y)
min_coord_y <- min(polygoncoord$y)
necessary_coord_x <- as.integer((max_coord_x - min_coord_x)/ 0.001)
necessary_coord_y <- as.integer((max_coord_y - min_coord_y)/ 0.001)
xcoordinates <- c(min_coord_x)
for(i in 1:necessary_coord_x){
  xcoordinates <- append(xcoordinates, (min_coord_x + (0.001*i)))
}
ycoordinates <- c(min_coord_y)
for(i in 1:necessary_coord_y){
  ycoordinates <- append(ycoordinates, (min_coord_y + (0.001*i)))
}

Warsawcentroids <- matrix(data = 0, nrow = (length(xcoordinates)*length(ycoordinates)), ncol =  3)
Warsawcentroids <- as.data.frame(Warsawcentroids)
colnames(Warsawcentroids) <- c("id", "x", "y")
s<-0
for(i in 1:length(xcoordinates)){
  for(n in 1:length(ycoordinates)){
    Warsawcentroids$x[(s+n)] <- xcoordinates[i]
    Warsawcentroids$y[(s+n)] <- ycoordinates[n]
  }
  s <- s + n
}
coordinates(Warsawcentroids)= ~x+y
proj4string(Warsawcentroids)=CRS("+proj=longlat +datum=WGS84")
plot(Warsawcentroids)

Warsawcentroids <- point.in.poly(Warsawcentroids, Warsaw_urbanarea, sp = T)
Warsawcentroids <- Warsawcentroids[!is.na(Warsawcentroids$area),]
Warsawcentroids$ORIG_FID <- c(seq(1:nrow(Warsawcentroids)))
plot(Warsawcentroids)


setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
writeOGR(Warsawcentroids, dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw",layer="Warsawcentroids",driver="ESRI Shapefile")
Warsawcentroids <- readOGR(dsn="C:/Dokumente/Utrecht/Master_these/Data/Warsaw",layer="Warsawcentroids")


######################################################################################################
################## Joining the points with Worldpop data and the metrostation shapefile #########
######################################################################################################
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/worldpop")
Warsaw_pop2018 <- raster("pol_ppp_2018.tif")
Warsaw_pop2012 <- raster("pol_ppp_2012.tif")

x_coord <- c(21.6, 21.6, 20.0, 20.0)
y_coord <- c(51.5, 53.0, 53.0, 51.5)
pol <- cbind(x_coord, y_coord)
poly <- Polygon(pol)
ps = Polygons(list(poly),1)
Clip_extent = SpatialPolygons(list(ps))
plot(Clip_extent)

Warsaw_pop2018 <- crop(Warsaw_pop2018, extent(Clip_extent))
Warsaw_pop2012 <- crop(Warsaw_pop2012, extent(Clip_extent))

idw(values, coords, grid, method = "Shepard", p = 2, R = 2, N = 15,
    distFUN = geo.dist, ...)

Warsaw_pop2018 <- as.integer(Warsaw_pop2018)
Warsaw_pop2012 <- as.integer(Warsaw_pop2012)
writeRaster(Warsaw_pop2012, "Warsaw_pop2012.tif", format = "GTiff")
writeRaster(Warsaw_pop2018, "Warsaw_pop2018.tif", format = "GTiff")
Warsaw_pop2018 <- rasterToPolygons(Warsaw_pop2018, dissolve= F)
Warsaw_pop2012 <- rasterToPolygons(Warsaw_pop2012, dissolve= F)
Warsaw_pop2018$pop_2018 <- Warsaw_pop2018$layer
Warsaw_pop2012$pop_2012 <- Warsaw_pop2012$layer


Warsawcentroid_data <- point.in.poly(Warsawcentroid_data, Warsaw_pop2018, sp = T, duplicate = T)
Warsawcentroid_data <- point.in.poly(Warsawcentroid_data, Warsaw_pop2012, sp = T, duplicate = T)

remove(Warsaw_pop2012, Warsaw_pop2018)
Warsawcentroid_data$popchange <- Warsawcentroid_data$pop_2018 - Warsawcentroid_data$pop_2012

###
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw/transport")
Warsaw_subwaysystem =read.csv("metrosystemWarsaw.csv", header = T)
coordinates(Warsaw_subwaysystem)= ~lon+lat
proj4string(Warsaw_subwaysystem)=CRS("+proj=longlat +datum=WGS84") # set it to UTM

plot(Warsaw_urbanarea)
plot(Warsaw_subwaysystem, pch = 20, col = "purple", add = TRUE)

d <- pointDistance(Warsawcentroids, Warsaw_subwaysystem, lonlat=TRUE)
r <- apply(d, 1, which.min)
Warsawcentroid_data <- as.data.frame(Warsawcentroids)
Warsawcentroid_data$metroid <- r
Warsaw_subwaysystem =read.csv("metrosystemWarsaw.csv", header = T)
Warsaw_subwaysystem$metroid <- c(seq(1,nrow(Warsaw_subwaysystem)))
Warsawcentroid_data <- merge(Warsawcentroid_data, Warsaw_subwaysystem, by= "metroid", all = T)
Warsawcentroid_data <- Warsawcentroid_data[!is.na(Warsawcentroid_data$ORIG_FID),]
Warsawcentroid_data <- Warsawcentroid_data[order(Warsawcentroid_data$ORIG_FID),]
for(i in 1:nrow(Warsawcentroid_data)){
  Warsawcentroid_data$NEAR_DIST[i] <- d[i, Warsawcentroid_data$metroid[i]]
}

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
write.csv(as.data.frame(Warsawcentroid_data), "Warsawcentroid_data.csv")
Warsawcentroid_data <- read.csv("Warsawcentroid_data.csv", header = T)

##############################################################################################################
################ aggregating population data in the 300m buffers of the respective centroid #################
##############################################################################################################

Warsawcentroids <- readOGR(dsn = "C:/Dokumente/Utrecht/Master_these/Data/Warsaw", layer = "Warsawcentroids")
Warsawcentroids <- spTransform(Warsawcentroids, CRS("+proj=longlat +datum=WGS84" ))

make_GeodesicBuffer <- function(pts, width) {
  ### A) Construct buffers as points at given distance and bearing
  # a vector of bearings (fallows a circle)
  dg <- seq(from = 0, to = 360, by = 5)
  
  # Construct equidistant points defining circle shapes (the "buffer points")
  buff.XY <- geosphere::destPoint(p = pts, 
                                  b = rep(dg, each = length(pts)), 
                                  d = width)
  
  ### B) Make SpatialPolygons
  # group (split) "buffer points" by id
  buff.XY <- as.data.frame(buff.XY)
  id  <- rep(1:length(pts), times = length(dg))
  lst <- split(buff.XY, id)
  
  # Make SpatialPolygons out of the list of coordinates
  poly   <- lapply(lst, sp::Polygon, hole = FALSE)
  polys  <- lapply(list(poly), sp::Polygons, ID = NA)
  spolys <- sp::SpatialPolygons(Srl = polys, 
                                proj4string = CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84")))
  # Disaggregate (split in unique polygons)
  spolys <- sp::disaggregate(spolys)
  return(spolys)
}

Warsawcentroids_buff300m  <- make_GeodesicBuffer(Warsawcentroids, width=300)
Warsawcentroids_buff300m$ORIG_FID_buff <- Warsawcentroids$ORIG_FID

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
Warsawcentroid_data <- read.csv("Warsawcentroid_data.csv", header = TRUE)
Warsaw_centr <- Warsawcentroid_data[,c("ORIG_FID", "centroid_lat", "centroid_lon", "pop_2012", "pop_2018")]
coordinates(Warsaw_centr) =~ centroid_lat + centroid_lon
proj4string(Warsaw_centr)=CRS("+proj=longlat +datum=WGS84")
Warsaw_centr <- point.in.poly(Warsaw_centr, Warsawcentroids_buff300m, sp = T, duplicate = T)


nrow(Warsaw_centr)/nrow(Warsawcentroid_data)
for(i in 1:nrow(Warsawcentroid_data)){
  Warsawcentroid_data$pop_2012_buffer[Warsawcentroid_data$ORIG_FID == i] <- sum(Warsaw_centr$pop_2012[Warsaw_centr$ORIG_FID_buff == i])
  Warsawcentroid_data$pop_2018_buffer[Warsawcentroid_data$ORIG_FID == i] <- sum(Warsaw_centr$pop_2018[Warsaw_centr$ORIG_FID_buff == i])
}

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")



##############################################################################################################
################ creating buffers around the centroids and full joining it with the Fsq Data #################
##############################################################################################################


Warsawcentroids <- readOGR(dsn = "C:/Dokumente/Utrecht/Master_these/Data/Warsaw", layer = "Warsawcentroids")
Warsawcentroids <- spTransform(Warsawcentroids, CRS("+proj=longlat +datum=WGS84" ))
plot(Warsawcentroids)

Warsaw_Fsq <- readOGR(dsn = "C:/Dokumente/Utrecht/Master_these/Data/Foursquare/Warsaw", layer = "Warsaw_Fsq")
names(Warsaw_Fsq) <- c("id", "categoryna", "categoryid", "venuename",  "venueid", "Latitude", "Longitude", "Metacatego", "openyear", "openmonths", "NA")


make_GeodesicBuffer <- function(pts, width) {
  ### A) Construct buffers as points at given distance and bearing
  # a vector of bearings (fallows a circle)
  dg <- seq(from = 0, to = 360, by = 5)
  
  # Construct equidistant points defining circle shapes (the "buffer points")
  buff.XY <- geosphere::destPoint(p = pts, 
                                  b = rep(dg, each = length(pts)), 
                                  d = width)
  
  ### B) Make SpatialPolygons
  # group (split) "buffer points" by id
  buff.XY <- as.data.frame(buff.XY)
  id  <- rep(1:length(pts), times = length(dg))
  lst <- split(buff.XY, id)
  
  # Make SpatialPolygons out of the list of coordinates
  poly   <- lapply(lst, sp::Polygon, hole = FALSE)
  polys  <- lapply(list(poly), sp::Polygons, ID = NA)
  spolys <- sp::SpatialPolygons(Srl = polys, 
                                proj4string = CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84")))
  # Disaggregate (split in unique polygons)
  spolys <- sp::disaggregate(spolys)
  return(spolys)
}

Warsawcentroids_buff300m  <- make_GeodesicBuffer(Warsawcentroids, width=300)
Warsawcentroids_buff300m$ORIG_FID <- Warsawcentroids$ORIG_FID

plot(Warsawcentroids_buff300m)
plot(Warsaw_subwaysystem, pch = 20, col = "purple", add = TRUE)
plot(Warsaw_Fsq, add = T, col = "purple")
plot(Sofiacentroids_buff300m, add = T)

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
Warsaw_bufferjoin_300m <- point.in.poly(Warsaw_Fsq, Warsawcentroids_buff300m, sp = T, duplicate = T)
write.csv(as.data.frame(Warsaw_bufferjoin_300m), "Warsawvenuehexagonbuffers_tot.csv")


#############################################################################################################################
#################### Transforming into bufferunit - venuecategory matrix for before and after the station openings ##########
############################################################################################################################

Warsawvenuehexagonbuffers_tot <- read.csv("Warsawvenuehexagonbuffers_tot.csv", header = TRUE)
Warsawvenuehexagonbuffers_tot <- subset(Warsawvenuehexagonbuffers_tot, select = -c( X, Field1, Metacatego, coords.x1, coords.x2))
Warsawvenuehexagonbuffers_tot <- distinct(Warsawvenuehexagonbuffers_tot)

venuecategories <- subset(Warsawvenuehexagonbuffers_tot, select = c("categoryid", "categoryna"))
venuecategories <- unique(venuecategories)
categorynames <- venuecategories$categoryid

columnnames <- paste("cat_" , categorynames, sep = "")
print(columnnames)


Warsawcat_matrix <- matrix(data = 0, nrow = nrow(Warsawcentroids), ncol = length(columnnames))
colnames(Warsawcat_matrix) <- columnnames

Warsawvenuecat_matrix_pre2018 <- Warsawcat_matrix
subset2018 <- Warsawvenuehexagonbuffers_tot[Warsawvenuehexagonbuffers_tot$creationye < 2018,]
subset2018 <- subset2018[!is.na(subset2018$ORIG_FID),]
for(i in 1:length(categorynames)){
  x <- subset2018[subset2018$categoryid == categorynames[i],]
  y <- x %>% count(ORIG_FID)
  Warsawvenuecat_matrix_pre2018[y$ORIG_FID, i] <- y$n
}
write.csv(Warsawvenuecat_matrix_pre2018,"Warsawvenuecat_matrix_pre2018.csv")


Warsawvenuecat_matrix_pre2015 <- Warsawcat_matrix
subset2015 <- Warsawvenuehexagonbuffers_tot[Warsawvenuehexagonbuffers_tot$creationye < 2015,]
subset2015 <- subset2015[!is.na(subset2015$ORIG_FID),]
for(i in 1:length(categorynames)){
  x <- subset2015[subset2015$categoryid == categorynames[i],]
  y <- x %>% count(ORIG_FID)
  Warsawvenuecat_matrix_pre2015[y$ORIG_FID, i] <- y$n
}
write.csv(Warsawvenuecat_matrix_pre2015,"Warsawvenuecat_matrix_pre2015.csv")


#########################################################################################################
################### density and basic SHannon Wiener entropy ############################################
#########################################################################################################
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")

entropy <- function(mat) {
  freqs <- mat/rowSums (mat)
  entropy <- - rowSums (freqs * log2(freqs+0.000000001))
  entropy <- round (entropy, digits = 3)
  return (entropy)
}

Warsawvenuecat_matrix_pre2018 <- read.csv("Warsawvenuecat_matrix_pre2018.csv", header = T)
Warsawvenuecat_matrix_pre2018 <- subset(Warsawvenuecat_matrix_pre2018, select = -c(X))
Warsawvenuecat_matrix_pre2015 <- read.csv("Warsawvenuecat_matrix_pre2015.csv", header = T)
Warsawvenuecat_matrix_pre2015 <- subset(Warsawvenuecat_matrix_pre2015, select = -c(X))

Warsawvenuecat_matrix_pre2015$entropy <- entropy(Warsawvenuecat_matrix_pre2015)
Warsawvenuecat_matrix_pre2015$rowsum <- rowSums(Warsawvenuecat_matrix_pre2015[, 1:length(columnnames)])
Warsawcentroid_data$entropy_14 <- Warsawvenuecat_matrix_pre2015$entropy[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$rowsum_14 <- Warsawvenuecat_matrix_pre2015$rowsum[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$entropy_14[Warsawcentroid_data$entropy_14 == "NaN"] <- 0

Warsawvenuecat_matrix_pre2018$entropy <- entropy(Warsawvenuecat_matrix_pre2018)
Warsawvenuecat_matrix_pre2018$rowsum <- rowSums(Warsawvenuecat_matrix_pre2018[, 1:length(columnnames)])
Warsawcentroid_data$entropy18 <- Warsawvenuecat_matrix_pre2018$entropy[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$rowsum18 <- Warsawvenuecat_matrix_pre2018$rowsum[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$entropy18[Warsawcentroid_data$entropy18 == "NaN"] <- 0

Warsawcentroid_data$entropy_diff <- (Warsawcentroid_data$entropy18 - Warsawcentroid_data$entropy_14)
Warsawcentroid_data$density_increase <- (Warsawcentroid_data$rowsum18 - Warsawcentroid_data$rowsum_14)

write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")


#######################################################################################################
########################## creating a venue similarity matrix #########################################
########################################################################################################

library(igraph)
library(networkD3)
library(rgexf)
library(tidytree)
library(data.tree)
library(DiagrammeR)

setwd("C:/Dokumente/Utrecht/Master_these/Data")
venue_cat_edges <- read.csv("venue categories_edge list_ids.csv", header = TRUE)
venue_cat_nodes <- read.csv("venue_cat_nodes.csv", header = TRUE)

cat_taxonomy= graph_from_data_frame(d=venue_cat_edges, vertices = venue_cat_nodes, directed = FALSE)
plot.igraph(cat_taxonomy)
venue_cat_degrees<- degree(cat_taxonomy)
write.csv(as.data.frame(venue_cat_degrees), "degrees venue categories.csv", row.names=TRUE)

jpeg("Venue Categories Degrees.jpeg")
plot(density(venue_cat_degrees, from = 0, to= 150, na.rm = TRUE), main = "Foursquare Venue Category Degree centrality")
dev.off()

dist<- distances(cat_taxonomy, v = V(cat_taxonomy), to = V(cat_taxonomy), mode = "all", weights = NULL, algorithm = "automatic")
write.csv(as.data.frame(dist), "distance venue categories.csv", row.names=TRUE)

setwd("C:/Dokumente/Utrecht/Master_these/Data")
venuesimilarity <- read.csv("distance venue categories.csv", header = TRUE)
venuesimilarity$categoryna <- venuesimilarity$X
venuecategories <- read.csv("venuecategories.csv", header = T)
colnames(venuecategories) <- c("categoryna", "categoryid")
venuecategories <- distinct(venuecategories)

venuecategories$categoryna <- as.character(venuecategories$categoryna)
venuesimilarity2 <- merge(venuesimilarity, venuecategories, by = "categoryna", all.x = TRUE)
rowcolumnselection <- venuesimilarity2$rownum[!is.na(venuesimilarity2$categoryid)]
totvenuesimilarity <- venuesimilarity2[!is.na(venuesimilarity2$categoryid), ]
rowcolumnselection <- append(rowcolumnselection, c(-2, 939))
totvenuesimilarity <- totvenuesimilarity[,(rowcolumnselection+3)]

for(i in 1:nrow(totvenuesimilarity)){
  totvenuesimilarity$columnnames[i]<- paste("cat_" , totvenuesimilarity$categoryid[i], sep = "")
}

columnnames <- totvenuesimilarity$columnnames
columnnames <- append(columnnames, c("categoryna", "categoryid", "columnnames"))
colnames(totvenuesimilarity) <- columnnames

write.csv(totvenuesimilarity, "totvenuesimilarity.csv")


#########################################################################################################
########################## Calculating Multifunctionality ###############################################
#########################################################################################################
setwd("C:/Dokumente/Utrecht/Master_these/Data")
totvenuesimilarity <- read.csv("totvenuesimilarity.csv", header=TRUE)
totvenuesimilarity <- subset(totvenuesimilarity, select = -c(X))


entropy <- function(mat) {
  freqs <- mat/rowSums (mat)
  entropy <- - rowSums (freqs * log2(freqs+0.000000001))
  entropy <- round (entropy, digits = 3)
  return (entropy)
}


setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
########### before subway expansion

Warsawvenuecat_matrix_pre2015 <- read.csv("Warsawvenuecat_matrix_pre2015.csv", header = TRUE)
Warsawvenuecat_matrix_pre2015 <- subset(Warsawvenuecat_matrix_pre2015, select = -c(X))
Warsawdissimilardiversity_matrix_2015 <- Warsawvenuecat_matrix_pre2015
Warsawdissimilardiversity_matrix_2015_between <- Warsawvenuecat_matrix_pre2015

Warsawdissimilardiversity_matrix_2015$ncat <- ""
Warsawdissimilardiversity_matrix_2015$entropysum <- ""
Warsawdissimilardiversity_matrix_2015$rowentropy <- ""


colnameset <- colnames(Warsawvenuecat_matrix_pre2015)
for(n in 1:nrow(Warsawvenuecat_matrix_pre2015)){
  localvenuesubset <- c()
  for(m in 1:ncol(Warsawvenuecat_matrix_pre2015)){
    if(Warsawvenuecat_matrix_pre2015[n,m] != 0){
      localvenuesubset <- append(localvenuesubset, m)
    }
  }
  colnameset_local <- colnameset[localvenuesubset]
  localvenuesubset_entries <- Warsawvenuecat_matrix_pre2015[n, localvenuesubset]
  if(length(localvenuesubset) > 0){
    for(g in 1:length(localvenuesubset)){
      a <- which(totvenuesimilarity$columnnames == colnameset_local[g])
      for(p in 1:length(localvenuesubset)){
        s <- which(colnames(totvenuesimilarity) == colnameset_local[p])
        pairsimilarity <- 10/((totvenuesimilarity[a, s]+2)/2)
        Warsawdissimilardiversity_matrix_2015_between[n,localvenuesubset[p]] <-  (localvenuesubset_entries[p] * pairsimilarity)
      }
      Warsawdissimilardiversity_matrix_2015[n, localvenuesubset[g]] <- entropy(Warsawdissimilardiversity_matrix_2015_between[n,])
    }
    Warsawdissimilardiversity_matrix_2015$ncat[n] <- length(localvenuesubset)
    Warsawdissimilardiversity_matrix_2015$entropysum[n] <- sum(Warsawdissimilardiversity_matrix_2015[n, 1:ncol(Warsawvenuecat_matrix_pre2015)])
    Warsawdissimilardiversity_matrix_2015$rowentropy[n] <- as.numeric(Warsawdissimilardiversity_matrix_2015$entropysum[n])/ as.numeric(Warsawdissimilardiversity_matrix_2015$ncat[n])
  }
}

write.csv(Warsawdissimilardiversity_matrix_2015, "Warsawdissimilardiversity_matrix_2015.csv")

########## after subway expansion

Warsawvenuecat_matrix_pre2018 <- read.csv("Warsawvenuecat_matrix_pre2018.csv", header = TRUE)
Warsawvenuecat_matrix_pre2018 <- subset(Warsawvenuecat_matrix_pre2018, select = -c(X))
Warsawdissimilardiversity_matrix_2018 <- Warsawvenuecat_matrix_pre2018
Warsawdissimilardiversity_matrix_2018_between <- Warsawvenuecat_matrix_pre2018


Warsawdissimilardiversity_matrix_2018$ncat <- ""
Warsawdissimilardiversity_matrix_2018$entropysum <- ""
Warsawdissimilardiversity_matrix_2018$rowentropy <- ""


colnameset <- colnames(Warsawvenuecat_matrix_pre2018)
for(n in 1:nrow(Warsawvenuecat_matrix_pre2018)){
  localvenuesubset <- c()
  for(m in 1:ncol(Warsawvenuecat_matrix_pre2018)){
    if(Warsawvenuecat_matrix_pre2018[n,m] != 0){
      localvenuesubset <- append(localvenuesubset, m)
    }
  }
  colnameset_local <- colnameset[localvenuesubset]
  localvenuesubset_entries <- Warsawvenuecat_matrix_pre2018[n, localvenuesubset]
  if(length(localvenuesubset) > 0){
    for(g in 1:length(localvenuesubset)){
      a <- which(totvenuesimilarity$columnnames == colnameset_local[g])
      for(p in 1:length(localvenuesubset)){
        s <- which(colnames(totvenuesimilarity) == colnameset_local[p])
        pairsimilarity <- 10/((totvenuesimilarity[a, s]+2)/2)
        Warsawdissimilardiversity_matrix_2018_between[n,localvenuesubset[p]] <-  (localvenuesubset_entries[p] * pairsimilarity)
      }
      Warsawdissimilardiversity_matrix_2018[n, localvenuesubset[g]] <- entropy(Warsawdissimilardiversity_matrix_2018_between[n,])
    }
    Warsawdissimilardiversity_matrix_2018$ncat[n] <- length(localvenuesubset)
    Warsawdissimilardiversity_matrix_2018$entropysum[n] <- sum(Warsawdissimilardiversity_matrix_2018[n, 1:ncol(Warsawvenuecat_matrix_pre2018)])
    Warsawdissimilardiversity_matrix_2018$rowentropy[n] <- as.numeric(Warsawdissimilardiversity_matrix_2018$entropysum[n])/ as.numeric(Warsawdissimilardiversity_matrix_2018$ncat[n])
  }
}

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
write.csv(Warsawdissimilardiversity_matrix_2018, "Warsawdissimilardiversity_matrix_2018.csv")

Warsawcentroid_data$rowentropy_pre18 <- Warsawdissimilardiversity_matrix_2018$rowentropy[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$rowentropy_pre15 <- Warsawdissimilardiversity_matrix_2015$rowentropy[Warsawcentroid_data$ORIG_FID]
Warsawcentroid_data$rowentropy_pre18[is.na(Warsawcentroid_data$rowentropy_pre18)] <- 0
Warsawcentroid_data$rowentropy_pre15[is.na(Warsawcentroid_data$rowentropy_pre15)] <- 0
Warsawcentroid_data$multifunctionality_diff <- (Warsawcentroid_data$rowentropy_pre18 - Warsawcentroid_data$rowentropy_pre15)

write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")

######################################################################################################################################
########Defining the city center based on maximum amenity density and calculate the points Euclidean distance to center###############
######################################################################################################################################

Warsawcentroid_data$center_x <- Warsawcentroid_data$coords.x1[which(Warsawcentroid_data$rowsum18 == max(na.omit(Warsawcentroid_data$rowsum18)))]
Warsawcentroid_data$center_y <- Warsawcentroid_data$coords.x2[which(Warsawcentroid_data$rowsum18 == max(na.omit(Warsawcentroid_data$rowsum18)))]

Warsawcentroid_data$centroid_dist_center <- spDists(as.matrix(Warsawcentroid_data[, c("coords.x1", "coords.x2")]), as.matrix(Warsawcentroid_data[1,c("center_x", "center_y")]), longlat = T)
Warsawcentroid_data$centroid_dist_center <- Warsawcentroid_data$centroid_dist_center *1000 # in meters



######################################################################################################
################## Rearranging the venuecategorymatrix to typecounts of Metacategory #################
######################################################################################################

setwd("C:/Dokumente/Utrecht/Master_these/Data")
### This is a edge list of metacategories and all their subcategories
metacategory_assignment <- read.csv("metacategory_assignment.csv", header = T)

Food_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Food"),]
Food_colnames <- Food_colnames$columnnames

Shops_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Shop & Service"),]
Shops_colnames <- Shops_colnames$columnnames

ArtsEnter_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Arts & Entertainment"),]
ArtsEnter_colnames <- ArtsEnter_colnames$columnnames

Nightlife_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Nightlife Spot"),]
Nightlife_colnames <- Nightlife_colnames$columnnames

Outdoors_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Outdoors & Recreation"),]
Outdoors_colnames <- Outdoors_colnames$columnnames

Proff_colnames <- metacategory_assignment[which(metacategory_assignment$finalmetacat == "Professional & Other Places"),]
Proff_colnames <- Proff_colnames$columnnames


######## regroup the counts of venuecategories into the Metacategories before and after the new metro stations opened
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")

Warsawvenuecat_matrix_pre2018 <- read.csv("Warsawvenuecat_matrix_pre2018.csv", header = T)
Warsawvenuecat_matrix_pre2015 <- read.csv("Warsawvenuecat_matrix_pre2015.csv", header = T)

Warsawvenuecat_matrix_pre2015$nr_Food_15 <- 0
Warsawvenuecat_matrix_pre2015$nr_ShopsServ_15 <- 0
Warsawvenuecat_matrix_pre2015$nr_ArtsEnter_15 <- 0
Warsawvenuecat_matrix_pre2015$nr_Nightlife_15 <- 0
Warsawvenuecat_matrix_pre2015$nr_Proff_15 <- 0
Warsawvenuecat_matrix_pre2015$nr_Outdoors_15 <- 0

matrix_colnames <- colnames(Warsawvenuecat_matrix_pre2015)
for(i in 1:ncol(Warsawvenuecat_matrix_pre2015)){
  if(is.element(matrix_colnames[i], Food_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_Food_15 <- (Warsawvenuecat_matrix_pre2015$nr_Food_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
  else if(is.element(matrix_colnames[i], Shops_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_ShopsServ_15 <- (Warsawvenuecat_matrix_pre2015$nr_ShopsServ_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
  else if(is.element(matrix_colnames[i], ArtsEnter_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_ArtsEnter_15 <- (Warsawvenuecat_matrix_pre2015$nr_ArtsEnter_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
  else if(is.element(matrix_colnames[i], Nightlife_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_Nightlife_15 <- (Warsawvenuecat_matrix_pre2015$nr_Nightlife_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
  else if(is.element(matrix_colnames[i], Proff_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_Proff_15 <- (Warsawvenuecat_matrix_pre2015$nr_Proff_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
  else if(is.element(matrix_colnames[i], Outdoors_colnames)){
    Warsawvenuecat_matrix_pre2015$nr_Outdoors_15 <- (Warsawvenuecat_matrix_pre2015$nr_Outdoors_15 + Warsawvenuecat_matrix_pre2015[, i])
  }
}


Warsawcentroid_data$nr_Food_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_Food_15 
Warsawcentroid_data$nr_ShopsServ_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_ShopsServ_15
Warsawcentroid_data$nr_ArtsEnter_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_ArtsEnter_15
Warsawcentroid_data$nr_Nightlife_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_Nightlife_15
Warsawcentroid_data$nr_Proff_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_Proff_15
Warsawcentroid_data$nr_Outdoors_15_300mbuff <- Warsawvenuecat_matrix_pre2015$nr_Outdoors_15

Warsawvenuecat_matrix_pre2018$nr_Food18 <- 0
Warsawvenuecat_matrix_pre2018$nr_ShopsServ18 <- 0
Warsawvenuecat_matrix_pre2018$nr_ArtsEnter18 <- 0
Warsawvenuecat_matrix_pre2018$nr_Nightlife18 <- 0
Warsawvenuecat_matrix_pre2018$nr_Proff18 <- 0
Warsawvenuecat_matrix_pre2018$nr_Outdoors18 <- 0

matrix_colnames <- colnames(Warsawvenuecat_matrix_pre2018)
for(i in 1:ncol(Warsawvenuecat_matrix_pre2018)){
  if(is.element(matrix_colnames[i], Food_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_Food18 <- (Warsawvenuecat_matrix_pre2018$nr_Food18 + Warsawvenuecat_matrix_pre2018[, i])
  }
  else if(is.element(matrix_colnames[i], Shops_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_ShopsServ18 <- (Warsawvenuecat_matrix_pre2018$nr_ShopsServ18 + Warsawvenuecat_matrix_pre2018[, i])
  }
  else if(is.element(matrix_colnames[i], ArtsEnter_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_ArtsEnter18 <- (Warsawvenuecat_matrix_pre2018$nr_ArtsEnter18 + Warsawvenuecat_matrix_pre2018[, i])
  }
  else if(is.element(matrix_colnames[i], Nightlife_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_Nightlife18 <- (Warsawvenuecat_matrix_pre2018$nr_Nightlife18 + Warsawvenuecat_matrix_pre2018[, i])
  }
  else if(is.element(matrix_colnames[i], Proff_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_Proff18 <- (Warsawvenuecat_matrix_pre2018$nr_Proff18 + Warsawvenuecat_matrix_pre2018[, i])
  }
  else if(is.element(matrix_colnames[i], Outdoors_colnames)){
    Warsawvenuecat_matrix_pre2018$nr_Outdoors18 <- (Warsawvenuecat_matrix_pre2018$nr_Outdoors18 + Warsawvenuecat_matrix_pre2018[, i])
  }
}


Warsawcentroid_data$nr_Food_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_Food18 
Warsawcentroid_data$nr_ShopsServ_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_ShopsServ18
Warsawcentroid_data$nr_ArtsEnter_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_ArtsEnter18
Warsawcentroid_data$nr_Nightlife_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_Nightlife18
Warsawcentroid_data$nr_Proff_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_Proff18
Warsawcentroid_data$nr_Outdoors_18_300mbuff <- Warsawvenuecat_matrix_pre2018$nr_Outdoors18

Warsawcentroid_data$nr_Food_1518_300mbuff <- Warsawcentroid_data$nr_Food_18_300mbuff - Warsawcentroid_data$nr_Food_15_300mbuff
Warsawcentroid_data$nr_ShopsServ_1518_300mbuff <- Warsawcentroid_data$nr_ShopsServ_18_300mbuff - Warsawcentroid_data$nr_ShopsServ_15_300mbuff 
Warsawcentroid_data$nr_ArtsEnter_1518_300mbuff <- Warsawcentroid_data$nr_ArtsEnter_18_300mbuff - Warsawcentroid_data$nr_ArtsEnter_15_300mbuff
Warsawcentroid_data$nr_Nightlife_1518_300mbuff <- Warsawcentroid_data$nr_Nightlife_18_300mbuff - Warsawcentroid_data$nr_Nightlife_15_300mbuff
Warsawcentroid_data$nr_Proff_1518_300mbuff <- Warsawcentroid_data$nr_Proff_18_300mbuff  - Warsawcentroid_data$nr_Proff_15_300mbuff 
Warsawcentroid_data$nr_Outdoors_1518_300mbuff <- Warsawcentroid_data$nr_Outdoors_18_300mbuff - Warsawcentroid_data$nr_Outdoors_15_300mbuff 



###############################################################################################################
###################### Computing Street network accessability to the center as control variable ###############
###############################################################################################################
setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
Warsaw_Fsq <- Warsaw_Fsq[!is.na(Warsaw_Fsq$lat), ]
kmeans_Warsaw <- as.data.frame(matrix(data = 0, ncol = 2, nrow = 50))
colnames(kmeans_Warsaw) <- c("cluster_number", "tot_withinss")
kmeans_Warsaw$cluster_number <- c(seq(1,50))
k.max <- 50
for(k in 1:k.max){
  kmeans_Warsaw$tot_withinss[k] <- kmeans(cbind(Warsaw_Fsq$lat, Warsaw_Fsq$lon), centers = k )$tot.withinss
}

jpeg("Warsaw k-means Clustering - Optimal Cluster Number.jpeg", width = 8, height = 7, units = "in", res = 1000)
plot(kmeans_Warsaw$cluster_number, kmeans_Warsaw$tot_withinss,
     type="b", pch = 7, frame = FALSE, 
     main = "Warsaw k-means Clustering - Optimal Cluster Number",
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
dev.off()

##select the elbow number
km <- kmeans(cbind(Warsaw_Fsq$lat, Warsaw_Fsq$lon), centers = 11)
jpeg("Warsaw k-means Clustering Cluster-Map.jpeg", width = 8, height = 8, units = "in", res = 1000)
plot(Warsaw_Fsq$lat, Warsaw_Fsq$lon, col = km$cluster, pch = 20, 
    main = "Warsaw k-means Clustering (Optimal Cluster Number = 11)")
dev.off()
nrow(km$centers)
install.packages("concaveman")
library(concaveman)

Warsaw_Fsq$cluster <- km$cluster
pnts <- Warsaw_Fsq[Warsaw_Fsq$cluster == 1,] %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326)
polygon <- concaveman(pnts)
polygon$clusterid <- 1

for(i in 2:nrow(km$centers)){
  pnts <- Warsaw_Fsq[Warsaw_Fsq$cluster == i,] %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
  polygon2 <- concaveman(pnts)
  polygon2$clusterid <- i
  polygon <- rbind(polygon, polygon2)
}

polygon <- as_Spatial(polygon)
proj4string(polygon)=CRS("+proj=longlat +datum=WGS84") # set it to UTM

setwd("C:/Dokumente/Utrecht/Master_these/Data/Warsaw")
Warsawcentroid_data <- read.csv("Warsawcentroid_data.csv", header = TRUE)
coordinates(Warsawcentroid_data)= ~centroid_lat+ centroid_lon
proj4string(Warsawcentroid_data)=CRS("+proj=longlat +datum=WGS84") # set it to UTM
Warsawcentroid_data <- point.in.poly(Warsawcentroid_data, polygon, sp = T, duplicate = F)
subcenters_Warsaw <- as.data.frame(km$centers)
Warsawcentroid_data$rowsum18[is.na(Warsawcentroid_data$rowsum18)] <- 0


subcenter <- Warsawcentroid_data$pid1
Warsawcentroid_data <- read.csv("Warsawcentroid_data.csv", header = TRUE)
Warsawcentroid_data$subcenterid <- subcenter

subcenters <- matrix(data=NA, nrow = nrow(km$centers), ncol = 5)
subcenters <- as.data.frame(subcenters)
colnames(subcenters) <- c("subcenterid", "subcenter_ORIG_FID", "subcenter_lat", "subcenter_lon", "subcenter_amen_density")
for(i in 1:nrow(km$centers)){
  subcentral_subset <-Warsawcentroid_data[Warsawcentroid_data$subcenterid == i,]
  subcentral_subset <- subcentral_subset[!is.na(subcentral_subset$ORIG_FID),]
  subcenters[i,1] <- as.numeric(i)
  subcenters[i,2] <- subcentral_subset$ORIG_FID[which(subcentral_subset$rowsum15 == max(subcentral_subset$rowsum15))][1]
  subcenters[i,3] <- subcentral_subset$centroid_lat[subcentral_subset$ORIG_FID == subcenters[i,2]]
  subcenters[i,4] <- subcentral_subset$centroid_lon[subcentral_subset$ORIG_FID == subcenters[i,2]]
  subcenters[i,5] <- subcentral_subset$rowsum15[subcentral_subset$ORIG_FID == subcenters[i,2]]
}


Warsaw_centroids <- Warsawcentroid_data[,c( "centroid_lon", "centroid_lat")]
rownames(Warsaw_centroids) <- c(Warsawcentroid_data$ORIG_FID)
colnames(Warsaw_centroids) <- c("lon", "lat")
subcenters_Warsaw <- subcenters[,c("subcenter_lon", "subcenter_lat")]
colnames(subcenters_Warsaw) <- c("lon", "lat")
dis_matrix <- pointDistance( Warsaw_centroids, subcenters_Warsaw, lonlat=TRUE)
dis_matrix <- as.data.frame(dis_matrix)
colnames(dis_matrix)<-row.names(subcenters_Warsaw)
row.names(dis_matrix)<-row.names(Warsaw_centroids)
for(i in 1:nrow(dis_matrix)){
  dis_matrix$subcenterID[i] <- which(dis_matrix[i, 1:nrow(km$centers)] == min(dis_matrix[i, 1:nrow(km$centers)]))
}
Warsawcentroid_data$subcenterID <- dis_matrix$subcenterID
subcenters_Warsaw$subcenterID <- c(seq(1,nrow(subcenters_Warsaw)))
colnames(subcenters_Warsaw) <- c("subcenter_lon", "subcenter_lat", "subcenterID")
Warsawcentroid_data <- merge(Warsawcentroid_data, subcenters_Warsaw, by = "subcenterID", all.x = T)


write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")

###traveltime	tt	Indicates whether the travel time information should be provided in summary entries
###costfactor	cf	Indicates whether the CostFactor information should be returned in summary entries.
###distance	di	Indicates whether distance information should be returned in summary entries.
###walking traveltime and distance are equivalent, probably based on the assumption that it takes 1 second to walk 1 meter

Warsawaccessability_matrix <- subset(Warsawcentroid_data, select = c(ORIG_FID, centroid_lon, centroid_lat, subcenter_lon, subcenter_lat))
names(Warsawaccessability_matrix) <- c("ORIG_FID", "centroid_lat", "centroid_lon", "center_x", "center_y")

Warsawaccessability_matrix$traveltime <- ""
Warsawaccessability_matrix$streetdistance <- ""

Warsawaccessability_matrix$centroid_lat <- as.numeric(Warsawaccessability_matrix$centroid_lat)
Warsawaccessability_matrix$center_y <- as.numeric(Warsawaccessability_matrix$center_y)



apikey <- 'ZFoIDwMr_B2_gQtL9ojl89gCuOhSi1tZA8Af3sUOT1Q'
apikey <- 'rdpNQnG-8sphTD7RkU154p4Zs6TOGEulfr0ap0m5jMA'
apikey <- '3w3-gzPjO4fkeQYA94-hKMaj9zEEDBZCLKZD_diWcT4'

TransportMode <- 'bicycle'
TransportMode <- 'pedestrian'

TransportMode <- 'car'
for(m in 1:nrow(Warsawaccessability_matrix)){
  url_geocod <-  paste("https://matrix.route.ls.hereapi.com/routing/7.2/calculatematrix.xml?apiKey=", as.character(apikey),
                       "&start0=geo!", as.character(Warsawaccessability_matrix$centroid_lat[m]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m]),
                       "&start1=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+1]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+1]), 
                       "&start2=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+2]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+2]), 
                       "&start3=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+3]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+3]), 
                       "&start4=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+4]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+4]), 
                       "&start5=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+5]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+5]), 
                       "&start6=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+6]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+6]), 
                       "&start7=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+7]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+7]), 
                       "&start8=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+8]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+8]), 
                       "&start9=geo!",  as.character(Warsawaccessability_matrix$centroid_lat[m+9]), ",", as.character(Warsawaccessability_matrix$centroid_lon[m+9]), 
                       "&destination0=geo!",  as.character(Warsawaccessability_matrix$center_x[m]), ",", as.character(Warsawaccessability_matrix$center_y[m]), 
                       "&mode=fastest;", as.character(TransportMode),";traffic:default",
                       "&summaryAttributes=distance,traveltime",
                       sep = "")
  webpage <- read_html(as.character(url_geocod))
  Warsawaccessability_matrix$traveltime[m] <- webpage[9]
  Warsawaccessability_matrix$streetdistance[m] <-webpage[8]
  Warsawaccessability_matrix$traveltime[m+1] <- webpage[14]
  Warsawaccessability_matrix$streetdistance[m+1] <-webpage[13]
  Warsawaccessability_matrix$traveltime[m+2] <- webpage[19]
  Warsawaccessability_matrix$streetdistance[m+2] <-webpage[18]
  Warsawaccessability_matrix$traveltime[m+3] <- webpage[24]
  Warsawaccessability_matrix$streetdistance[m+3] <-webpage[23]
  Warsawaccessability_matrix$traveltime[m+4] <- webpage[29]
  Warsawaccessability_matrix$streetdistance[m+4] <-webpage[28]
  Warsawaccessability_matrix$traveltime[m+5] <- webpage[34]
  Warsawaccessability_matrix$streetdistance[m+5] <-webpage[33]
  Warsawaccessability_matrix$traveltime[m+6] <- webpage[39]
  Warsawaccessability_matrix$streetdistance[m+6] <-webpage[38]
  Warsawaccessability_matrix$traveltime[m+7] <- webpage[44]
  Warsawaccessability_matrix$streetdistance[m+7] <-webpage[43]
  Warsawaccessability_matrix$traveltime[m+8] <- webpage[49]
  Warsawaccessability_matrix$streetdistance[m+8] <-webpage[48]
  Warsawaccessability_matrix$traveltime[m+9] <- webpage[54]
  Warsawaccessability_matrix$streetdistance[m+9] <-webpage[53]
  m <- m+9
}

write.csv(Warsawaccessability_matrix, "Warsawaccessability_matrix.csv")

colnames(subcenters) <-  c("subcenterID" ,"subcenter_ORIG_FID", "subcenter_lat", "subcenter_lon", "subcenter_amen_density")
subcenters <- subcenters[,c("subcenterID", "subcenter_amen_density")]
Warsawcentroid_data <- merge(Warsawcentroid_data, subcenters, by = "subcenterID")
Warsawcentroid_data <- merge(Warsawcentroid_data, Warsawaccessability_matrix, by = "ORIG_FID", all.x = T)

write.csv(Warsawcentroid_data, "Warsawcentroid_data.csv")
