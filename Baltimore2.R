# Moran -> Outliers (per year?)
# Circular cartogam
# Take as a model on a network, define some neighbourhood structure, test correlation, etc.
# Permutation test on a simple mpp
# Kernel density?

# Limitation: Distance is not Eulerian



# library(dismo) # Better than ggmap for proper CRS --- yeah but now doesn't work because of Google API changes
library(dplyr)
library(tidyr)
library(ggmap)
# library(data.table)

library(sf)
library(tmap)
library(ggplot2)
library(raster)

library(spatstat)


######### Loading the data ##########

getwd()
setwd("C:/gibbs/Baltimore")


## Get map of baltimore - need to use Google API
# register_google(key = "AIzaSyBHr6A2shDMVd8w[]DOcwBDYnmmaNEuKqBRg")
coords <- geocode("baltimore")
map_raw <- get_map(location = "baltimore", zoom = 11,
                   maptype = "satellite", source = "google")
ggmap(map_raw)
map <- map_raw
#map <- gmap('Baltimore', type='terrain',zoom=11) # gmap from dismo doesn't seem to work, probably because of the change in google API


## Convert ggmap object to true RasterBrick
m <- as.matrix(map)
mrgb <- col2rgb(as.vector(m))
mr <- matrix(mrgb[1,], ncol = ncol(m), nrow = nrow(m))
mg <- matrix(mrgb[2,], ncol = ncol(m), nrow = nrow(m))
mb <- matrix(mrgb[3,], ncol = ncol(m), nrow = nrow(m))
map <- brick(raster(mr),raster(mg),raster(mb))
rm(m,mrgb,mr,mg,mb)

projection(map) <- CRS("+init=epsg:3857")  # The grid is projected in pseudomercator
bb <- as(extent(unlist(attr(map_raw,"bb"))[c(2,4,1,3)]),"SpatialPolygons")
projection(bb) <- CRS("+init=epsg:4326") # While the bounding box is in 4326
bb <- spTransform(bb, CRS("+init=epsg:3857"))@bbox
extent(map) <- c(bb[1,],bb[2,])

crs(map)


## Read the Baltimore City Line shapefile
balt_city <- st_read('cityshp/geo_export_a293f06e-391f-4b66-bd60-94813f31fef3.shp')
st_crs(balt_city)
balt_city <- st_transform(balt_city, crs = crs(map, asText=TRUE))

# Crop map to Baltimore City
map <- crop(map,as(st_buffer(balt_city,2000),"Spatial" ))

tm_shape(map)+
  tm_rgb()+
  tm_shape(balt_city)+
  tm_lines()
  
  
## Read the neighbourhood structure 
neighborhoods <- st_read('neighborhoods/geo_export_a7cfb788-4ca3-4a43-8949-2529800d0a9f.shp')
neighborhoods <- neighborhoods %>% dplyr::select(population,name)
st_crs(neighborhoods)
neighborhoods <- st_transform(neighborhoods, crs = crs(map, asText=TRUE))

tm_shape(map)+
 tm_rgb()+
 tm_shape(neighborhoods)+
 tm_polygons(alpha = 0.5)+
 tm_shape(balt_city)+
 tm_lines(col='red')


## Load the flat data
balt <- read.csv('pair_correlation_dataset.csv', header= TRUE, sep=",", check.names = FALSE)

# Convert to sf file
balt <- st_as_sf(balt, coords = c("longitude","latitude"), crs = 4326)
balt <- st_transform(balt, crs = crs(map, asText=TRUE))

# 0 - Inhabited
# 1 - Vacant
head(balt)
summary(balt)


# Plot everything to see if it fits
tm_shape(map)+
  tm_rgb()+
tm_shape(balt_city)+
  tm_lines(col='red', lwd=2)+
tm_shape(neighborhoods)+
  tm_polygons()+
tm_shape(balt)+
  tm_dots()





######### Data Manipulation #########
# Duplicates?
anyDuplicated(st_coordinates(balt)) # No two observations share a locaion

## Use only one neighborhoods for simplicity and speed
neighborhood <- neighborhoods %>% filter(name == "Barclay")
balt_EW <- balt[st_intersects(neighborhood, balt, sparse=TRUE)[[1]], ]



# Tidy the data
balt_EW.tidy <- balt_EW %>% gather(Year, Vacant, `2004`:`2015`)
# Take only positive
balt_EW.vac <- balt_EW.tidy %>% filter(Vacant == 1) %>% dplyr::select(-Vacant)

# Marked pp with totals
cat(paste(2004:2015,sep=""),sep="`+`")
cat(paste(2004:2015,sep=""),sep="`,`")
balt_EW.total <- balt_EW %>% 
  mutate(total = `2004`+`2005`+`2006`+`2007`+`2008`+`2009`+`2010`+`2011`+`2012`+`2013`+`2014`+`2015`) %>% 
  dplyr::select(-c(`2004`,`2005`,`2006`,`2007`,`2008`,`2009`,`2010`,`2011`,`2012`,`2013`,`2014`,`2015`))



# Convert to spatstat's ppp for later
neighborhood.spat <- as(neighborhood, Class="Spatial")
neighborhood.lines <- neighborhood.spat@polygons[[1]]@Polygons[[1]]@coords

#balt_city.spat <- as(balt_city, Class="Spatial")
#balt_city.lines <- balt_city.spat@lines[[1]]@Lines[[1]]@coords

coords <- data.frame(st_coordinates(balt_EW.total))
balt_EW.total.ppp <- ppp(coords$X, coords$Y, poly=list(x=rev(neighborhood.lines[,1]), y=rev(neighborhood.lines[,2])), marks = balt_EW.total$total)


coords <- data.frame(st_coordinates(balt_EW.vac))
balt_EW.ppp <- ppp(coords$X, coords$Y, poly=list(x=rev(neighborhood.lines[,1]), y=rev(neighborhood.lines[,2])), marks = balt_EW.vac$Year)



# Output standardized data for splancs
standardize <- function(v, pad = 0.01){
  return( ( (v - min(v)) / (max(v) - min(v)) )*(1-pad) + pad/2    )
}
coords.std <- apply(coords,2,standardize)
times.std <- standardize(as.numeric(balt_EW.vac$Year))

saveRDS(list(coords = coords.std, times= times.std), "balt_standardized.rds")



######### Exploratory Analysis #########
# Sums per year
as.data.frame(balt_EW) %>% dplyr::select(-one_of(c('geometry','total'))) %>% summarise_all(sum)

# Map of vacant by year 
tm_shape(neighborhood, bbox = st_bbox(neighborhood)) + tm_polygons(alpha=0.3) + tm_shape(balt_EW.vac) + tm_dots() + tm_facets("Year")


# Histogram of totals over all years
ggplot(data=balt_EW.total, aes(total)) + 
  geom_histogram() 
# Total > 0
balt_EW.total %>% 
  filter(total>0) %>% 
  ggplot(.,aes(total)) + 
  geom_histogram(bins=10) 


# Histogram of vacant per year
balt_EW.vac %>% ggplot(.,aes(Year)) +
  geom_bar()



## Vacant runs
maxRun <- function(x){
  runs <- rle(x)
  if (1 %in% runs$values) {
    return(with(runs, max(lengths[values==1])))
  }
  else {
    return(0)
  }
}


st_dropgeom <- function(df) {
  as.data.frame(df) %>% 
    dplyr::select(-geometry)  
}


balt_EW.nogeom <- st_dropgeom(balt_EW)
runs <- apply(balt_EW.nogeom, 1, maxRun)
hist(runs)
hist(runs[runs>0])


## Plot of totals
balt_EW.total$ftotal = cut(balt_EW.total$total, breaks=c(0,1,2,4,6,8,10,12), include.lowest=TRUE)
darkblue = "#0a1628"
tmap_options(bg.color=darkblue, legend.text.color = "white", inner.margins = c(0.05,0.1,0.05,0.05))
tm_shape(neighborhood) + 
  tm_polygons(col=darkblue) + 
  tm_shape(balt_EW.total) + 
  tm_dots(col="ftotal", scale=2, palette="OrRd") 
#  +tm_layout(bg.color="black", legend.text.color = "white", inner.margins = c(0.05,0.1,0.05,0.05))




######## Testing ###########
# Permutation test on totals
# Correlation for some graph structure
# Diggle surface








## Work with splancs
library(splancs)


spl.xy <- as.matrix(coords(ppsmall))
spl.t <- as.numeric(as.character(ppsmall$marks))
spl.win <- as.matrix(as.data.frame(Window(ppsmall)))
tlimits <- c(2004,2015)
s <- seq(0.001,0.15,0.005)
tm <- seq(2004,2015,1)

khat <- stkhat(spl.xy,spl.t,spl.win,tlimits,s,tm)
se <- stsecal(spl.xy, spl.t, spl.win, tlimits, s, tm)
stdiagn(spl.xy,khat,se, Dzero=TRUE)

mc <- stmctest(spl.xy,spl.t,spl.win,tlimits,s,tm,10)
stdiagn(spl.xy,khat,se,mc)











######### Graveyard #########
## Join neighborhoodss
joined <- st_join(neighborhoods, balt.tidy, left=TRUE)
neigh_all <- joined %>% group_by(name,Year) %>% summarize(Homes = n(), Vacant = sum(Vacant))


aggregate(joined, Vacant ~ name + Year, FUN="sum" )


# Random subset to make EDA faster
balt_s <- balt %>% sample_frac(0.1)
balt_s.tidy <- balt_s %>% gather(Year, Vacant, `2004`:`2015`)
