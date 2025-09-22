library("terra")
library("sf")
library("landscapemetrics")
library("viridis")
library("tidyverse")



# Data setup --------------------------------------------------------------


rast2007 <- rast("spatial/predicted2007.grd")
rast2016 <- rast("spatial/predicted2016.grd")


fire_boundary <- rast("spatial/FIRE_SEV09.grd")
fire_boundary <- crop(fire_boundary, rast2007)

rast2007 <- mask(rast2007, fire_boundary)
rast2016 <- mask(rast2016, fire_boundary)

df <- c(rast2007, rast2016, fire_boundary)
df <- as.data.frame(df, xy = TRUE)
df$transition <- paste0(df$str2007, " to ", df$str2016)

severe_rast <- df %>%
  filter(!is.na(str2007), !is.na(CLASS)) %>%
  mutate(code = case_when(transition == "1 to 1" ~ NA,
                          transition == "1 to 2" ~ NA,
                          transition == "1 to 3" ~ 1,
                          transition == "2 to 1" ~ NA,
                          transition == "2 to 2" ~ NA,
                          transition == "2 to 3" ~ 1))
severe_rast <- rast(severe_rast[c(1,2,7)], type = "xyz")

transition_rast <- df %>%
  filter(!is.na(str2007), !is.na(CLASS)) %>%
  mutate(code = case_when(transition == "1 to 1" ~ 11,
                          transition == "1 to 2" ~ 12,
                          transition == "1 to 3" ~ 13,
                          transition == "2 to 1" ~ 21,
                          transition == "2 to 2" ~ 22,
                          transition == "2 to 3" ~ 23))
transition_rast <- rast(transition_rast[c(1,2,7)], type = "xyz")


binary_rast <- df %>%
  filter(!is.na(str2007), !is.na(CLASS)) %>%
  mutate(code = case_when(transition == "1 to 1" ~ 0,
                          transition == "1 to 2" ~ 1,
                          transition == "1 to 3" ~ 1,
                          transition == "2 to 1" ~ 0,
                          transition == "2 to 2" ~ 0,
                          transition == "2 to 3" ~ 1))
binary_rast <- rast(binary_rast[c(1,2,7)], type = "xyz")

# 1 ha crown loss ---------------------------------------------------------


#Create a raster with canopy cover by value
canopy2007 <- mask(rast2007, fire_boundary)
canopy2007[canopy2007 ==2] <- .5
canopy2007[is.na(canopy2007)] <- 0

canopy2016 <- mask(rast2016, fire_boundary)
canopy2016[canopy2016 ==2] <- .5
canopy2016[canopy2016 ==3] <- 0
canopy2016[is.na(canopy2016)] <- 0

#Aggregate to get canopy cover
canopy2007 <- aggregate(canopy2007, fact = 7, fun = sum)
canopy2007[canopy2007 == 0] <- NA

canopy2016 <- aggregate(canopy2016, fact = 7, fun = sum)
canopy2016 <- mask(canopy2016, canopy2007)

#Max canopy cover
coverage <- rast2007
coverage[coverage >0] <- 1
coverage <- aggregate(coverage, fact = 7, fun = sum, na.rm = TRUE)

#Divide by maximum canopy coverage per aggregated pixel
canopy2007 <- canopy2007/coverage
canopy2016 <- canopy2016/coverage

plot(canopy2007)
plot(canopy2016)


canopyloss <- canopy2016-canopy2007
canopyloss[canopyloss > 0] <- 0
plot(canopyloss)


canopyloss <- as.data.frame(canopyloss, xy = TRUE)

canopyloss <- canopyloss %>%
  mutate(treatment = case_when(str2016 > -0.3 ~ 1, #No harvesting
                               str2016 <= -0.3 & str2016 > -0.5 ~ 2, #thinning
                               str2016 <= -0.5 & str2016 > -0.8 ~3, #retention
                               str2016 <= -0.8 ~ 4)) #CF/seed tree

canopyzones <- rast(canopyloss, type = "xyz")

plot(canopyzones[[2]])

zone1 <- canopyzones[[2]]
zone2 <- canopyzones[[2]]
zone3 <- canopyzones[[2]]
zone4 <- canopyzones[[2]]


zone1[zone1!=1] <- NA
zone2[zone2!=2] <- NA
zone3[zone3!=3] <- NA
zone4[zone4!=4] <- NA

# Agg vs dispersed retention ----------------------------------------------

zone2_patches <- patches(zone2, directions = 4)
zone2_sf <- as.polygons(zone2_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone2_sf <- st_as_sf(zone2_sf)
zone2_sf$area <- st_area(zone2_sf)/10000

zone3_patches <- patches(zone3, directions = 4)
zone3_sf <- as.polygons(zone3_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone3_sf <- st_as_sf(zone3_sf)
zone3_sf$area <- st_area(zone3_sf)/10000

zone4_patches <- patches(zone4, directions = 4)
zone4_sf <- as.polygons(zone4_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone4_sf <- st_as_sf(zone4_sf)
zone4_sf$area <- st_area(zone4_sf)/10000


# get areas of each -------------------------------------------------------

tab <- freq(canopyzones[[2]])
total <- sum(tab$count)
tab$count/total



# Agg or disp? ------------------------------------------------------------

zone3_poly <- as.polygons(zone3, aggregate = FALSE)
zone3_poly <- st_as_sf(zone3_poly)
zone3_poly$ID <- 1:nrow(zone3_poly)


all_rasters <- list()

for(i in 1:nrow(zone3_poly)){
  
  local <- zone3_poly[i,]
  landscape <- crop(binary_rast, local)
  
  all_rasters[[i]] <- landscape
}

lsm_list2 <- list()

for(i in 1:nrow(zone3_poly)){
  
  current_rast <- all_rasters[[i]]
  
  lsm <- calculate_lsm(current_rast, #input raster
                       what = c("lsm_c_ca",
                                "lsm_c_lpi",
                                "lsm_c_np"), #mean euclidean nearest neighbour), 
                       directions = 4) #This is the neighborhood rule we are using
  
  lsm$row <- i
  lsm_list2[[i]] <- lsm
  
}


lsm_table <- do.call(rbind, lsm_list2)
lsm_table1 <- lsm_table %>%
  filter(class == 1)

ca <- lsm_table1 %>%
  filter(metric == "ca")

lpi <- lsm_table1 %>%
  filter(metric == "lpi")

np <- lsm_table1 %>%
  filter(metric == "np")

combined <- cbind(ca[,6], lpi[,6], np[,6:7])
names(combined) <- c("ca", "lpi", "np", "ID")

agg <- combined %>%
  filter(np == 1) %>%
  mutate(type = "agg")

disp <- combined %>%
  filter(np > 1) %>%
  mutate(type = "disp")


summary(agg$ca)
summary(disp$ca)

new_col <- rbind(disp[,4:5], agg[,4:5])



lsm_table2 <- do.call(rbind, lsm_list2)
lsm_table2 <- lsm_table2 %>%
  filter(class == 1)

ca <- lsm_table2 %>%
  filter(metric == "ca")

lpi <- lsm_table2 %>%
  filter(metric == "lpi")

np <- lsm_table2 %>%
  filter(metric == "np")

combined <- cbind(ca[,6], lpi[,6], np[,6:7])
names(combined) <- c("ca", "lpi", "np", "ID")



agg <- combined %>%
  filter(np == 1) %>%
  mutate(type = 5)

disp <- combined %>%
  filter(np > 1) %>%
  mutate(type = 6)


summary(agg$ca)
summary(disp$ca)

new_col <- rbind(disp[,4:5], agg[,4:5])
names(new_col) <- c("patches", "type")

new_zone3 <- left_join(zone3_sf, new_col)
new_zone3_rast <- rasterize(new_zone3, zone3, field = "type")


agg_rast <- new_zone3_rast
agg_rast[agg_rast == 6] <- NA

disp_rast <- new_zone3_rast
disp_rast[disp_rast == 5] <- NA


agg_poly <- patches(agg_rast, directions = 4)
agg_poly <- as.polygons(agg_poly, aggregate = TRUE, values = TRUE, na.rm = TRUE)
agg_poly <- st_as_sf(agg_poly)
agg_poly$area <- st_area(agg_poly)/10000


disp_poly <- patches(disp_rast, directions = 4)
disp_poly <- as.polygons(disp_poly, aggregate = TRUE, values = TRUE, na.rm = TRUE)
disp_poly <- st_as_sf(disp_poly)
disp_poly$area <- st_area(disp_poly)/10000

new_class <- merge(zone1, zone2, zone4, agg_rast, disp_rast)
new_class[new_class == 5] <- 3
new_class[new_class == 4] <- 5
new_class[new_class == 6] <- 4

plot(new_class)


#Proportions
tab <- freq(new_class)
sum(tab$count) - 15203
tab$count/11977

writeRaster(new_class, "spatial/silv_equiv.grd")

# Topography --------------------------------------------------------------

DEM <- rast("spatial/DEM.grd")
northness <- cos(DEM[["aspect"]] * pi / 180)
eastness <- sin(DEM[["aspect"]] * pi / 180)
twi <- rast("spatial/twi.grd")
twi <- project(twi, DEM)

crs(new_class) <- crs(canopy2016)


stack <- c(DEM, northness, eastness, twi)
stack <- aggregate(stack, fact = 7, fun = "mean")

canopyzones <- resample(new_class, stack, method = "near")

stack <- mask(stack, canopyzones)
stack <- c(canopyzones, stack)


zones_df <- as.data.frame(stack, xy = TRUE)
names(zones_df) <- c("x", "y", "treatment", "elevation", "slope", "aspect", "tpi", "tri", "nothness", "eastness", "twi")

ggplot(zones_df, aes(x = factor(treatment), y = elevation)) +
  geom_boxplot()
ggplot(zones_df, aes(x = factor(treatment), y = slope)) +
  geom_boxplot()
ggplot(zones_df, aes(x = factor(treatment), y = nothness)) +
  geom_boxplot()
ggplot(zones_df, aes(x = factor(treatment), y = eastness)) +
  geom_boxplot()
ggplot(zones_df, aes(x = factor(treatment), y = twi)) +
  geom_boxplot()



# Patch sizes etc ---------------------------------------------------------

mets <- list_lsm(level = "class")

lsm <- calculate_lsm(new_class, #input raster
                     what = c("lsm_c_area_mn",
                              "lsm_c_area_sd"), #mean euclidean nearest neighbour), 
                     directions = 4)


zone2 <- new_class
zone3 <- new_class
zone4 <- new_class
zone5 <- new_class


zone2[zone2!=2] <- NA
zone3[zone3!=3] <- NA
zone4[zone4!=4] <- NA
zone5[zone5!=5] <- NA

zone2_patches <- patches(zone2, directions = 4)
zone2_sf <- as.polygons(zone2_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone2_sf <- st_as_sf(zone2_sf)
zone2_sf$area <- st_area(zone2_sf)/10000

zone3_patches <- patches(zone3, directions = 4)
zone3_sf <- as.polygons(zone3_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone3_sf <- st_as_sf(zone3_sf)
zone3_sf$area <- st_area(zone3_sf)/10000

zone4_patches <- patches(zone4, directions = 4)
zone4_sf <- as.polygons(zone4_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone4_sf <- st_as_sf(zone4_sf)
zone4_sf$area <- st_area(zone4_sf)/10000

zone5_patches <- patches(zone5, directions = 4)
zone5_sf <- as.polygons(zone5_patches, aggregate = TRUE, values = TRUE, na.rm = TRUE)
zone5_sf <- st_as_sf(zone5_sf)
zone5_sf$area <- st_area(zone5_sf)/10000


np_zone2 <- nrow(zone2_sf)
np_zone3 <- nrow(zone3_sf)
np_zone4 <- nrow(zone4_sf)
np_zone5 <- nrow(zone5_sf)

min_zone2 <- min(zone2_sf$area)
min_zone3 <- min(zone3_sf$area)
min_zone4 <- min(zone4_sf$area)
min_zone5 <- min(zone5_sf$area)

max_zone2 <- max(zone2_sf$area)
max_zone3 <- max(zone3_sf$area)
max_zone4 <- max(zone4_sf$area)
max_zone5 <- max(zone5_sf$area)

mn_zone2 <- mean(zone2_sf$area)
mn_zone3 <- mean(zone3_sf$area)
mn_zone4 <- mean(zone4_sf$area)
mn_zone5 <- mean(zone5_sf$area)

med_zone2 <- median(zone2_sf$area)
med_zone3 <- median(zone3_sf$area)
med_zone4 <- median(zone4_sf$area)
med_zone5 <- median(zone5_sf$area)

np <- rbind(np_zone2, np_zone3, np_zone4, np_zone5)
min <- rbind(min_zone2, min_zone3, min_zone4, min_zone5)
max <- rbind(max_zone2, max_zone3, max_zone4, max_zone5)
mn <- rbind(mn_zone2, mn_zone3, mn_zone4, mn_zone5)
med <- rbind(med_zone2, med_zone3, med_zone4, med_zone5)

tab2 <- as.data.frame(cbind(np,min,max,mn,med))
names(tab2) <- c("np", "min", "max", "mn", "md")
