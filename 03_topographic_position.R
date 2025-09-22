library("terra")
library("sf")
library("landscapemetrics")
library("viridis")
library("tidyverse")
library("patchwork")
library("scales")
library("exactextractr")


set.seed(406)


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


fire_patches <- get_patches(severe_rast, directions = 8)
fire_patches <- rast(fire_patches$layer_1$class_1)

severe_poly <- as.polygons(fire_patches, aggregate = TRUE)
crs(severe_poly) <- crs(rast2016)
severe_poly <- st_as_sf(severe_poly)
severe_poly$ID <- c(1:nrow(severe_poly))
names(severe_poly) <- c("geometry", "ID")
head(severe_poly)


loghist <- st_read("spatial/LASTLOG25.shp")
loghist <- loghist %>% 
  mutate(year = as.numeric(substr(SEASON, 1 , 4))+1, 
         area = as.numeric(st_area(loghist)/10000),
         history = "harvest")
logarea <- loghist %>%
  filter(year >= 1970, X_FORETYPE == "Mountain Ash",
         X_SILVSYS == "Clearfelling" |
           X_SILVSYS == "Seed Tree (includes retained overwood)")
logarea <- logarea[,26:27]
logarea <- vect(logarea)
logarea <- project(logarea, rast2016)
logarea <- crop(logarea, rast2016)
logarea <- st_as_sf(logarea)
logarea$ID <- 1:nrow(logarea)



DEM <- rast("spatial/DEM.grd")
northness <- cos(DEM[["aspect"]] * pi / 180)
eastness <- sin(DEM[["aspect"]] * pi / 180)
twi <- rast("spatial/twi.grd")

elev <- exact_extract(DEM, severe_poly, fun = "mean", append_cols = "ID")

east <- exact_extract(eastness, severe_poly, fun = "mean", append_cols = "ID")
names(east) <- c("ID", "mean.eastness")

north <- exact_extract(northness, severe_poly, fun = "mean", append_cols = "ID")
names(north) <- c("ID", "mean.northness")

twi_fire <- exact_extract(twi, severe_poly, fun = "mean", append_cols = "ID")
names(twi_fire) <- c("ID", "mean.twi")

topo_fire <- left_join(elev, east)
topo_fire <- left_join(topo_fire, north)
topo_fire <- left_join(topo_fire, twi_fire)
topo_fire$history <- "fire"

write.csv(topo_fire, "data/fire_patches_topo.csv")


elev <- exact_extract(DEM, logarea, fun = "mean", append_cols = "ID")

east <- exact_extract(eastness, logarea, fun = "mean", append_cols = "ID")
names(east) <- c("ID", "mean.eastness")

north <- exact_extract(northness, logarea, fun = "mean", append_cols = "ID")
names(north) <- c("ID", "mean.northness")

twi_harv <- exact_extract(twi, logarea, fun = "mean", append_cols = "ID")
names(twi_harv) <- c("ID", "mean.twi")


topo_harv <- left_join(elev, east)
topo_harv <- left_join(topo_harv, north)
topo_harv <- left_join(topo_harv, twi_harv)
topo_harv$history <- "harvest"

write.csv(topo_harv, "data/harv_patches_topo.csv")


# Read files --------------------------------------------------------------

topo_fire <- read.csv("data/fire_patches_topo.csv")
topo_harv <- read.csv("data/harv_patches_topo.csv")

# Plot it out! ------------------------------------------------------------
set.seed(406)

topo_fire <- na.omit(topo_fire)

topo_sample <- topo_fire %>%
  slice_sample(n = 3450)
all_data <- rbind(topo_sample, topo_harv)


a <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.elevation, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C")  +
  scale_colour_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean elevation", y = "Proportion")



b <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.slope, color = history, fill = history), alpha = 0.5, position = "identity", show.legend =  FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean slope", y = "Proportion")



c <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.eastness, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean eastness", y = "Proportion")


d <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.northness, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C")  +
  scale_colour_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean northness", y = "Proportion")

e <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.twi, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C")  +
  scale_colour_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean TWI", y = "Proportion")


legend_plot <- ggplot(all_data)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=mean.twi, color = history, fill = history), alpha = 0.5, position = "identity")+
  scale_fill_viridis(discrete = TRUE, option = "C")  +
  scale_colour_viridis(discrete = TRUE, option = "C")  +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean TWI", y = "Proportion")

legend <- cowplot::get_legend(legend_plot)

layout <- "
AABBCC
DDEEF#
"

a + b + c + d + e + legend + plot_layout(design = layout) + theme(legend.position = "none")
