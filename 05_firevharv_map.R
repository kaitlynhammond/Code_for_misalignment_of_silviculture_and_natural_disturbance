library(terra)
library(sf)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggspatial)
library(patchwork)
library(ggnewscale)
library(spatialEco)


# Data --------------------------------------------------------------------
bbox <- st_read("spatial/bbox.shp")
mtnash <- rast("spatial/mtnash_extent.grd") 
silv_equiv <- rast("Spatial/silv_equiv.grd")
crs(silv_equiv) <- crs(mtnash)
bbox <- st_as_sfc(st_bbox(bbox, crs = st_crs(mtnash)))


silv_equiv <- as.polygons(silv_equiv, dissolve=TRUE)
silv_equiv <- st_as_sf(silv_equiv)
silv_equiv <- st_transform(silv_equiv, st_crs(bbox))
silv_equiv <- st_crop(silv_equiv, bbox)
silv_equiv <- silv_equiv %>%
  mutate(System = case_when(treatment == 1 ~ "Unharvested",
                            treatment == 2 ~ "Selective harvesting",
                            treatment == 3 ~ "Dispersed retention",
                            treatment == 4 ~ "Aggregated retention",
                            treatment == 5 ~ "Clearfell/Seed tree"),
         System = factor(System, levels = c("Unharvested", "Selective harvesting","Dispersed retention", "Aggregated retention", "Clearfell/Seed tree")))



firesev <- vect("Spatial/FIRE_SEV09_POLY.shp")
firesev <- project(firesev, crs(mtnash))
firesev <- st_as_sf(firesev) %>%
  mutate(fire_class = factor(case_when(CLASS == 1 ~ "Crown Burn",
                                       CLASS == 2 ~ "Crown Scorch",
                                       CLASS == 3 ~ "Moderate Crown Scorch",
                                       CLASS == 4 ~ "Understorey Burnt",
                                       CLASS == 5 ~ "Understorey Burnt"), 
                             levels = c("Crown Burn","Crown Scorch", "Moderate Crown Scorch", "Understorey Burnt"))) %>%
  filter(!is.na(fire_class))
firesev <- st_transform(firesev, st_crs(bbox))
firesev <- st_crop(firesev, bbox)


lidar2015 <- vect("Spatial/Lidar2015_boundary_name.shp")
lidar2015 <- project(lidar2015, crs(mtnash))
lidar2015 <- rasterize(lidar2015, mtnash)
lidar2015 <- as.polygons(lidar2015, dissolve=TRUE)
lidar2015 <- st_as_sf(lidar2015)

firesev <- st_intersection(firesev, silv_equiv)

australia <- read_sf("Spatial/STE_2021_AUST_GDA2020.shp")
australia <- australia %>%
  mutate(is_vic = ifelse(STE_NAME21 == "Victoria", 1, 0)) %>%
  filter((STE_NAME21 != "Other Territories" & STE_NAME21 != "Outside Australia"))

victoria <- australia %>% filter(STE_NAME21 == "Victoria")

melbourne <- st_as_sf(data.frame(y = -37.813628, x = 144.963058), coords = c("x", "y"), crs = "WGS84")
melbourne$name = "Melbourne"



# Inset --------------------------------------------------------------------

inset <- ggplot()+
  geom_sf(data = australia, aes(fill = is_vic),lwd = 1)+
  scale_fill_gradient(low = "grey90", high = "grey40")+
  theme_void() +
  theme(legend.position="none")+
  coord_sf(crs = 4326)


vicmap <- ggplot()+
  geom_sf(data = victoria, fill = "grey90", lwd = 1)+
  geom_sf(data = lidar2015, fill = "grey40") +
  geom_sf(data = melbourne, shape = 15)+
  geom_sf_text(data = melbourne, aes(label = name), size = 5, nudge_y = .3, nudge_x = -.6) +
  theme_void() +
  theme(legend.position="none")

inset2 <- inset <- ggplot()+
  geom_sf(data = australia, aes(fill = is_vic),lwd = 1)+
  geom_sf(data =  st_as_sfc(st_bbox(lidar2015)), color = "red", fill = NA, lwd = 1) +
  scale_fill_gradient(low = "grey90", high = "grey40")+
  theme_void() +
  theme(legend.position="none")+
  coord_sf(crs = 4326)


# Fire severity mapping ---------------------------------------------------


firemap <- ggplot()+
  geom_sf(data = firesev, aes(fill = fire_class), color = NA) +
  scale_fill_brewer(palette = "YlOrRd", direction = -1)+
  guides(fill = guide_legend(nrow = 4))+
  theme(legend.title=element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        text=element_text(size=18),
        legend.spacing.x = unit(0, 'cm'),
        legend.margin = margin(l=-2, unit="cm"),
        plot.title = element_text(size=22),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90))+
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_orienteering,
                         pad_x = unit(.5,"cm"),
                         pad_y = unit(1,"cm"))+
  annotation_scale(location = 'bl', style = 'ticks', text_col = 'black',
                   pad_x = unit(.5,"cm"),
                   pad_y = unit(.5,"cm"),
                   text_cex = 1) +
  labs(title="a) Fire severity")

firemap



# LidarCoverage -----------------------------------------------------------

harvmap <- ggplot()+
  geom_sf(data = silv_equiv, aes(fill = System), color = NA) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme(legend.title=element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "transparent"),
        legend.background = element_blank(),
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(l=-2, unit="cm"),
        text=element_text(size=18),
        plot.title = element_text(size=22),
        legend.text = element_text(size = 16),
  axis.text.x = element_text(angle = 90))+
  annotation_north_arrow(location = "bl", which_north = "true",
                         style = north_arrow_orienteering,
                         pad_x = unit(.5,"cm"),
                         pad_y = unit(1,"cm"))+
  annotation_scale(location = 'bl', style = 'ticks', text_col = 'black',
                   pad_x = unit(.5,"cm"),
                   pad_y = unit(.5,"cm"),
                   text_cex = 1) +
  labs(title="b) Silvicultural equivalents for fire severity")

harvmap

# Compose plots -----------------------------------------------------------


layout <- "
A####
BBBBB
BBBBB
CCCCC
CCCCC
"

layout2 <-"
AAAA
AAAA
BBBB
BBBB"

map1 <- inset / vicmap

map2 <- firemap + harvmap + plot_layout(design = layout2)

map3 <- inset2 + firemap + harvmap + plot_layout(design = layout)

map2
map3

library(magick)
new <- image_trim(image_read("rplot06.png"))
image_write(new, path = "inset.png", format = "png")



inset2
