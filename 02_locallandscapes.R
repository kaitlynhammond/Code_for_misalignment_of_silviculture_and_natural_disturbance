library("landscapemetrics")
library("terra")
library("sf")
library("viridis")
library("tidyverse")
library("patchwork")
library("scales")
library("corrplot")


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

sev_rast <- df %>%
  filter(!is.na(str2007), !is.na(CLASS)) %>%
  mutate(code = case_when(transition == "1 to 1" ~ 0,
                          transition == "1 to 2" ~ 0,
                          transition == "1 to 3" ~ 1,
                          transition == "2 to 1" ~ 0,
                          transition == "2 to 2" ~ 0,
                          transition == "2 to 3" ~ 1))
sev_rast <- rast(sev_rast[c(1,2,7)], type = "xyz")

# Spatial pattern ---------------------------------------------------------


buffer1000 <- st_read("spatial/buffer1000.shp")
buffer1000 <- buffer1000 %>%
  select(geometry)


all_rasters <- list()

for(i in 1:nrow(buffer1000)){
  
  local <- buffer1000[i,]
  landscape <- crop(rast2016, local)
  landscape[is.na(landscape)] <- 0
  
  all_rasters[[i]] <- landscape
}

metrics <- list_lsm(level = "class")


# Landscape metrics local landscapes --------------------------------------
lsm_list <- list()

for(i in 1:nrow(buffer1000)){
  
  current_rast <- all_rasters[[i]]
  current_rast[current_rast >= 2] <- 2
  
  lsm <- calculate_lsm(current_rast, #input raster
                       what = c("lsm_c_ai", #aggregation index
                                "lsm_c_clumpy",#Clumpiness index
                                "lsm_c_lsi",#landscape shape index
                                "lsm_c_np", #number of patches
                                "lsm_c_pd", #patch density
                                "lsm_c_ed", #Edge density
                                "lsm_c_lpi", #largest patch index
                                "lsm_c_enn_mn",#mean euclidean nearest neighbour)
                                "lsm_c_area_mn",#mean patch area,
                                "lsm_c_te", #total edge
                                "lsm_c_frac_mn", #mean fractal dimension index
                                "lsm_c_shape_mn"), #mean shape index
                       directions = 8) #This is the neighborhood rule we are using
  
  lsm$row <- i
  lsm_list[[i]] <- lsm
  
  
}



lsm_table <- do.call(rbind, lsm_list)


class_area <- lsm_table %>%
  filter(class == 1, metric == "ai") %>%
  select(row, value)
names(class_area) <- c("id", "ai")

lsi <- lsm_table %>%
  filter(class == 1, metric == "lsi") %>%
  select(row, value)
names(lsi) <- c("id", "lsi")

npatch <- lsm_table %>%
  filter(class == 1, metric == "np") %>%
  select(row, value)
names(npatch) <- c("id", "np")

pd <- lsm_table %>%
  filter(class == 1, metric == "pd") %>%
  select(row, value)
names(pd) <- c("id", "pd")

ed <- lsm_table %>%
  filter(class == 1, metric == "ed") %>%
  select(row, value)
names(ed) <- c("id", "ed")

lpi <- lsm_table %>%
  filter(class == 1, metric == "lpi") %>%
  select(row, value)
names(lpi) <- c("id", "lpi")

area_mn <- lsm_table %>%
  filter(class == 1, metric == "area_mn") %>%
  select(row, value)
names(area_mn) <- c("id", "area_mn")

te <- lsm_table %>%
  filter(class == 1, metric == "te") %>%
  select(row, value)
names(te) <- c("id", "te")

frac_mn <- lsm_table %>%
  filter(class == 1, metric == "frac_mn") %>%
  select(row, value)
names(frac_mn) <- c("id", "frac_mn")

enn_mn <- lsm_table %>%
  filter(class == 1, metric == "enn_mn") %>%
  select(row, value)
names(enn_mn) <- c("id", "enn_mn")

shape <- lsm_table %>%
  filter(class == 1, metric == "shape_mn") %>%
  select(row, value)
names(shape) <- c("id", "shape")


fire_metrics <- full_join(class_area, lsi)
fire_metrics <- full_join(fire_metrics, npatch)
fire_metrics <- full_join(fire_metrics, pd)
fire_metrics <- full_join(fire_metrics, ed)
fire_metrics <- full_join(fire_metrics, lpi)
fire_metrics <- full_join(fire_metrics, area_mn)
fire_metrics <- full_join(fire_metrics, te)
fire_metrics <- full_join(fire_metrics, frac_mn)
fire_metrics <- full_join(fire_metrics, enn_mn)
fire_metrics <- full_join(fire_metrics, shape)

fire_metrics$history <- "fire"




# Harvest lsm -------------------------------------------------------------

loghist <- st_read("spatial/LASTLOG25.shp")
loghist <- loghist %>% 
  mutate(year = as.numeric(substr(SEASON, 1 , 4))+1, 
         area = as.numeric(st_area(loghist)/10000),
         history = "harvest")
logarea <- loghist %>%
  filter(year >= 1974, X_FORETYPE == "Mountain Ash",
         X_SILVSYS == "Clearfelling" |
           X_SILVSYS == "Seed Tree (includes retained overwood)")

logarea <- logarea[,26:27]
logarea <- vect(logarea)
logarea <- project(logarea, rast2016)
logarea <- crop(logarea, rast2016)


gridrast <- rast(extent = ext(rast2016), crs = crs(rast2016), resolution = 1000, vals = 1)
gridpoly <- as.polygons(gridrast, aggregate = FALSE)
gridpoly <- mask(gridpoly, logarea)
gridpoly <- st_as_sf(gridpoly)
gridpoly$area <- st_area(gridpoly)
plot(gridpoly)


logarea <- rasterize(logarea, rast2016, background = "0")


harv_patches <- lsm_p_area(logarea, directions = 8) %>%
  filter(class == 1)


all_rasters <- list()

for(i in 1:nrow(gridpoly)){
  
  local <- gridpoly[i,]
  landscape <- crop(logarea, local)
  landscape[is.na(landscape)] <- 0
  
  all_rasters[[i]] <- landscape
}




# Landscape metrics local landscapes --------------------------------------
lsm_list <- list()

for(i in 1:nrow(gridpoly)){
  
  current_rast <- all_rasters[[i]]
  
  lsm <- calculate_lsm(current_rast, #input raster
                       what = c("lsm_c_ai", #aggregation index
                                "lsm_c_clumpy",#Clumpiness index
                                "lsm_c_lsi",#landscape shape index
                                "lsm_c_np", #number of patches
                                "lsm_c_pd", #patch density
                                "lsm_c_ed", #Edge density
                                "lsm_c_lpi", #largest patch index
                                "lsm_c_enn_mn",#mean euclidean nearest neighbour)
                                "lsm_c_area_mn",#mean patch area,
                                "lsm_c_te", #total edge
                                "lsm_c_frac_mn", #mean fractal dimension index
                                "lsm_c_shape_mn"), #mean shape index
                       directions = 8)  #This is the neighborhood rule we are using
  
  lsm$row <- i
  lsm_list[[i]] <- lsm
  
  
}


lsm_table <- do.call(rbind, lsm_list)


class_area <- lsm_table %>%
  filter(class == 1, metric == "ai") %>%
  select(row, value)
names(class_area) <- c("id", "ai")

lsi <- lsm_table %>%
  filter(class == 1, metric == "lsi") %>%
  select(row, value)
names(lsi) <- c("id", "lsi")

npatch <- lsm_table %>%
  filter(class == 1, metric == "np") %>%
  select(row, value)
names(npatch) <- c("id", "np")

pd <- lsm_table %>%
  filter(class == 1, metric == "pd") %>%
  select(row, value)
names(pd) <- c("id", "pd")

ed <- lsm_table %>%
  filter(class == 1, metric == "ed") %>%
  select(row, value)
names(ed) <- c("id", "ed")

lpi <- lsm_table %>%
  filter(class == 1, metric == "lpi") %>%
  select(row, value)
names(lpi) <- c("id", "lpi")

area_mn <- lsm_table %>%
  filter(class == 1, metric == "area_mn") %>%
  select(row, value)
names(area_mn) <- c("id", "area_mn")

te <- lsm_table %>%
  filter(class == 1, class == 1, metric == "te") %>%
  select(row, value)
names(te) <- c("id", "te")

frac_mn <- lsm_table %>%
  filter(class == 1, metric == "frac_mn") %>%
  select(row, value)
names(frac_mn) <- c("id", "frac_mn")

enn_mn <- lsm_table %>%
  filter(class == 1, metric == "enn_mn") %>%
  select(row, value)
names(enn_mn) <- c("id", "enn_mn")

shape <- lsm_table %>%
  filter(class == 1, metric == "shape_mn") %>%
  select(row, value)
names(shape) <- c("id", "shape")


harv_metrics <- full_join(class_area, lsi)
harv_metrics <- full_join(harv_metrics, npatch)
harv_metrics <- full_join(harv_metrics, pd)
harv_metrics <- full_join(harv_metrics, ed)
harv_metrics <- full_join(harv_metrics, lpi)
harv_metrics <- full_join(harv_metrics, area_mn)
harv_metrics <- full_join(harv_metrics, te)
harv_metrics <- full_join(harv_metrics, frac_mn)
harv_metrics <- full_join(harv_metrics, enn_mn)
harv_metrics <- full_join(harv_metrics, shape)

harv_metrics$history <- "harvest"





# Combine and plot -------------------------------------------------------

harv_subset <- harv_metrics %>%
  filter(!is.na(enn_mn))
fire_subset <- fire_metrics %>%
  filter(!is.na(enn_mn))


all_metrics <- rbind(fire_metrics, harv_metrics) %>%
  group_by(history)%>%
  slice_sample(n = 800)

all_metrics2 <- rbind(fire_subset, harv_subset) %>%
  group_by(history) %>%
  slice_sample(n = 700)


# Correlation -------------------------------------------------------------

corMat<-cor(all_metrics[,c(2,4,6,7,8,11,12)], use="complete.obs", method = "pearson")
corrplot(corMat, 
         method="shade",
         type="lower",
         diag = FALSE,
         addCoef.col = "black")


# Plots -------------------------------------------------------------------

a <- ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=ai, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Aggregation index", y = "Proportion")



b <- ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=lsi, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Landscape shape index", y = "Proportion")



c <- ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=np, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Number of patches", y = "Proportion")



d <- ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=pd, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C")+
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Patch density", y = "Proportion")



e <- ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=ed, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Edge density (m/ha)", y = "Proportion")



f <-  ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=lpi, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Largest patch index", y = "Proportion")

g <-  ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=area_mn, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Mean patch area (ha)", y = "Proportion")

h <-  ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=te, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  scale_x_log10(breaks= c(0,1,10,100,1000,10000), labels = label_number(), expand = c(0,0)) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Total edge", y = "Proportion")


i <-  ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=frac_mn, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_colour_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = str_wrap("Mean fractal dimension index", width = 20) , y = "Proportion", label_wrap(15))

j <-  ggplot(all_metrics2)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=enn_mn, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
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
  labs(x = str_wrap("Mean euclidean nearest neighbour distance (m)", width = 30), y = "Proportion")

k <-  ggplot(all_metrics)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=shape, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
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
  labs(x = "Mean shape index", y = "Proportion")


l <- cowplot::get_legend(ggplot(all_metrics)+
                           geom_histogram(aes(x=shape, color = history, fill = history), alpha = 0.5, position = "identity")+
                           scale_fill_viridis(discrete = TRUE, option = "C")  +
                           scale_colour_viridis(discrete = TRUE, option = "C") +
                           theme(legend.title=element_blank()))

layout <- "
ABCD
EFGH
"

c + g + f + a +plot_layout(design = layout, axis_titles = "collect_y")
 
a + c + e + f + g + j + k + l+ plot_layout(design = layout, axis_titles = "collect_y")

c + g + f + e + a + k + j + l+ plot_layout(design = layout, axis_titles = "collect_y")

