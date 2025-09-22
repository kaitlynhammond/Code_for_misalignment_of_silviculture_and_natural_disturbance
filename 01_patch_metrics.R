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

sev_rast  <- df %>%
  filter(!is.na(str2007), !is.na(CLASS)) %>%
  mutate(code = case_when(transition == "1 to 1" ~ 0,
                          transition == "1 to 2" ~ 0,
                          transition == "1 to 3" ~ 1,
                          transition == "2 to 1" ~ 0,
                          transition == "2 to 2" ~ 0,
                          transition == "2 to 3" ~ 1))
sev_rast <- rast(sev_rast[c(1,2,7)], type = "xyz")
plot(sev_rast)
# Harvest -----------------------------------------------------------------

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


harv_rast <- rasterize(logarea, rast2016, background = "0")

harv_rast[is.na(harv_rast)] <- 0
plot(harv_rast)



# Patch-level metrics -----------------------------------------------------
metrics <- list_lsm(level = "patch")

sev_rast[is.na(sev_rast)] <- 0
harv_rast[is.na(harv_rast)] <- 0

lsm_fire <- calculate_lsm(sev_rast, #input raster
                     what = c("lsm_p_area", #patch area
                              "lsm_p_cai",#core area index
                              "lsm_p_enn",#Euclidean nearest neighbour distance
                              "lsm_p_perim", #patch perimeter
                              "lsm_p_frac", #fractal dimesion index
                              "lsm_p_shape"), #shape index
                     directions = 8)

lsm_fire$history <- "fire"
lsm_fire <- lsm_fire %>%
  filter(class == 1)



lsm_harv <- calculate_lsm(harv_rast, #input raster
                          what = c("lsm_p_area", #patch area
                                   "lsm_p_cai",#core area index
                                   "lsm_p_enn",#Euclidean nearest neighbour distance
                                   "lsm_p_perim", #patch perimeter
                                   "lsm_p_frac", #fractal dimesion index
                                   "lsm_p_shape"), #shape index
                          directions = 8)

lsm_harv$history <- "harvest"
lsm_harv <- lsm_harv %>%
  filter(class ==1)




# group and make frames ---------------------------------------------------
harv_ids <- unique(lsm_harv$id)
harv_ids <- sample(harv_ids, 1000, replace = FALSE, prob = NULL)


harv_subset <- lsm_harv %>%
  filter(id %in% harv_ids)

fire_ids <- unique(lsm_fire$id)
fire_ids <- sample(fire_ids, 1000, replace = FALSE, prob = NULL)


fire_subset <- lsm_fire %>%
  filter(id %in% fire_ids)


all_lsm <- rbind(fire_subset, harv_subset)

patch_area <- all_lsm %>%
  filter(metric == "area") %>%
  select(id, value, history)
names(patch_area) <- c("id", "area", "history")

core_area_i <- all_lsm %>%
  filter(metric == "cai") %>%
  select(id, value, history)
names(core_area_i) <- c("id", "cai", "history")

enn <- all_lsm %>%
  filter(metric == "enn") %>%
  select(id, value, history)
names(enn) <- c("id", "enn", "history")

perim <- all_lsm %>%
  filter(metric == "perim") %>%
  select(id, value, history)
names(perim) <- c("id", "perim", "history")

frac <- all_lsm %>%
  filter(metric == "frac") %>%
  select(id, value, history)
names(frac) <- c("id", "frac", "history")

shape <- all_lsm %>%
  filter(metric == "shape") %>%
  select(id, value, history)
names(shape) <- c("id", "shape", "history")


# Correlation -------------------------------------------------------------
check_cor <- left_join(patch_area, core_area_i)
check_cor <- left_join(check_cor, enn)
check_cor <- left_join(check_cor, perim)
check_cor <- left_join(check_cor, frac)
check_cor <- left_join(check_cor, shape)

corMat<-cor(check_cor[,c(4,5,6)], use="complete.obs", method = "pearson")
corrplot(corMat, 
         method="shade",
         type="lower",
         diag = FALSE,
         addCoef.col = "black")


# Plots -------------------------------------------------------------------


a <- ggplot(patch_area)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=area, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_x_log10(breaks= c(0.01,.1,1,10,100,1000), limits = c(0.01, 10000), labels = label_number(), expand = c(0,0)) +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Patch area", y = "Proportion")

b <- ggplot(core_area_i)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=cai, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Core area index", y = "Proportion")

c <- ggplot(enn)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=enn, color = history, fill = history), alpha = 0.5, position = "identity")+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_x_log10(breaks= c(0,1,10,100,1000,10000,100000), labels = label_number()) +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = str_wrap("Euclidean nearest neighbour distance (m)", width = 20), y = "Proportion")

d <- ggplot(perim)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=perim, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  scale_x_log10(breaks= c(0,1,10,100,1000,10000,100000), labels = label_number()) +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Patch perimeter (m)", y = "Proportion")

e <- ggplot(frac)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=frac, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Fractal dimension index", y = "Proportion")

f <- ggplot(shape)+
  geom_histogram(aes(y = after_stat(count / sum(count)),x=shape, color = history, fill = history), alpha = 0.5, position = "identity", show.legend = FALSE)+
  scale_fill_viridis(discrete = TRUE, option = "C") +
  scale_color_viridis(discrete = TRUE, option = "C") +
  annotation_logticks(sides = "b",
                      outside = TRUE) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank()) +
  labs(x = "Shape index", y = "Proportion")


layout <- "
ABCD
"



a + b + d + c + plot_layout(design = layout, axis_titles = "collect")

