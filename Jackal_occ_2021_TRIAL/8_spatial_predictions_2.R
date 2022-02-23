## 
#===========    SPATIAL EXPLICIT OCCUPANCY ANALYSIS OF PAPER:

# The influence of the landscape context around reserves on black backed jackal occupancy across South Africa
# Lain E. Pardo, Lourens Swanepoel, Goncalo Curveira-Santos, and Jan A. Venter


# this script is to build the spatial predictions based on the best model

#1. estimate occ --------------------------------


library(readr)
library(unmarked)
library(ggplot2)
library(tidyverse)  
library(unmarked)
library(readr)

#bring original values, not scaled, better do it here, means and sd will be needed
covars  <- read.csv("data_in/siteCovs_NO_scaled_fin.csv")
jackal_DH <- read_csv("data_in/+jackal_DH_without_sites_no_temp_rep_no_outli.csv")

covars2 <- covars[,-1]

str(covars2)


# Create unmarkedFrame
covars2$Tree.Cover.Cam <- scale(covars2$Tree.Cover.Cam)  #careful, just tree cover, because that´s the only model I need to use based on previous model sel
attr(covars2$Tree.Cover.Cam, "scaled:center")
names(covars2)

y <- jackal_DH
umf = unmarkedFrameOccu(y = y, 
                        siteCovs = covars2)
head(umf)



# predictions -------------------------------------------------------------

psi_tree <-occu(~ Flash 
                ~ Tree.Cover.Cam + (1 | Season), umf) 

occ_pred <- predict(psi_tree, 
                    newdata = covars2,  
                    type = "state")#,appendData=T)

str(occ_pred) # prediction +SE and CI


------------------------------
  # 2. using a raster as input -------------------------------------------------
------------------------------
#based on Chandler_2029_Modeling and mapping species distributions (Chandler_2019_Model_mapping_sp_distributions_vignettes_steps)
# output map is trellis class not raster

library(raster)
library(rgdal)
library(RColorBrewer)
library(sf)

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

#set projection
map_proj <- st_crs(4326)

#data
tree.r <- raster("data_in/MOD44B.006_Percent_Tree_Cover_doy2018065_aid0001.tif")

GDALinfo("data_in/MOD44B.006_Percent_Tree_Cover_doy2018065_aid0001.tif")

nlayers(tree.r)
crs(tree.r)
#plot(tree.r)


tree.r.prj3 <- projectRaster(tree.r,
                             crs = crs("+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs")) #WGS 84 / UTM zone 34S

scale_fill_viridis_c(na.value = 'deeppink') #checkin NA

#double check here
res(tree.r) 
res(tree.r.prj3) 
crs(tree.r.prj3)  

# aggregate that to have a coarser scale, similar to the grid size 2.3 km 
tree.r.low <- aggregate(tree.r.prj3, fact = 10, fun = mean) 
res(tree.r.low) #now 2320 2350 so similar to our grid
#plot(tree.r.low)

# plot the map
tree.r.low #
e <- extent(-60000,2173000,6000000,7554550) #set another extend to get rid of white spaces, and have them writen in the borders of layout

tree.fin <- crop(tree.r.low,e)
plot(tree.fin, axes = TRUE) #
res(tree.fin) # plot it
crs(tree.fin) # check proj



# Standardize map values (of raster) --------------------------------------


# Since we standardized the covariates during the model fitting process, we need to transform the
# country-wide data using the same values. Note, we don't want to use the mean and SD of the rasters
# themselves, we want to use the mean and SD of the original covariates used to fit the models, which
# are stored as attributes of the sc object.
 


attr(covars2$Tree.Cover.Cam, "scaled:center")  # original mean: [1] 5.176564  

attr(covars2$Tree.Cover.Cam, "scaled:scale")  #original sd

mean(covars$Tree.Cover.Cam) #
sd(covars$Tree.Cover.Cam)  #[1] 7.267338,, ok

tree.s <- (tree.fin-5.176564)/ 7.267338

#tree.s.plot <- plot(tree.s, col= rev(terrain.colors(50)))  


# predict
(beta <- coef(psi_tree, type="state"))

logit.psi <- beta[1] + beta[2]*tree.s 
psi <- exp(logit.psi) / (1 + exp(logit.psi))
psi 

#plot(psi, col=terrain.colors(100))

pred_fin <-print(spplot(psi, col.regions= rev(terrain.colors(50)))) 
class(pred_fin)  #trellis

------------------------------
  # 3. using a raster as input not trellis map but raster -------------------------------------------------
------------------------------
  

class(tree.s)
ef <- stack(tree.s)            # need to transform the scaled map into a rasterstack
names(ef) <- "Tree.Cover.Cam"  # the raster vars must be named,
crs(tree.s) #ok
crs(ef)#ok

psi.tree.only <-occu(~ Flash 
                     ~ Tree.Cover.Cam , umf) #Cannot work with the random effect in raster

#this takes time!
#E.psi <- predict(psi.tree.only, type="state", newdata= ef) #rasterstack object
# with hash to save time when rendering this script

#need to project the predictio too then?
crs(E.psi) <- tree.r.prj3 
crs(E.psi)

raster_plot <- plot(E.psi, axes= FALSE, col=rev(terrain.colors(50))) # working fine! 4 graphs (SE...)

crs(E.psi) 

#save it as raster
#writeRaster(E.psi, "data_out/pred_raster_proj_2.tif") #


# now put the shapes ------------------------------------------------------
#load spatial pred.
E.psi <- raster("data_out/pred_raster_proj_2.tif")
raster_plot <- plot(E.psi, axes= FALSE, col=rev(terrain.colors(50)))

#get only one layer of pred and save

nlayers(E.psi)
nbands(E.psi)
crs(E.psi) 
crs(tree.s) 


# work with the raster ----------------------------------------------------

#we need to convert to a data frame first.
E.psi_df  <- as.data.frame(E.psi, xy = TRUE, na.rm=TRUE) #we had 475,475 rows, without NA now we have 229,739

str(E.psi_df)

pplot1 <- ggplot() +
  geom_histogram(data = E.psi_df, aes(pred_raster_proj_2))
pplot1 

#raster plot of prediction

plot2a <- ggplot() +
  geom_raster(data = E.psi_df,
              aes(x = x, y = y, fill = pred_raster_proj_2)) +
  scale_fill_gradientn(name = "Occ Pred", colors = rev(terrain.colors(100))) +
  labs(x = "", y = "") +
  coord_equal() 
plot2a 
#ggsave("data_out/final_results/occ_pred_map_07-08-2021.jpg")


# now include reserves on map  --------------------------------
# library(rgdal)
reserves.shp = readOGR("data_in/15_Reserves_jackal_across_SA_fin.shp")#SpatialPolygonsDataFrame


# transform shp
reserves.shp.trans <- spTransform(reserves.shp,
                                  CRSobj = crs(tree.s))
summary(reserves.shp.trans)
class(reserves.shp.trans)


library(sf)
res_sf <- reserves.shp.trans %>%    
  st_as_sf(coords = c("Long_X", "Lat_Y"), crs = crs)# %>%  #no

#plot res  
plot_res <- ggplot() +
  geom_sf(data = res_sf) #

#plot both (df of psi + sf of res)
  ggplot() +
  geom_raster(data = E.psi_df, 
              aes(x = x, y = y, fill = pred_raster_proj_2)) + 
  scale_fill_gradientn(name = "Occ pred", colors = rev(terrain.colors(50))) +
  geom_sf(data = res_sf, color = "blue", fill = NA) +
  labs(x = "", y = "") +
  coord_sf()  #reserves overlying the prediction map
  
  #ggsave("data_out/final_results/occ_pred_map_21-08-2021_geographic.jpg")

# END ---------------------------------------------------------------------

 

# rmarkdown::render("8_spatial_predictions_2.R") 








