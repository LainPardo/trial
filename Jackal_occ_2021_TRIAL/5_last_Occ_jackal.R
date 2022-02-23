## 
#===========    OCCUPANCY ANALYSIS OF PAPER:

# The influence of the landscape context around reserves on black backed jackal occupancy across South Africa
# Lain E. Pardo, Lourens Swanepoel, Goncalo Curveira-Santos, and Jan A. Venter


# NOTES
# Data preparation is not shown in this repository



library(readr)
library(unmarked)
library(tidyverse)
library(AICcmodavg)


---------------------------
####  2. occupancy analysis ------------------------------------------------
----------------------------

# data


site_covs <- read.csv("data_in/siteCovs_No_scaled_fin.csv")
jackal_DH <- read_csv("data_in/+jackal_DH_without_sites_no_temp_rep_no_outli.csv")

### scale covariates ------------------------------------------------------

names(site_covs)
site_covs$Reserve.Area <-scale(site_covs$Reserve.Area) 
site_covs$Tree.Cover.Cam <-scale(site_covs$Tree.Cover.Cam)
site_covs$Human.Density.Buffer <-scale(site_covs$Human.Density.Buffer)
site_covs$Agriculture.Buffer <-scale(site_covs$Agriculture.Buffer)
site_covs$Build.up.Buffer <- scale(site_covs$Build.up.Buffer)
site_covs$Area.PA.Buffer <-scale(site_covs$Area.PA.Buffer)
site_covs$Dist.Border.Res <-scale(site_covs$Dist.Border.Res)
site_covs$Cam.Days <-scale(site_covs$Cam.Days)

#before convert establish year to age
site_covs$Stablishment.Year <-2021- (site_covs$Stablishment.Year)
site_covs$Stablishment.Year <-scale(site_covs$Stablishment.Year)
site_covs$Season <-as.factor(site_covs$Season)
site_covs$Code.Loc <-as.factor(site_covs$Code.Loc)

str(site_covs) 
siteCovs <- site_covs[,-1]

str(siteCovs)

# naive quick
naive <- apply(X = jackal_DH,
               MARGIN = 1,
               FUN = "max", na.rm = TRUE)

mean(naive) #[1] 0.32 

# Create unmarkedFrame

y <- jackal_DH
data.umf = unmarkedFrameOccu(y = y, 
                             siteCovs = siteCovs)
                              

summary(data.umf) #

unique(siteCovs$Code.Loc)
unique(siteCovs$RSA.Vegetation)
unique(siteCovs$RSA.Biome)

# Data summaries

yy <- rowSums(y, na.rm=TRUE)  # Number of times spp were detected per site (row):sum of 1?s basically
yy
nn <- rowSums(!is.na(y))      # basically the total occassions per site,
nn  

cbind(yy, nn) #just to check y and n, per site
plot(data.umf)

-----------------------------------------------
# 3. model fitting --------------------------------------------------------
-----------------------------------------------
#Modeling selection approach = selecting for psi with general model for p (Mackenzie)...
  

# 1. constant MODEL

M0 <- occu(~1 ~1, data.umf)
backTransform(M0, type="state")  # 0.381 SE 0.0335
backTransform(M0, type="det") # 0.27 SE 0.019

# null
psi_nul <-occu(~Flash + Cam.Days + Tree.Cover.Cam
                 ~ 1, data.umf)


## 2. reserve level (macro scale): reserve features - same covariate value to all cameras within a reserve ---------------
# season is acting as a random factor

# global out: disturbance
class(siteCovs$Season)
class(siteCovs$Code.Loc)
names(siteCovs)

global_out <-occu(~Flash + Cam.Days + Tree.Cover.Cam
                  ~ Human.Density.Buffer  + Agriculture.Buffer +  Build.up.Buffer +
                    (1 | Season), data.umf) 

global_out


#global inside: if variables inside are more important than outside


global_in_res <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                  ~ Reserve.Area  + Stablishment.Year + Tree.Cover.Res + (1 | Season), data.umf) 
global_in_res

# inside: reserve level individual

psi_area_res <-occu(~Flash + Cam.Days + Tree.Cover.Cam
                    ~ Reserve.Area + (1 | Season), data.umf) 
psi_area_res

#
psi_year <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                ~ Stablishment.Year +  (1 | Season), data.umf) 
psi_year

#
psi_tree.res <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                    ~ Tree.Cover.Res +  (1 | Season), data.umf) 
psi_tree.res


# Buffer outside: reserve level individual
psi_human <-occu(~ Flash + Cam.Days +Tree.Cover.Cam
                 ~ Human.Density.Buffer + (1 | Season), data.umf)  
psi_human

#

psi_agri <-occu(~ Flash + Cam.Days +Tree.Cover.Cam
                ~ Agriculture.Buffer + (1 | Season), data.umf)  
psi_agri

#

psi_build <-occu(~ Flash + Cam.Days +Tree.Cover.Cam
                 ~ Build.up.Buffer + (1 | Season), data.umf)  
psi_build

#

psi_area_PA <-occu(~ Flash + Cam.Days +Tree.Cover.Cam
                   ~ Area.PA.Buffer + (1 | Season), data.umf)  
psi_area_PA




# 3. site level:  features at camera scale - different value per site

global_in_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                     ~ Tree.Cover.Cam + Dist.Border.Res + (1 | Season), data.umf) 

global_in_site


# individual factors
psi_tree_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                ~ Tree.Cover.Cam + (1 | Season), data.umf)
psi_tree_site

#
psi_border <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                  ~ Dist.Border.Res + (1 | Season), data.umf)
psi_border



# 4. carnivores models at site level ---------------------------------------

#
psi_lion_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                 ~ Lion.Site + (1 | Season), data.umf)
psi_lion_site

# 
psi_leopard_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                 ~ Leopard.Site + (1 | Season), data.umf)
psi_leopard_site

#
psi_hyena_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                 ~ Hyena.Site + (1 | Season), data.umf)
psi_hyena_site


# psi_all_pred_site
psi_all_pred_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                       ~Lion.Site + Leopard.Site + Hyena.Site + (1 | Season), data.umf)
psi_all_pred_site 

# pred richness site
psi_richness_site <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                         ~ Pred.Rich.Site + (1 | Season), data.umf)
summary(psi_richness_site)  


# 5. carnivore models at reserve level ------------------------------------

psi_lion_reserve <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                     ~ Lion.Reserve + (1 | Season), data.umf)
psi_lion_reserve

# 
psi_leopard_reserve <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                        ~ Leopard.Reserve + (1 | Season), data.umf)
psi_leopard_reserve

#
psi_hyena_reserve <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                      ~ Hyena.Reserve + (1 | Season), data.umf)
psi_hyena_reserve

# psi_all_pred_reserve
psi_all_pred_res <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                         ~Lion.Reserve + Leopard.Reserve + Hyena.Reserve + (1 | Season), data.umf)
psi_all_pred_res 

# richness reserve
psi_richness_reserve <-occu(~ Flash + Cam.Days + Tree.Cover.Cam
                            ~ Pred.Rich.Reserve + (1 | Season), data.umf)
summary(psi_richness_reserve) 


# model selection all models -----------------------


#reserve features: models 1-9
# site level features: models 10-12
# community site 13-16
# community reserve 17-20

Cand.mods1 <- list("p(Flash + Cam.Days + tree)psi(global.out)" = global_out,
                   "p(Flash + Cam.Days + tree)psi(global.in.res)" = global_in_res,
                   "p(Flash + Cam.Days + tree)psi(human.out)" = psi_human,
                   "p(Flash + Cam.Days + tree)psi(agri.out)" = psi_agri,
                   "p(Flash + Cam.Days + tree)psi(buildup.out)" = psi_build,    
                   "p(Flash + Cam.Days + tree)psi(area.res.out)" = psi_area_PA,
                   "p(Flash + Cam.Days + tree)psi(area.res)" = psi_area_res,
                   "p(Flash + Cam.Days + tree)psi(year.res)" = psi_year,
                   "p(Flash + Cam.Days + tree)psi(tree.cov.res)" = psi_tree.res,
                   
                   "p(Flash + Cam.Days + tree)psi(global.in.site)" = global_in_site,
                   "p(Flash + Cam.Days + tree)psi(tree.cov.site)" = psi_tree_site,
                   "p(Flash + Cam.Days + tree)psi(dist.border.site)" = psi_border,
                  
                   "p(Flash + Cam.Days + tree)psi(lion.site)" = psi_lion_site,
                   "p(Flash + Cam.Days + tree)psi(leopard.site)" = psi_leopard_site,
                   "p(Flash + Cam.Days + tree)psi(hyena.site)" = psi_hyena_site,
                   "p(Flash + Cam.Days + tree)psi(pred.rich.site)" = psi_richness_site,
                  
                   "p(Flash + Cam.Days + tree)psi(lion.res)" = psi_lion_reserve,
                   "p(Flash + Cam.Days + tree)psi(leopard.res)" = psi_leopard_reserve,
                   "p(Flash + Cam.Days + tree)psi(hyena.res)" = psi_hyena_reserve,
                  "p(Flash + Cam.Days + tree)psi(pred.rich.res)" = psi_richness_reserve,
                 
                   
                   "p(Flash + Cam.Days + tree)psi(.)" = psi_nul,
                   "p(.)psi(.)" = M0)

                
(model_sel_all <- aictab(Cand.mods1))
(model_sel_all <- aictab(Cand.mods1,rank="QAICc", chat=1.34))  
(model_sel_table <- as.data.frame(model_sel_all))
#write.csv(model_sel_all, "data_out/final_results/model_sel_psi_last_QAICc.csv")

# MODEL SELECTION FOR P BASED ON THESE RESULTS only two models below 2AIC
flash_p <-occu(~Flash 
               ~ Tree.Cover.Cam + (1 | Season) , data.umf)
camdays_p <- occu(~Cam.Days 
                  ~ Tree.Cover.Cam + (1 | Season) , data.umf)
tree_p <- occu(~Tree.Cover.Cam 
               ~ Tree.Cover.Cam + (1 | Season) , data.umf)

null_p <- occu(~1 
               ~ Tree.Cover.Cam + (1 | Season) , data.umf)

#for second mod
flash_p2 <- occu(~ Flash 
                    ~ Tree.Cover.Cam + Dist.Border.Res + (1 | Season), data.umf) 

camdays_p2 <- occu(~ Cam.Days 
                    ~ Tree.Cover.Cam + Dist.Border.Res + (1 | Season), data.umf) 
tree_p2 <- occu(~ Tree.Cover.Cam 
                    ~ Tree.Cover.Cam + Dist.Border.Res + (1 | Season), data.umf) 
null_p2 <- occu(~ 1
           ~ Tree.Cover.Cam + Dist.Border.Res + (1 | Season), data.umf) 


Cand.mods_p <- list("p(Flash)psi(tree.cover.site)" = flash_p,
                    "p(Cam.Days)psi(tree.cover.site)" =camdays_p,
                    "p(tree)psi(tree.cover.site)" = tree_p, 
                    
                    
                    "p(Flash)psi(global.in.site)" = flash_p2,
                    "p(Cam.Days)psi(global.in.site)" =camdays_p2,
                    "p(tree)psi(global.in.site)" = tree_p2, 
                    
                    "p(.)psi(tree.cover.site)" = null_p,
                    "p(.)psi(global.in.site)" = null_p2,
                    "p(.)psi(.)" = M0)
                    
(model_sel_all2 <- aictab(Cand.mods_p))  
(model_sel_all2 <- aictab(Cand.mods_p,rank="QAICc", chat=1.34 ))  
(model_sel_table <- as.data.frame(model_sel_all))
#write.csv(model_sel_all2, "data_out/model_sel_with_p.csv")



#godness of fit ---------------
coef(global_out)
global_out@TMB$par
global_out@TMB$starts_order
 
#gofdet1 <- mb.gof.test(global_out, nsim = 10000, plot.hist=FALSE) #
#gofdet1
  
# Chi-square statistic = 3905.671 
# Number of bootstrap samples = 10000
# P-value = 0.1238
#   
# Quantiles of bootstrapped statistics:
#     0%    25%    50%    75%   100% 
#   558   1177   1534   2304 499110 
#   
#  Estimate of c-hat = 1.34

 ------------------------------------
  # 4. PREDICTIONS -------------------------------------------------------------
  -----------------------------------

# BEST MODEL NOW:

flash_p <- occu(~ Flash
                ~ Tree.Cover.Cam + (1 | Season), data.umf)

range(siteCovs$Tree.Cover.Cam)
  
newData  <-  data.frame(Tree.Cover.Cam = seq(-0.7123054, 4.6376593, length = 8), 
                          Season = c("autunm","autunm-winter", "spring", "spring-summer", "summer", "summer-autunm", "winter", 
                                   "winter-spring")) 
                          
  pred <- predict(flash_p,  type  =  "state",  newdata  =  newData,  appendData  =  TRUE)
  pred


# PLOT IT
  
  f1 <- ggplot(pred,aes(Tree.Cover.Cam, Predicted, ymin=lower,ymax=upper))
  f1 + geom_line(colour = "blue") + 
    geom_ribbon(alpha=0.2,colour=NA) + 
    labs(x = "Tree cover (std)", y = "Prob. occupancy") +
    theme_bw() +
    theme(axis.text=element_text(size=15),axis.line=element_line(colour="black"))+
    ylim(0,1)
  
  
  #ggsave("figures/pred_plot_with_season.jpg")
  
# alternative figure without random factor season 
  
  flash_NS <-occu(~Flash 
                 ~ Tree.Cover.Cam, data.umf)
  newData1  <-  data.frame(Tree.Cover.Cam = seq(-0.7123054, 4.6376593, by = 0.1)) 
                         
  pred2 <- predict(flash_NS,  type  =  "state",  newdata  =  newData1,  appendData  =  TRUE)
  pred2
  
   
  # PLOT IT
  
  f2 <- ggplot(pred2,aes(Tree.Cover.Cam, Predicted, ymin=lower,ymax=upper))
  f2 + geom_line(colour = "blue") + 
    geom_ribbon(alpha=0.2,colour=NA) + 
    labs(x = "Tree cover (std)", y = "Prob. occupancy") +
    theme_bw() +
    theme(axis.text=element_text(size=15),axis.line=element_line(colour="black"))+
    ylim(0,1)
  
  
  #ggsave("figures/pred_plot_NO_season_std.jpg")



-------------------------------
# average prediction per reserve ---------------------------------------------------
-------------------------------
 
#predict again but put all columns with appendData
#using season as random, model flash_p

class(siteCovs) 
siteCovs <- as.data.frame(siteCovs)
    
predmap <- predict(flash_p, newdata = siteCovs, type = "state", appendData=T)   
      
pred_all <- as.data.frame(predmap)
   
   
library("plotrix") #to calc se
   
#add col with code
avg_pred <-   pred_all %>% 
     group_by(Code.Loc) %>% 
     summarise(psi_mean = mean(Predicted),
               psi_sd = sd(Predicted),
               psi_se = std.error(Predicted), 
               min =  psi_mean - psi_sd,
               max =  psi_mean + psi_sd)
avg_pred
# write.csv(avg_pred,"data_out/avg_occ_persite.csv")
   

# ordered catterpillar plot
    
str(avg_pred)
occgraph <- avg_pred
   
occgraph$Code.Loc <- factor(occgraph$Code.Loc, levels = occgraph$Code.Loc, ordered = TRUE)
   
str(occgraph)
   
rr4 <- transform(occgraph,Reserve = reorder(Code.Loc,psi_mean))
theme_set(theme_bw())
   
ggplot(rr4,aes(Reserve,psi_mean,ymin=min,ymax=max))+
     geom_pointrange(data = rr4, aes(ymin =psi_mean - psi_se, ymax = psi_mean + psi_se), 
                     alpha = 2, width = 4, color = "blue") +
     #c +
     # xlim=c(-1.5,1.5) +
     xlab("Reserves") +
     ylab("Average estimated Occupancy probability") + 
     theme(axis.text=element_text(size=11, color =  "black")) +
     coord_flip() +
     theme(legend.position='none') + 
     geom_hline(yintercept= 0.5, linetype="dashed", 
                color = "grey", size=0.7)
   
#ggsave("figures/catterpillar_avg_psi_res1.jpg")
   
    

# END ---------------------------------------------------------------------


# knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) 

# rmarkdown::render("5_last_Occ_jackal.R") not working...Error in -`*tmp*` : invalid argument to unary operator
# (Ctrl+Shift+K)
   