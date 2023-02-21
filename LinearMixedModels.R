###Patterns and drivers of climatic niche dynamics during biological invasions of insular tetrapods GLMMs 

pkgs <- c("viridis","glmm","rcompanion","lmtest", "spdep","lme4", 
          "Hmisc","reshape2","nlme","r2glmm","MuMIn", "sjPlot", "sjmisc", 
          "glmmTMB", "ggpubr", "ggridges", "ggplot2","plotly",
          "ggcorrplot", "dplyr")
sapply(pkgs, require, character.only = TRUE)

#Load data
data<-read.csv("~/Data_Clusters.csv")
hist(data$mean_NMI)

###Remove outliers 
Q <- quantile(data$mean_NMI, probs=c(.05, .95), na.rm = FALSE)
iqrNMI<-IQR(data$mean_NMI)
dataNoOut<- subset(data, data$mean_NMI > (Q[1] - 1.5*iqrNMI) & data$mean_NMI < (Q[2]+1.5*iqrNMI))
clean.dat<-dataNoOut ##112 obs

###Check for multicolinearity in the drivers
###Drivers
drivers<-na.omit(data[,c(6,7,8,10,11,12,13,14)])
names(drivers)<-c("AltRange", "Roughness", "Remoteness", "Dist_Nat-Alien", "PD_Rec_comm","BodyMass", "Nat_RangeSize", "ED_aliensp")

###CORRELATION ANALYSIS for the drivers tested.
# Melt the correlation matrix
library(reshape2)
cormat <- round(cor(drivers),2)
head(cormat)
melted_cormat <- melt(cormat)
head(melted_cormat)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "dodgerblue3", high = "firebrick3", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
  panel.grid.major = element_blank(),panel.border = element_blank(),
  panel.background = element_blank(),axis.ticks = element_blank(),
  legend.justification = c(1, 0),legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+guides(fill = guide_colorbar(barwidth = 6.5, 
  barheight = 0.5,title.position = "top", title.hjust = 0.5))

####Linear Mixed Models
dat.scaled<-clean.dat %>%
  mutate_if(is.numeric, scale) ###SCALE ALL DATA
dat.scaled %>%
  summarise_if(is.numeric, ~ near(mean(.x, na.rm=T), y = 0)) ###CHECK ALL MEANS=0
dat.scaled  %>%
  summarise_if(is.numeric, ~ near(sd(.x, na.rm=T), y = 1)) ######CHECK ALL SD=1

library(lmerTest)
##RUN FULL OVERALL MODEL (ALL DATA+ALL VARS)
full.model <- lmer(mean_NMI~ AltRange + Rough + Remot +  RecInsul + DistNatAlien + PD + BodyMass + NatRangeSize + ED.Alien  + (1|class/species), data=na.omit(dat.scaled))
summary(full.model)

#Results table
res.full.mod.all<-round(summary(full.model)$coefficients, 2)

###Plot all full model
par(mfrow=c(2,2))
effsize_all<-plot_model(full.model, type = "std", show.intercept =FALSE, grid.breaks = 1, ci.style = "whisker",title = "All taxa", line.size = 1.1, dot.size = 3)+ theme_cleveland()
effsize_all

###AMPHIBIANS
amph<-na.omit(clean.dat[clean.dat$class=="Amphibia",])
amph.scaled<-na.omit(clean.dat[clean.dat$class=="Amphibia",]) %>%
  mutate_if(is.numeric, scale) 
amph.scaled %>%
  summarise_if(is.numeric, ~ near(mean(.x, na.rm=T), y = 0)) ###CHECK ALL MEANS=0
amph.scaled  %>%
  summarise_if(is.numeric, ~ near(sd(.x, na.rm=T), y = 1)) ######CHECK ALL SD=1

#Model
full.model_A <- lmer(mean_NMI~AltRange + Rough + Remot +  RecInsul + DistNatAlien + PD + BodyMass + NatRangeSize + ED.Alien +  (1|species), data=amph.scaled)
summary(full.model_A)
#Results table
res.full.mod.A<-round(summary(full.model_A)$coefficients, 2)

#Plot
effsize_amphs<-plot_model(full.model_A, type = "std", show.intercept =FALSE, grid.breaks = 12, ci.style = "whisker",title = "Amphibians", line.size = 1.1, dot.size = 3) + theme_cleveland()
effsize_amphs

###BIRDS
bird.scaled<-na.omit(clean.dat[clean.dat$class=="Aves",]) %>%
  mutate_if(is.numeric, scale) 
bird.scaled %>%
  summarise_if(is.numeric, ~ near(mean(.x, na.rm=T), y = 0)) ###CHECK ALL MEANS=0
bird.scaled  %>%
  summarise_if(is.numeric, ~ near(sd(.x, na.rm=T), y = 1)) ######CHECK ALL SD=1

#Model
full.model_B <- lmer(mean_NMI~AltRange + Rough + Remot +  RecInsul + DistNatAlien + PD + BodyMass + NatRangeSize + ED.Alien  +  (1|species), data=bird.scaled)
summary(full.model_B)
#Results table
res.full.mod.B<-round(summary(full.model_B)$coefficients, 2)

##Plot Birds full model
effsize_birds<-plot_model(full.model_B, type = "std", show.intercept =FALSE, grid.breaks = 1, ci.style = "whisker",title = "Birds",line.size = 1.1, dot.size = 3) + theme_cleveland()
effsize_birds

###REPTILES
rept.scaled<-na.omit(clean.dat[clean.dat$class=="Reptilia",]) %>%
  mutate_if(is.numeric, scale) 
rept.scaled %>%
  summarise_if(is.numeric, ~ near(mean(.x, na.rm=T), y = 0)) ###CHECK ALL MEANS=0
rept.scaled  %>%
  summarise_if(is.numeric, ~ near(sd(.x, na.rm=T), y = 1)) ######CHECK ALL SD=1

#Model
full.model_R <- lmer(mean_NMI~AltRange + Rough + Remot +  RecInsul + DistNatAlien + PD + BodyMass + NatRangeSize + ED.Alien +  (1|species), data=rept.scaled)
summary(full.model_R)

#Table results
res.full.mod.R<-round(summary(full.model_R)$coefficients, 2)

#Plot repts full model
effsize_repts<-plot_model(full.model_R, type = "std", show.intercept =FALSE, grid.breaks = 1, ci.style = "whisker",title = "Reptiles", line.size = 1.1, dot.size = 3) + theme_cleveland()
effsize_repts

##Plot all models
ggarrange(effsize_all,effsize_amphs,effsize_birds, effsize_repts, ncol = 2, nrow = 2)
