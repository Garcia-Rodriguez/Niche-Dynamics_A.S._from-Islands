rm(list=ls())
gc()

# Analyses:Estimation of niche dynamics (COUE protocol) 
# Paper: Patterns and drivers of climatic niche dynamics in biological invasions of insular tetrapods
# Author: Adrian Garcia 

#Load packages
PKGs <- c("viridis","raster","ecospat","rlang", "rgbif", "dplyr","scrubr", "rgdal","maps", "maptools", "exactextractr","sp","spocc","plyr", "tidyverse", "pbapply", "usdm", "geosphere", "pbapply", "sf", "ggplot2", "rgeos", "ade4")

sapply(PKGs, require, character.only = TRUE)


setwd("C:/Users/garciaa/Dropbox/PROJECTS_VIENNA/CLIMATIC DYNAMICS/AHORASI/")
setwd("C:/Users/Lenovo/Dropbox/PROJECTS_VIENNA/CLIMATIC DYNAMICS/AHORASI/")

#Load Climatic Data
clim<-raster::getData('worldclim', var='bio', res=2.5)[[c(1,2,4,10,11,16:18)]] 

# Load Occs Data
occs<-read.csv("C:/Users/Lenovo/Dropbox/PROJECTS_VIENNA/CLIMATIC DYNAMICS/FOR SUBMISSION/GCB_Submission/DATA&CODES/All_Occs.csv", sep=";")
occs.nat<- occs[occs$Status=="NATIVE",]
occs.alien<- occs[occs$Status=="ALIEN",]

# Keep only clusts with at least 5 records
All.clusts<-count(occs.alien$Cluster)
Plus5Clusts<-All.clusts[All.clusts$freq>4, 1]
occs.alienOK<-occs.alien[occs.alien$Cluster %in%Plus5Clusts, ] 
clust<-unique(occs.alienOK$Cluster) ###Define alien clusters

##Add to each cluster the respective species native records
data.nat.alien<-pblapply(1:length(clust), function(x){
  alien.clust<- occs.alienOK[occs.alienOK$Cluster==clust[x],]
  native.occs<-occs.nat[occs.nat$Species==unique(alien.clust$Species),]
  spX.clustX<-rbind(alien.clust, native.occs)
})

##COUE analysis
#Create folder to save resulting plots
dir.create("RESULTS/COUE/PlotsVarContrib")
dir.create("RESULTS/COUE/DensityKernells")

#1.Extract climatic niche for all points
NicheDynamics<-pblapply(sample(1:10), function(x) #run example for 10 
{
  EnvAllPts<-raster::extract(clim,data.nat.alien[[x]][,4:5])
  EnvAllPts<-cbind(data.nat.alien[[x]],EnvAllPts)
  EnvAllPts<-na.omit(EnvAllPts)
  
  #2.Extract climatic data for NATIVE POINTS
  NatOccs<-data.nat.alien[[x]][data.nat.alien[[x]]$Status=="NATIVE",]
  EnvNatPts<-raster::extract(clim,NatOccs[,4:5])
  EnvNatPts<-cbind(NatOccs,EnvNatPts)
  EnvNatPts<-na.omit(EnvNatPts)
  
  #3. Extract climatic data for ALIEN POINTS
  InvOccs<-data.nat.alien[[x]][data.nat.alien[[x]]$Status=="ALIEN",]
  EnvInvPts<-raster::extract(clim,InvOccs[,4:5])
  EnvInvPts<-cbind(InvOccs,EnvInvPts)
  EnvInvPts<-na.omit(EnvInvPts)
  
  #4.Extract climatic data for Native BG
  dist=50000
  NatOccs<-data.nat.alien[[x]][data.nat.alien[[x]]$Status=="NATIVE",]  
  Natpnts_SF <- st_as_sf(x = NatOccs[,4:5], coords = c("decimalLongitude",
                                                       "decimalLatitude"), crs = crs(clim[[1]]))
  
  chNat <- st_convex_hull(st_union(Natpnts_SF))
  buffs.Nat<-chNat %>% st_buffer(dist)
  buffs.Nat<-as(buffs.Nat, Class = "Spatial")
  Nat.BG<-crop(clim,buffs.Nat)
  Nat.BG<-raster::mask(Nat.BG,buffs.Nat)
  Nat.BG<-crop(Nat.BG,buffs.Nat)
  Nat.BG<-na.omit(as.data.frame(Nat.BG, xy=T))
  
  #5.Extract climatic data for Native BG
  InvOccs<-data.nat.alien[[x]][data.nat.alien[[x]]$Status=="ALIEN",]  
  Inv.pnts_SF <- st_as_sf(x = InvOccs[,4:5], coords = c("decimalLongitude", "decimalLatitude"),
                          crs = crs(clim[[1]]))
  chInv <- st_convex_hull(st_union(Inv.pnts_SF)) 
  buffs.Inv<-chInv %>%st_buffer(dist)%>% st_union()
  buffs.Inv<-as_Spatial(buffs.Inv)
  Inv.BG<-crop(clim,buffs.Inv)
  Inv.BG<-raster::mask(Inv.BG,buffs.Inv)
  Inv.BG<-crop(Inv.BG,buffs.Inv)
  Inv.BG<-na.omit(as.data.frame(Inv.BG, xy=T))
  
  #6. Climatic data all BG
  All.BG<-rbind(Inv.BG,Nat.BG)
  
  #7. PCA-ENVIRONMENT
  pca.env <- dudi.pca(All.BG[,3:10],scannf=F,nf=2)
  
  ##Plot and save variable contributions
  mypath <- file.path("RESULTS/COUE/PlotsVarContrib/",
                      paste("Contributions", unique(data.nat.alien[[x]]$Spp_Country),
                            ".pdf", sep = ""))
  
  pdf(file=mypath)
  ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
  dev.off()
  
  ##PCA scores for the whole study area
  scores.globclim <- pca.env$li
  
  ##PCA scores for the species native points
  scores.sp.nat <-suprow(pca.env,EnvNatPts[,8:15])$li
  
  # PCA scores for the species alien points
  scores.sp.inv <-suprow(pca.env,EnvInvPts[,8:15])$li
  
  # PCA scores for the whole native study area
  scores.clim.nat <- suprow(pca.env,Nat.BG[, 3:10])$li
  
  # PCA scores for the whole invaded study area
  scores.clim.inv <- suprow(pca.env,Inv.BG[, 3:10])$li
  
  ####Calculate the Occurrence Densities Grid with ecospat.grid.clim.dyn()
  # gridding the native niche
  grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.nat,
                                         sp= scores.sp.nat, R=100,th.sp=0.1)
  
  # gridding the invasive niche
  grid.clim.inv <-  ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.inv,
                                          sp=scores.sp.inv, R=100,th.sp=0.1)
  
  #Calculate Niche Overlap with ecospat.niche.overlap()
  # Compute Schoener's D, index of niche overlap
  D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D
  
  ##Delimiting niche categories and quantifying niche dynamics in analogue climates
  ##with ecospat.niche.dyn.index()
  niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = 0.1)
  expansion<-niche.dyn$dynamic.index.w[1]
  stability<-niche.dyn$dynamic.index.w[2]
  unfilling<-niche.dyn$dynamic.index.w[3]
  metrics<-c(unique(InvOccs$Species), unique(InvOccs$Country),D.overlap,expansion, stability, unfilling)
  names(metrics)<-c("species","cluster","D.overlap","expansion", "stability", "unfilling")
  
  ###Kernell Plots Overlap
  mypath2 <- file.path("RESULTS/COUE/DensityKernells/",paste("Niche Overlap", unique(InvOccs$Species),"-", unique(InvOccs$Country), ".pdf", sep = " "))
  
  pdf(file=mypath2, width = 10, height = 10)
  ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.2, interest=1,
                         title= paste("Niche Overlap", unique(InvOccs$Species), "in",
                                      unique(InvOccs$Country)), name.axis1="PC1",
                         name.axis2="PC2")
  dev.off()
  
  return(metrics)
})

Output_metrics<-do.call(rbind.data.frame, NicheDynamics)
names(Output_metrics)<-c("species","country","d.overlap", "expansion", "stability", "unfilling")

print(Output_metrics)
