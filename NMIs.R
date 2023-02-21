rm(list=ls())
gc()

# Analyses:Estimation of niche dynamics (Niche Margin Index) 
# Paper: Patterns and drivers of climatic niche dynamics in biological invasions of insular tetrapods

Pckgs<-c("raster", "ade4","rgeos","rgdal","ecospat", "ks", "R2jags", "runjags", "ggplot2", 
         "bayesplot", "boa", "mcmcplots", "viridis", "reshape2", "dplyr", "pbapply")

sapply(Pckgs, require, character.only = TRUE)

## Load Functions
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/NMI_function.R")
source("https://raw.githubusercontent.com/ecospat/NMI/master/functions/mve_function.R")

### choose the parameters of the analysis
grain=2.5 # resolution of the climatic data. In Broennimann et al. 2020 NCOM: 10,30,60
envelope="kde" # choise of the envelop methode, "kde" or "mve"
level=99 # level of inclusion of rare climatic conditions for kde and mve. In Broennimann et al. 2020 NCOM: 99,95,90

#setwd("C:/Users/Lenovo/Dropbox/PROJECTS_VIENNA/CLIMATIC DYNAMICS/AHORASI/")
#Load Records
# Load records
occs<-read.csv("C:/Users/Lenovo/Dropbox/PROJECTS_VIENNA/CLIMATIC DYNAMICS/FOR SUBMISSION/GCB_Submission/DATA&CODES/All_Occs.csv", sep=";")

# Generate climatic data at the chosen grain
clim<-raster::getData('worldclim', var='bio', res=2.5)[[c(1,2,4,10,11,16:18)]] #select bio2,4,10:11,16:18
clim<-aggregate(clim,grain/2.5)
clim.df<-na.exclude(getValues(clim))

### run PCA on all background points (PCA-env sensu Broennimann et al 2012)
pca<- dudi.pca(clim.df,scannf = FALSE, nf = 2)

# extract scores of the background in PCA space
bkg.scores<-pca$li

###NICHE MARGIN INDEX
###*****Example for amphibians. Repeat for the other taxa
#Load shapefile of native distributions 
#*We obtained range maps from IUCN and BirdLige, then for each of the species in the study we subset the full shape to keep only those polygons depicting the native range. These polygons are flagged as origin=1 in the case of IUCN range maps.

# import shapefile of native distributions
amph.shps<-readOGR(dsn="DATA/SPATIAL/RANGES/AMPHIBIANS", layer="NativeRngsAmphs")

#Occs categorized
amphs.nat<- occs[occs$Status=="NATIVE" & occs$Class=="Amphibia",]
amphs.alien<- occs[occs$Status=="ALIEN"& occs$Class=="Amphibia",]
clust.amph<-amphs.alien[,c(2,7)]###Define alien clusters for amphibians (total 44)
clust.amph<-unique(clust.amph)


# extract scores of the native niche in PCA space (from all native ranges)
NMIs.amph<-pblapply(1:5, function (x)
{
  sp.shp<- amph.shps[which(amph.shps$binomial==clust.amph$Species[x])]
  extract<-raster::extract(clim,sp.shp)
  extract<- do.call(rbind,extract)
  sp.scores<-na.omit(suprow(pca,extract)$li)
  
  # extract scores of introductions in PCA space (for all clusters)
  alien<-amphs.alien[which(amphs.alien$Cluster==clust.amph$Cluster[x]),]
  alien.scores<-na.omit(suprow(pca,raster::extract(clim,alien[,4:5]))$li)
  
  ### delineate background and niche margins (Run only once as is the same for all spp)
  fhat<-kde(bkg.scores,compute.cont=TRUE)
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=0.0001)
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  bkg.pol=SpatialPolygons(l99)
  bkg.pol<-aggregate(bkg.pol)
  
  fhat<-kde(sp.scores,compute.cont=TRUE) 
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,
                    level=fhat$cont[length(fhat$cont)])
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  sp.pol=SpatialPolygons(l99)
  sp.pol<-aggregate(sp.pol)
  
  ### calculate NMI values
  bkg.grid<-raster(nrows=100,ncols=100,extent(bkg.pol))
  bkg<- rasterize(bkg.pol,bkg.grid) #(Run only once as is the same for all spp)
  #plot(bkg)
  
  bkg.pts<-data.frame(coordinates(bkg))
  coordinates(bkg.pts) <- cbind(bkg.pts$x , bkg.pts$y)
  bkg.NMI<-NMI(foc.pop = bkg.pts,niche=sp.pol)
  
  coordinates(alien.scores) <- cbind(alien.scores$Axis1 , alien.scores$Axis2)
  alien.NMI<-NMI(foc.pop = alien.scores,niche=sp.pol)
  alien.NMI2<-data.frame(alien.NMI$NMI)
  alien.NMI2$cluster<-clust.amph$Cluster[x]
  
  ### plot Fig and Save
  ##Plot and save variable contributions
  mypath <- file.path("RESULTS/NMI/Plots",
                      paste("NMI_", clust.amph$Cluster[x],
                            ".pdf", sep = ""))
  pdf(file=mypath)
  # plot axes
    plot(coordinates(bkg.pts),type="n",xlab="PC1",ylab="PC2",main=clust.amph$Cluster[x])
  
  #plot NMI
  in.bkg<-!is.na(values(bkg)) #pixels in bkg
  pos<-bkg.NMI$NMI>=0 #pixels with positive values
  
  #plot innerness
  inner<-bkg  
  inner[in.bkg&pos]<-bkg.NMI$NMI[in.bkg&pos] 
  inner[inner==1]<-NA
  Pal.inner <- colorRampPalette(c('white','olivedrab4'))
  plot(inner,col=Pal.inner(100),add=T, legend=FALSE)
  
  #plot outerness
  outer<-bkg  
  outer[in.bkg&!pos]<-bkg.NMI$NMI[in.bkg&!pos]
  outer[outer==1]<-NA
  Pal.outer <- colorRampPalette(c('#283350','deepskyblue4','gray89'))
  plot(outer,col=Pal.outer(100),add=T, legend=FALSE, main=clust.amph$Cluster[x])
  
  # accessible area
  plot(bkg.pol,add=T,lty=4) 
  
  # native climatic niche margin
  plot(sp.pol,add=T,lwd=3,border="olivedrab4") 
  
  #plot intros
  points(alien.scores,col="black",pch=20, cex=1.4)
  points(alien.scores,col="red",pch=20,cex=1)
  dev.off()
  return(alien.NMI2)
})

NMIs.amph<-do.call(rbind,NMIs.amph)
NMIs.amph$class<-"Amphibia"
print(NMIs.amph)

######REPTILES
# import shapefile of native distributions (REPTILES)
rept.shps<-readOGR(dsn="DATA/SPATIAL/RANGES/REPTILES", layer="NativeRngsRepts")

###Occs clasified
repts.nat<- occs[occs$Status=="NATIVE" & occs$Class=="Reptilia",]
repts.alien<- occs[occs$Status=="ALIEN"& occs$Class=="Reptilia",]

clust.repts<-repts.alien[,c(2,7)]###Define alien clusters for reptiles (total 36)

spps<-intersect(clust.repts$Species,rept.shps$binomial) ##double check spp with map
clust.repts<-clust.repts[clust.repts$Species%in%spps,] ##clusters for species with maps (35)
clust.repts<-unique(clust.repts)

# extract scores of the background in PCA space
bkg.scores<-pca$li

# extract scores of the native niche in PCA space (from all native ranges)
NMIs.repts<-pblapply(1:5, function (x)
{
  sp.shp<- rept.shps[which(rept.shps$binomial==clust.repts$Species[x]),]
  extract<-raster::extract(clim,sp.shp)
  extract<- do.call(rbind,extract)
  sp.scores<-na.omit(suprow(pca,extract)$li)
  
  # extract scores of introductions in PCA space (for all clusters)
  alien<-repts.alien[which(repts.alien$Cluster==clust.repts$Cluster[x]),]
  alien.scores<-na.omit(suprow(pca,raster::extract(clim,alien[,4:5]))$li)
  
  ### delineate background and niche margins (Run only once as is the same for all spp)
  fhat<-kde(bkg.scores,compute.cont=TRUE)
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=0.0001)
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  bkg.pol=SpatialPolygons(l99)
  bkg.pol<-aggregate(bkg.pol)
  
  fhat<-kde(sp.scores,compute.cont=TRUE) 
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,
                    level=fhat$cont[length(fhat$cont)])
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  sp.pol=SpatialPolygons(l99)
  sp.pol<-aggregate(sp.pol)
  
  ### calculate NMI values
  bkg.grid<-raster(nrows=100,ncols=100,extent(bkg.pol))
  bkg<- rasterize(bkg.pol,bkg.grid) #(Run only once as is the same for all spp)
  #plot(bkg)
  
  bkg.pts<-data.frame(coordinates(bkg))
  coordinates(bkg.pts) <- cbind(bkg.pts$x , bkg.pts$y)
  bkg.NMI<-NMI(foc.pop = bkg.pts,niche=sp.pol)
  
  coordinates(alien.scores) <- cbind(alien.scores$Axis1 , alien.scores$Axis2)
  alien.NMI<-NMI(foc.pop = alien.scores,niche=sp.pol)
  alien.NMI2<-data.frame(alien.NMI$NMI)
  alien.NMI2$cluster<-alien$Cluster
  
  ### plot Fig and Save
  ##Plot and save variable contributions
  mypath <- file.path("RESULTS/NMI/Plots",
                      paste("NMI_", clust.repts$Cluster[x],
                            ".pdf", sep = ""))
  pdf(file=mypath)
  # plot axes
  plot(coordinates(bkg.pts),type="n",xlab="PC1",ylab="PC2",main=clust.repts$Cluster[x])
  
  #plot NMI
  in.bkg<-!is.na(values(bkg)) #pixels in bkg
  pos<-bkg.NMI$NMI>=0 #pixels with positive values
  
  #plot innerness
  inner<-bkg  
  inner[in.bkg&pos]<-bkg.NMI$NMI[in.bkg&pos] 
  inner[inner==1]<-NA
  Pal.inner <- colorRampPalette(c('white','olivedrab4'))
  plot(inner,col=Pal.inner(100),add=T, legend=FALSE)
  
  #plot outerness
  outer<-bkg  
  outer[in.bkg&!pos]<-bkg.NMI$NMI[in.bkg&!pos]
  outer[outer==1]<-NA
  Pal.outer <- colorRampPalette(c('#283350','deepskyblue4','gray89'))
  plot(outer,col=Pal.outer(100),add=T, legend=FALSE, main=clust.repts$Cluster[x])
  
  # accessible area
  plot(bkg.pol,add=T,lty=4) 
  
  # native climatic niche margin
  plot(sp.pol,add=T,lwd=3,border="olivedrab4") 
  
  #plot intros
  points(alien.scores,col="black",pch=20, cex=1.4)
  points(alien.scores,col="red",pch=20,cex=1)
  dev.off()
  return(alien.NMI2)
})

NMIs.repts<-do.call(rbind,NMIs.repts)
NMIs.repts$class<-"Reptilia"
print(NMIs.repts)

######BIRDS
# import shapefile of native distributions (BIRDS)
birds.shps<-readOGR(dsn="DATA/SPATIAL/RANGES/BIRDS", layer="NativeRngsBirds")

#Occs classified
birds.nat<- occs[occs$Status=="NATIVE" & occs$Class=="Aves",]
birds.alien<- occs[occs$Status=="ALIEN"& occs$Class=="Aves",]

birds.spp<-unique(birds.alien$Species) #(19 spp)

spps.maps<-intersect(birds.shps$binomial,birds.spp) #(19 spp with map)
clust.birds<-birds.alien[,c(2,7)] #49 total clusters
clust.birds<-clust.birds[clust.birds$Species%in%spps.maps,] ##49 (all) clusters has map of nat range
clust.birds<-unique(clust.birds)

#extract scores of the background in PCA space
bkg.scores<-pca$li

# extract scores of the native niche in PCA space (from all native ranges)


NMIs.birds<-pblapply(1:5, function (x)
{
  sp.shp<- birds.shps[which(birds.shps$binomial==clust.birds$Species[x]),]
  extract<-terra::extract(clim,sp.shp)
  extract<- do.call(rbind,extract)
  sp.scores<-na.omit(suprow(pca,extract)$li)
  
  # extract scores of introductions in PCA space (for all clusters)
  alien<-birds.alien[which(birds.alien$Cluster==clust.birds$Cluster[x]),]
  alien.scores<-na.omit(suprow(pca,raster::extract(clim,alien[,4:5]))$li)
  
  ### delineate background and niche margins (Run only once as is the same for all spp)
  fhat<-kde(bkg.scores,compute.cont=TRUE)
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,level=0.0001)
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  bkg.pol=SpatialPolygons(l99)
  bkg.pol<-aggregate(bkg.pol)
  
  fhat<-kde(sp.scores,compute.cont=TRUE) 
  c99<-contourLines(fhat$eval.points[[1]],fhat$eval.points[[2]],fhat$estimate,
                    level=fhat$cont[length(fhat$cont)])
  l99<-list()        
  for(k in 1:length(c99))l99[[k]]=Polygons(list(Polygon(c99[[k]][-1])),ID=k)
  sp.pol=SpatialPolygons(l99)
  sp.pol<-aggregate(sp.pol)
  
  ### calculate NMI values
  bkg.grid<-raster(nrows=100,ncols=100,extent(bkg.pol))
  bkg<- rasterize(bkg.pol,bkg.grid) #(Run only once as is the same for all spp)
  #plot(bkg)
  
  bkg.pts<-data.frame(coordinates(bkg))
  coordinates(bkg.pts) <- cbind(bkg.pts$x , bkg.pts$y)
  bkg.NMI<-NMI(foc.pop = bkg.pts,niche=sp.pol)
  
  coordinates(alien.scores) <- cbind(alien.scores$Axis1 , alien.scores$Axis2)
  alien.NMI<-NMI(foc.pop = alien.scores,niche=sp.pol)
  alien.NMI2<-data.frame(alien.NMI$NMI)
  alien.NMI2$cluster<-clust.birds$Cluster[x]
  
  ### plot Fig and Save
  ##Plot and save variable contributions
  mypath <- file.path("RESULTS/NMI/Plots",
                      paste("NMI_", clust.birds$Cluster[x],
                            ".pdf", sep = ""))
  pdf(file=mypath)
  # plot axes
  plot(coordinates(bkg.pts),type="n",xlab="PC1",ylab="PC2",main=clust.birds$Cluster[x])
  
  #plot NMI
  in.bkg<-!is.na(values(bkg)) #pixels in bkg
  pos<-bkg.NMI$NMI>=0 #pixels with positive values
  
  #plot innerness
  inner<-bkg  
  inner[in.bkg&pos]<-bkg.NMI$NMI[in.bkg&pos] 
  inner[inner==1]<-NA
  Pal.inner <- colorRampPalette(c('white','olivedrab4'))
  plot(inner,col=Pal.inner(100),add=T, legend=FALSE)
  
  #plot outerness
  outer<-bkg  
  outer[in.bkg&!pos]<-bkg.NMI$NMI[in.bkg&!pos]
  outer[outer==1]<-NA
  Pal.outer <- colorRampPalette(c('#283350','deepskyblue4','gray89'))
  plot(outer,col=Pal.outer(100),add=T, legend=FALSE, main=clust.birds$Cluster[x])
  
  # accessible area
  plot(bkg.pol,add=T,lty=4) 
  
  # native climatic niche margin
  plot(sp.pol,add=T,lwd=3,border="olivedrab4") 
  
  #plot intros
  points(alien.scores,col="black",pch=20, cex=1.4)
  points(alien.scores,col="red",pch=20,cex=1)
  dev.off()
  return(alien.NMI2)
})

NMIs.birds<-do.call(rbind,NMIs.birds)
NMIs.birds$class<-"Aves"
print(NMIs.birds)

###All estimates
NMIs.all<-rbind(NMIs.amph, NMIs.repts, NMIs.birds)



meanNMI.clusts<- data.frame(NMIs.all %>%
                            group_by(cluster)%>%
                  dplyr::summarize(n = n(),mean_NMI = mean(alien.NMI.NMI, na.rm=T),
                                               sd_NMI = sd(alien.NMI.NMI, na.rm=T)))




