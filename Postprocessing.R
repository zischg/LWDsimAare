#===============================================================================

# POSTPROCESSING

# for detailed instructions please see the user manual: "LWDsimR: Simulation of Woody Debris Dynamics during floods"

#===============================================================================

# define working directory (Folder "Transportmodel")
wd<-getwd()
setwd(wd)

# define data directory (Folder "Project")
dd<-gsub("Transportmodel","Project",wd)
# ------------------------------------------------------------------

#load packages
library(rgdal)
library(raster)
library(RColorBrewer)
library(pracma)

#define projection (default = CH1903/LV03)
kbs<-CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs")

#define color ramp
colr<-brewer.pal(9,"Pastel1")

#load data
load(paste(dd,"/Results/Results.RData",sep="")) #loading the results may changing the working directory. 
                                                #Please assure that the working directory is correct before proceeding.

s.SB<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="Systemboundary_example")
s.ES<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="vegetation_mask_example")
#-------------------------------------------------------------------

#Identify trees
v.move<-which(a.res[,"xcoord",dim(a.res)[3]-1]!=a.res[,"xcoord",1])                              # ID of all mobilized trees
v.move_GH<-which(a.res[,"xcoord",dim(a.res)[3]-1]!=a.res[,"xcoord",1] & a.res[,"Status",1]==1)   # ID of all mobilized greenwood
v.move_TH<-which(a.res[,"xcoord",dim(a.res)[3]-1]!=a.res[,"xcoord",1] & a.res[,"Status",1]==2)   # ID of all mobilized dead wood
v.nomove<-which(a.res[,"xcoord",dim(a.res)[3]-1]==a.res[,"xcoord",1])                            # ID of trees which was not transported
v.dep<-which(a.res[,"xcoord",dim(a.res)[3]-1]!=a.res[,"xcoord",1] &                              # ID of all deposited trees
               a.res[,"TS_out",dim(a.res)[3]-1]==0 &
               a.res[,"TS_jam",dim(a.res)[3]-1]==0 &
               a.res[,"Status",dim(a.res)[3]-1]==2 )
v.vkl<-which(a.res[,"TS_jam",dim(a.res)[3]-1]>0)                                                 # ID of all entrapped trees
v.out<-which(a.res[,"TS_out",dim(a.res)[3]-1]>0 & a.res[,"TS_out",dim(a.res)[3]-1]!=1)           # ID of trees at LSB
v.Wood_vol<-d.Wood$DBH^2*d.Wood$Length*0.4                                                       # Volume of the trees
#-------------------------------------------------------------------

#Write meta information and volumes
m.volume<-matrix(ncol=3,nrow=14,data=NA)
rownames(m.volume)<-c(" ","Total Area          ","Share of Forest Area ","Total wood stock    "," ","Total volume of water","Duration of event"," ","Mobilized wood total","Mobilized living wood","Mobilized dead wood","Deposited Wood       ","Jammed Wood       ","Wood on LSB        ")
m.volume[1,1]<-"              "
m.volume[1,2]<-"Amount"
m.volume[1,3]<-"Unit"
m.volume[2,1]<- round(s.SB@polygons[[1]]@Polygons[[1]]@area/10000)
m.volume[2,2]<-"ha"
m.volume[3,1]<-round((sum(sapply(1:length(s.ES),function(x) s.ES@polygons[[x]]@Polygons[[1]]@area))/10000)/as.numeric(m.volume[2,1])*100)
m.volume[3,2]<-"%"
m.volume[4,1]<-round(sum(v.Wood_vol))
m.volume[4,2]<-"m3"
m.volume[6,1]<-round(trapz(v.Input[,1],v.Input[,2]))
m.volume[6,2]<-"m3"
m.volume[7,1]<-round(v.Input[length(v.Input[,1]),1])/60
m.volume[7,2]<-"min"
m.volume[9,1]<-round(sum(v.Wood_vol[v.move]))
m.volume[9,2]<-"m3"
m.volume[10,1]<-round(sum(v.Wood_vol[v.move_GH]))
m.volume[10,2]<-"m3"
m.volume[11,1]<-round(sum(v.Wood_vol[v.move_TH]))
m.volume[11,2]<-"m3"
m.volume[12,1]<-round(sum(v.Wood_vol[v.dep]))
m.volume[12,2]<-"m3"
m.volume[13,1]<-round(sum(v.Wood_vol[v.vkl]))
m.volume[13,2]<-"m3"
m.volume[14,1]<-round(sum(v.Wood_vol[v.out]))
m.volume[14,2]<-"m3"

write.table(m.volume,paste(dd,"/Results/Volume.txt",sep=""),na="",quote=F,sep="\t",dec=".",col.names = F)

#-------------------------------------------------------------------
#Temporal dependency of woody debris dynamics
#-------------------------------------------------------------------

v.volout<-sapply(1:TS,function(x)sum(v.Wood_vol[which(a.res[,"TS_out",dim(a.res)[3]-1]==x)]))
v.volout[1]<-0
v.volmobi<-sapply(1:TS,function(x)sum(v.Wood_vol[which(a.res[,"TS_mobi",dim(a.res)[3]-1]==x)]))
v.volmobiGH<-sapply(1:TS,function(x)sum(v.Wood_vol[which(a.res[,"TS_mobi",dim(a.res)[3]-1]==x & a.res[,"Status",1]==1)]))                                            
v.volmobiTH<-sapply(1:TS,function(x)sum(v.Wood_vol[which(a.res[,"TS_mobi",dim(a.res)[3]-1]==x & a.res[,"Status",1]==2)]))                                                          
v.volvkl<-sapply(1:TS,function(x)sum(v.Wood_vol[which(a.res[,"TS_jam",dim(a.res)[3]-1]==x)]))
v.voldep<-sapply(1:(TS*TSWT*(TSW/nsave)),f.voldepos)
v.voldep<-sapply(1:TS, function(x) sum(v.voldep[c((x*TSWT*(TSW/nsave)-(TSWT*(TSW/nsave)-1)):(x*TSWT*(TSW/nsave)))]))


#Plot temporal distribution

pdf(paste(dd,"/Results/Temporal_distribution.pdf",sep=""),width=11.69,height=8.27)

n.length<- length(v.voldep)
n.heigth<- ceiling(max(v.voldep+v.volvkl+v.volout)/20)*20
n.ratio<-  max(v.Input[,2])/n.heigth
  
layout(matrix(c(1,2),2,1,byrow=T))
par(mar=c(0,4.1,4.1,4.1),cex=1,cex.axis=1,cex.lab=1,cex.main=1.3)

barplot(rbind(v.volmobiGH,v.volmobiTH),col=c(colr[3],colr[5]),xlim=c(0,n.length*1.2),ylim=c(0,n.heigth),border=NA,main="Temporal distribution of the woody debris dynamics")
if(exists("v.Input")){
  points(c(0:n.length)*1.2,v.Input[,2]/n.ratio,t="l")
  points(c(0:n.length)*1.2,v.Output[,2]/n.ratio,t="l",lty=2)
}
axis(side=4,at=seq(0,n.heigth,20),labels=paste(round(seq(0,n.heigth,20)*n.ratio,digits=0)))
mtext("Discharge [m3]",side=4,line=2.5,adj=0.5)

par(mar=c(5.1,4.1,0,4.1))
barplot(rbind(-v.voldep,-v.volvkl,-v.volout),col=c(colr[2],colr[4],colr[1]),xlim=c(0,n.length*1.2),ylim=c(-n.heigth,0),border=NA,xlab="Time [min]")
axis(side=1,at=seq(0,n.length,20)*1.2,labels=paste(seq(0,n.length,20)))
mtext("Volume of wood [m3]",side=2,line=2.5,adj=1.3)
abline(h=0)

legend(x="bottomleft",
       xpd=T,ncol=1,
       cex=1.2,
       legend=c("Mobilized living wood","Mobilized dead wood","Deposited wood","Jammed wood","Wood reached LSB","Input Hydrograph","Output Hydrograph"),
       pch=c(15,15,15,15,15,NA,NA),
       lty=c(NA,NA,NA,NA,NA,1,2),
       pt.cex=2,
       col=c(colr[3],colr[5],colr[2],colr[4],colr[1],"black","black"))

dev.off()

#-------------------------------------------------------------------
#Spatial dependency of woody debris dynamics
#-------------------------------------------------------------------

#Write vegetation Shapefile at the end of the modelrun
d.res_tot<-data.frame(a.res[,,dim(a.res)[3]-1])
colnames(d.res_tot)<-c("xcoord2","ycoord2","Status2","TS_out2","Obs_Nr2","TS_in2","TS_mobi2","TS_jam2")
d.res_tot<-cbind(d.Wood,d.res_tot)
s.spatial_points<-SpatialPointsDataFrame(a.res[,c("xcoord","ycoord"),dim(a.res)[3]-1],d.res_tot,proj4string = kbs) 
writeOGR(s.spatial_points,paste(dd,"/Results",sep=""),layer="trees_end",driver="ESRI Shapefile", overwrite_layer = T)

#-------------------------------------------------------------------
#Calculate raster of woody debris dynamic
#Define size of raster
minx<-trunc(min(a.res[,"xcoord",],na.rm=T))
maxx<-trunc(max(a.res[,"xcoord",],na.rm=T)+1)
miny<-trunc(min(a.res[,"ycoord",],na.rm=T))
maxy<-trunc(max(a.res[,"ycoord",],na.rm=T)+1)

m.defraster<-matrix(ncol=length(seq(minx,maxx,1)),nrow=length(seq(miny,maxy,1)),data=0)
colnames(m.defraster)<-seq(minx,maxx,1)
rownames(m.defraster)<-rev(seq(miny,maxy,1))
#-------------------------------------------------------------------
#Flow path

agfac<-2 #define factor of aggrigation

#Amount of trees
m.path<-m.defraster
for(i in 1:(dim(a.res)[3]-2)){
  treeID<-which(a.res[,"Status",i]==3)
  if(length(treeID)>0){
    m.treecords<-trunc(a.res[treeID,c("ycoord","xcoord"),i])
    if(is.matrix(m.treecords)){
      for(j in 1:length(m.treecords[,1])){
        m.path[paste(m.treecords[j,1]),paste(m.treecords[j,2])]<-m.path[paste(m.treecords[j,1]),paste(m.treecords[j,2])]+1
      }} else {m.path[paste(m.treecords[1]),paste(m.treecords[2])]<-m.path[paste(m.treecords[1]),paste(m.treecords[2])]+1}
  }
}
r.path<-raster(m.path)
projection(r.path)<-kbs
extent(r.path)<-c(minx,maxx,miny,maxy)
r.path<-aggregate(r.path,fact=agfac,fun=sum)

writeRaster(r.path,paste(dd,"/Results/Flowpath_amount.tif",sep=""),format="GTiff",overwrite=T)


#Volume
m.path<-m.defraster
for(i in 1:(dim(a.res)[3]-2)){
  treeID<-which(a.res[,"Status",i]==3)
  if(length(treeID)>0){
    m.treecords<-trunc(a.res[treeID,c("ycoord","xcoord"),i])
    if(is.matrix(m.treecords)){
      for(j in 1:length(m.treecords[,1])){
        m.path[paste(m.treecords[j,1]),paste(m.treecords[j,2])]<-m.path[paste(m.treecords[j,1]),paste(m.treecords[j,2])]+v.Wood_vol[as.numeric(rownames(m.treecords)[j])]
      }} else {m.path[paste(m.treecords[1]),paste(m.treecords[2])]<-m.path[paste(m.treecords[1]),paste(m.treecords[2])]+v.Wood_vol[treeID]}
  }
}
r.path<-raster(m.path)
projection(r.path)<-kbs
extent(r.path)<-c(minx,maxx,miny,maxy)
r.path<-aggregate(r.path,fact=agfac,fun=sum)

writeRaster(r.path,paste(dd,"/Results/Flowpath_volume.tif",sep=""),format="GTiff",overwrite=T)
#-------------------------------------------------------------------

#Mobilized volume

agfac<-10 #define factor of aggrigation

m.mobisum<-m.defraster

for(i in 1:length(v.move)){
  m.mobisum[paste(trunc(a.res[v.move[i],"ycoord",1])),paste(trunc(a.res[v.move[i],"xcoord",1]))]<-m.mobisum[paste(trunc(a.res[v.move[i],"ycoord",1])),paste(trunc(a.res[v.move[i],"xcoord",1]))]+v.Wood_vol[v.move[i]]
}

r.mobisum<-raster(m.mobisum)
projection(r.mobisum)<-kbs
extent(r.mobisum)<-c(minx,maxx,miny,maxy)

r.mobisum<-aggregate(r.mobisum,fact=agfac,fun=sum)

writeRaster(r.mobisum,paste(dd,"/Results/Mobilisation_volume.tif",sep=""),format="GTiff",overwrite=T)
#-------------------------------------------------------------------

#Mobilized volume of wood on lower system boundary

agfac<-10 #define factor of aggrigation

m.outsum<-m.defraster

for(i in 1:length(v.out)){
  m.outsum[paste(trunc(a.res[v.out[i],"ycoord",1])),paste(trunc(a.res[v.out[i],"xcoord",1]))]<-m.outsum[paste(trunc(a.res[v.out[i],"ycoord",1])),paste(trunc(a.res[v.out[i],"xcoord",1]))]+v.Wood_vol[v.out[i]]
}


r.outsum<-raster(m.outsum)
projection(r.outsum)<-kbs
extent(r.outsum)<-c(minx,maxx,miny,maxy)

r.outsum<-aggregate(r.outsum,fact=agfac,fun=sum)

writeRaster(r.outsum,paste(dd,"/Results/MobilisationLSB_volume.tif",sep=""),format="GTiff",overwrite=T)
#-------------------------------------------------------------------

#Deposited volume

agfac<-10 #define factor of aggrigation

m.depsum<-m.defraster

for(i in 1:length(v.dep)){
  m.depsum[paste(trunc(a.res[v.dep[i],"ycoord",dim(a.res)[3]-2])),paste(trunc(a.res[v.dep[i],"xcoord",dim(a.res)[3]-2]))]<-m.depsum[paste(trunc(a.res[v.dep[i],"ycoord",dim(a.res)[3]-2])),paste(trunc(a.res[v.dep[i],"xcoord",dim(a.res)[3]-2]))]+v.Wood_vol[v.dep[i]]
}


r.depsum<-raster(m.depsum)
projection(r.depsum)<-kbs
extent(r.depsum)<-c(minx,maxx,miny,maxy)

r.depsum<-aggregate(r.depsum,fact=agfac,fun=sum)

writeRaster(r.depsum,paste(dd,"/Results/Deposition_volume.tif",sep=""),format="GTiff",overwrite=T)
#-------------------------------------------------------------------

#Jammed volume

agfac<-10 #define factor of aggrigation

m.vklsum<-m.defraster

for(i in 1:length(v.vkl)){
  m.vklsum[paste(trunc(a.res[v.vkl[i],"ycoord",dim(a.res)[3]-1])),paste(trunc(a.res[v.vkl[i],"xcoord",dim(a.res)[3]-1]))]<-m.vklsum[paste(trunc(a.res[v.vkl[i],"ycoord",dim(a.res)[3]-1])),paste(trunc(a.res[v.vkl[i],"xcoord",dim(a.res)[3]-1]))]+v.Wood_vol[v.vkl[i]]
}


r.vklsum<-raster(m.vklsum)
projection(r.vklsum)<-kbs
extent(r.vklsum)<-c(minx,maxx,miny,maxy)

r.vklsum<-aggregate(r.vklsum,fact=agfac,fun=sum)

writeRaster(r.vklsum,paste(dd,"/Results/Jammed_volume.tif",sep=""),format="GTiff",overwrite=T)


