#===============================================================================

# PREPROCESSING

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
library(FNN)
library(rgeos)

#define projection (default = CH1903/LV03)
kbs<-CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +units=m +no_defs")
# ------------------------------------------------------------------

#load data

#DTM and DSM
r.dtm<-raster(paste(dd,"/Preprocessing/DTM_example.tif",sep="")) #load DTM (example data: (c)Amt für Geoinformation des Kantons Bern)
r.dom<-raster(paste(dd,"/Preprocessing/DOM_example.tif",sep="")) #load DSM (example data: (c)Amt für Geoinformation des Kantons Bern)
crs(r.dtm)<-kbs
crs(r.dom)<-kbs

#calculate canopy layer (subtract DTM from DSM)
r.canopy<-r.dom-r.dtm

#calculate slope
r.slope<-terrain(r.dtm,opt="slope", unit="degrees",neighbors=8)

#load vegetation mask "Wald BE". (If possible allready clipped with the System boundary) 
s.ES<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="vegetation_mask_example")
s.ES<-s.ES[-which(s.ES@data$CODE==0),] 

#If vegetation mask "Wald BE" is not available, load alternative vegetation mask (polygon shapefile of the forested area)
#Execute the next line only if "Wald BE" is not available
s.ES<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="vegetation_mask_alternative")

#********************************************************************************
#Extract trees out of canopy layer
#********************************************************************************

#search local maxima
n.window_size<-7    #define size of the moving window [cells]
                    #the number of generated trees should be compared with forest inventories if possible
                    #decrease the size of the window if the number of generated trees is too low and increase it if the number is too high.
n.min<-3            #define minimum heigth of the trees [m]

f.max<-              function(x){max(x,na.rm=T)}
m.window<-           matrix(rep(1,n.window_size^2), nrow=n.window_size)
r.max<-              focal(r.canopy,w=m.window,fun=f.max,pad=T,padValue=T)
r.localmax<-         r.max==r.canopy & r.canopy>n.min
s.trees<-            xyFromCell(r.localmax,Which(r.localmax==1,cells=T),spatial=T)
s.trees<-            SpatialPointsDataFrame(s.trees,as.data.frame(matrix(nrow=length(s.trees))),proj4string = kbs)
projection(s.trees)<-kbs

#clip the extracted local maxima with forest mask
s.trees<-s.trees[c(is.na(s.trees%over%s.ES[,1])==F),]


# write shapefile "trees_temp.shp" as safety copy.
writeOGR(s.trees,paste(dd,"/Preprocessing",sep=""),layer="trees_temp",driver="ESRI Shapefile",overwrite=T)


#********************************************************************************
#Calculate attributes of the trees 
#********************************************************************************

#If the vegetation mask "Bestandesinformationen Wald BE" is not available,
#skip the following section and continue with line 167: "Calculate attribute of trees without "Bestandesinformation Wald BE"

#load trees
s.trees<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="trees_temp")
# ------------------------------------------------------------------

#Calculate vegetation density
v.tree_dichte_ES<-c() #Tree density in each polygon
v.tree_dichte<-c()    #Tree densiy in each polygon for every tree
for(i in 1:length(s.ES)){
  tree_in_polygon<-which(is.na(s.trees%over%s.ES[i,][,2])==F)
  v.tree_dichte_ES[i]<-length(tree_in_polygon)/(sum(sapply(slot(s.ES@polygons[[i]],"Polygons"),slot,"area"))/10000) 
  v.tree_dichte[c(tree_in_polygon)]<-v.tree_dichte_ES[i]
  print(paste(round(i/length(s.ES)*100,digits=1),"%",sep=""))
}

#Calculate slope at place of the tree
v.slope_tree<-extract(r.slope,s.trees,method='bilinear')

#Extract Entwicklungsstufe from the vegetation mask "Wald BE"
v.ES_tree<-(s.trees%over%s.ES)[,1]

#Calculate height of the tree
v.height_tree<-extract(r.canopy,s.trees,method='bilinear')

#Calculate altitude of the trees (optional)
v.elev_tree<-extract(r.dtm,s.trees,method='bilinear')

#Calculate DBH out of the ES
v.BHD<-c()
for(i in 1:length(s.trees)){
  if(v.ES_tree[i]==1){
    v.BHD[i]<-sample(c(5:12),1,replace=T)/100
  } else {if(v.ES_tree[i]==2){
    v.BHD[i]<-sample(c(12:20),1,replace=T)/100
  } else {if(v.ES_tree[i]==3){
    v.BHD[i]<-sample(c(21:30),1,replace=T)/100
  } else {if(v.ES_tree[i]==4){
    v.BHD[i]<-sample(c(31:40),1,replace=T)/100
  } else {if(v.ES_tree[i]==5){
    v.BHD[i]<-sample(c(41:50),1,replace=T)/100
  } else {if(v.ES_tree[i]==6){
    v.BHD[i]<-sample(c(50:100),1,replace=T)/100
  } else {
    v.BHD[i]<-0
  }}}}}}
}
# ------------------------------------------------------------------

#Classification of the vegetation according to Mazzorana et al. (2011)
n.grenze_dichte<-441  #define threshold for dense/fragmentary vegetation [tree/ha]
v.Waldzustand<-sapply(1:length(s.trees),function(x) 
  if(v.tree_dichte[x]>=n.grenze_dichte & v.ES_tree[x]==1){
    1
  } else{if(v.tree_dichte[x]<n.grenze_dichte & v.ES_tree[x]==1){
    2
  } else{if(v.tree_dichte[x]>=n.grenze_dichte & v.ES_tree[x]!=1){
    3
  } else{4}}}
)

#Classification of the slope
n.grenze_slope<-25 #define threshold for flat/steep [°]
v.slope<-sapply(1:length(s.trees),function(x)
  if(v.slope_tree[x]>=25){
    1
  } else {2}
)
# ------------------------------------------------------------------

#Create the vegetation layer

s.trees@data$V1<-            c(1:length(s.trees))  #ID
s.trees@data$"xcoord"<-      s.trees@coords[,1]
s.trees@data$"ycoord"<-      s.trees@coords[,2]
s.trees@data$"DBH"<-         v.BHD
s.trees@data$"Status"<-      1
s.trees@data$"Rootwad"<-     v.BHD*5
s.trees@data$"Length"<-      v.height_tree
s.trees@data$"Structure"<-   v.Waldzustand
s.trees@data$"Slope"<-       v.slope
s.trees@data$"TS_in"<-       1

names(s.trees@data)[1]<-"ID"


#********************************************************************************
#Calculate attribute of trees without "Bestandesinformation Wald BE"
#********************************************************************************

#run this section is only if the "Bestandesinformation Wald BE" is not available.
#Several parameters are set by default and the model complexity is therefore reduced.
#if you allready generated the attributes with the "Bestandesinformation Wald BE" skip this section and continue with line 252

#load trees
s.trees<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="trees_temp")
# ------------------------------------------------------------------

#Calculate slope
v.slope_tree<-extract(r.slope,s.trees,method='bilinear')

#Calculate height of the tree
v.height_tree<-extract(r.canopy,s.trees,method='bilinear')

#Calculate altitude of the trees (optional)
v.elev_tree<-extract(r.dtm_zulg,s.trees,method='bilinear')

#Define decile for DBH distribution
#The suggested deciles are calculatet from samples at the Aare river and can be adapted manually.
v.quantile_BHD_GH<-c(0.16,0.18,0.20,0.23,0.25,0.29,0.35,0.42,0.54)

#Calculate DBH out of the deciles
v.10tel<-seq(1,length(s.trees),round(length(s.trees)/10))[-1]
n.diff<-v.10tel[2]-v.10tel[1]
n.min<-0.15 #define minimum DBH [m]
n.max<-1    #define maximum DBH [m]

v.BHD<-c(1:length(s.trees))
for(i in 1:length(s.trees)){
  if(i<=v.10tel[1]){
    v.BHD[c(1:v.10tel[1]-1)]<-n.min
  } else{if(i<=v.10tel[2]){
    v.BHD[c(v.10tel[1]:v.10tel[2]-1)]<- sample(c((v.quantile_BHD_GH[1]*100):(v.quantile_BHD_GH[2]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[3]){
    v.BHD[c(v.10tel[2]:v.10tel[3]-1)]<- sample(c((v.quantile_BHD_GH[2]*100):(v.quantile_BHD_GH[3]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[4]){
    v.BHD[c(v.10tel[3]:v.10tel[4]-1)]<-  sample(c((v.quantile_BHD_GH[3]*100):(v.quantile_BHD_GH[4]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[5]){
    v.BHD[c(v.10tel[4]:v.10tel[5]-1)]<- sample(c((v.quantile_BHD_GH[4]*100):(v.quantile_BHD_GH[5]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[6]){
    v.BHD[c(v.10tel[5]:v.10tel[6]-1)]<-  sample(c((v.quantile_BHD_GH[5]*100):(v.quantile_BHD_GH[6]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[7]){
    v.BHD[c(v.10tel[6]:v.10tel[7]-1)]<- sample(c((v.quantile_BHD_GH[6]*100):(v.quantile_BHD_GH[7]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[8]){
    v.BHD[c(v.10tel[7]:v.10tel[8]-1)]<- sample(c((v.quantile_BHD_GH[7]*100):(v.quantile_BHD_GH[8]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[9]){
    v.BHD[c(v.10tel[8]:v.10tel[9]-1)]<- sample(c((v.quantile_BHD_GH[8]*100):(v.quantile_BHD_GH[9]*100)),n.diff,replace=T)/100
  } else{if(i>v.10tel[9]){
    v.BHD[c(v.10tel[9]:length(s.trees))]<- sample(c((v.quantile_BHD_GH[9]*100):(n.max*100)),n.diff-1,replace=T)/100
  } else{}}}}}}}}}}
}

#Classification of the slope
n.grenze_slope<-25 #define threshold for flat/steep [°]
v.slope<-sapply(1:length(s.trees),function(x)
  if(v.slope_tree[x]>=25){
    1
  } else {2}
)

#Classification of the vegetation according to Mazzorana et al. (2011)
woodstr<-3 #1=young+dense, 2= young+fragmentary, 3=even aged+dense, 4=even aged + fragmentary
v.Waldzustand<-rep(woodstr,length(s.trees))
# ------------------------------------------------------------------

#Create the vegetation layer

s.trees@data$V1<-            c(1:length(s.trees))  #ID
s.trees@data$"xcoord"<-      s.trees@coords[,1]
s.trees@data$"ycoord"<-      s.trees@coords[,2]
s.trees@data$"DBH"<-         v.BHD
s.trees@data$"Status"<-      1
s.trees@data$"Rootwad"<-     v.BHD*5
s.trees@data$"Length"<-      v.height_tree
s.trees@data$"Structure"<-   v.Waldzustand
s.trees@data$"Slope"<-       v.slope
s.trees@data$"TS_in"<-       1

names(s.trees@data)[1]<-"ID"


#********************************************************************************
#Add dead wood to the vegetation layer
#********************************************************************************

#The following section enables adding dead wood to the layer
#If no dead wood should be considered, skip this section an continue with line 350

#load shapefile of dead wood. The points can be generated randomly with the requested density in QGIS
s.TH<-readOGR(paste(dd,"/Preprocessing",sep=""),layer="dead_wood_example")

#Define decile for length and DBH distribution
#The suggested deciles are calculated from samples at the Aare river and can be adapted.
v.quantile_BHD<-c(0.150,0.160,0.170,0.182,0.200,0.230,0.280,0.350,0.490)
v.quantile_length<-c(3.500,4.500,5.000,6.000,7.890,9.000,10.100,13.532,18.360)

#Calculate DBH out of the deciles
v.10tel<-seq(1,length(s.TH),round(length(s.TH)/10))[-1]
n.diff<-v.10tel[2]-v.10tel[1]
n.min<-0.15 #define minimum DBH [m]
n.max<-0.8  #define maximum DBH [m] 

v.TH_BHD<-c(1:length(s.TH))
for(i in 1:length(s.TH)){
  if(i<=v.10tel[1]){
    v.TH_BHD[c(1:v.10tel[1]-1)]<-n.min
  } else{if(i<=v.10tel[2]){
    v.TH_BHD[c(v.10tel[1]:v.10tel[2]-1)]<- sample(c((v.quantile_BHD[1]*100):(v.quantile_BHD[2]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[3]){
    v.TH_BHD[c(v.10tel[2]:v.10tel[3]-1)]<- sample(c((v.quantile_BHD[2]*100):(v.quantile_BHD[3]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[4]){
    v.TH_BHD[c(v.10tel[3]:v.10tel[4]-1)]<-  sample(c((v.quantile_BHD[3]*100):(v.quantile_BHD[4]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[5]){
    v.TH_BHD[c(v.10tel[4]:v.10tel[5]-1)]<- sample(c((v.quantile_BHD[4]*100):(v.quantile_BHD[5]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[6]){
    v.TH_BHD[c(v.10tel[5]:v.10tel[6]-1)]<-  sample(c((v.quantile_BHD[5]*100):(v.quantile_BHD[6]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[7]){
    v.TH_BHD[c(v.10tel[6]:v.10tel[7]-1)]<- sample(c((v.quantile_BHD[6]*100):(v.quantile_BHD[7]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[8]){
    v.TH_BHD[c(v.10tel[7]:v.10tel[8]-1)]<- sample(c((v.quantile_BHD[7]*100):(v.quantile_BHD[8]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[9]){
    v.TH_BHD[c(v.10tel[8]:v.10tel[9]-1)]<- sample(c((v.quantile_BHD[8]*100):(v.quantile_BHD[9]*100)),n.diff,replace=T)/100
  } else{if(i>v.10tel[9]){
    v.TH_BHD[c(v.10tel[9]:length(s.TH))]<- sample(c((v.quantile_BHD[9]*100):(n.max*100)),n.diff-1,replace=T)/100
  } else{}}}}}}}}}}
}

#Calculate length out of the deciles
n.min<-3  #define minimum length [m]
n.max<-35 #define maximum length [m]

v.TH_length<-c(1:length(s.TH))
for(i in 1:length(s.TH)){
  if(i<=v.10tel[1]){
    v.TH_length[c(1:v.10tel[1]-1)]<-sample(c((n.min*100):(v.quantile_length[1]*100)),n.diff-1,replace=T)/100
  } else{if(i<=v.10tel[2]){
    v.TH_length[c(v.10tel[1]:v.10tel[2]-1)]<- sample(c((v.quantile_length[1]*100):(v.quantile_length[2]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[3]){
    v.TH_length[c(v.10tel[2]:v.10tel[3]-1)]<- sample(c((v.quantile_length[2]*100):(v.quantile_length[3]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[4]){
    v.TH_length[c(v.10tel[3]:v.10tel[4]-1)]<-  sample(c((v.quantile_length[3]*100):(v.quantile_length[4]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[5]){
    v.TH_length[c(v.10tel[4]:v.10tel[5]-1)]<- sample(c((v.quantile_length[4]*100):(v.quantile_length[5]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[6]){
    v.TH_length[c(v.10tel[5]:v.10tel[6]-1)]<-  sample(c((v.quantile_length[5]*100):(v.quantile_length[6]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[7]){
    v.TH_length[c(v.10tel[6]:v.10tel[7]-1)]<- sample(c((v.quantile_length[6]*100):(v.quantile_length[7]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[8]){
    v.TH_length[c(v.10tel[7]:v.10tel[8]-1)]<- sample(c((v.quantile_length[7]*100):(v.quantile_length[8]*100)),n.diff,replace=T)/100
  } else{if(i<=v.10tel[9]){
    v.TH_length[c(v.10tel[8]:v.10tel[9]-1)]<- sample(c((v.quantile_length[8]*100):(v.quantile_length[9]*100)),n.diff,replace=T)/100
  } else{if(i>v.10tel[9]){
    v.TH_length[c(v.10tel[9]:length(s.TH))]<- sample(c((v.quantile_length[9]*100):(n.max*100)),n.diff-1,replace=T)/100
  } else{}}}}}}}}}}
}

#Add rootwad. Only 10% of the sample dead wood had rootwads. The percentage can be adapted
perc<-10 #% of dead wood on which a rootwad should be added
rdiam<-5 #size of the rootwad in relation to the DBH
v.sample<-sample(c(1:length(s.TH)),trunc(length(s.TH)/100*perc))
v.root<-c(1:length(s.TH))
v.root[v.sample]<-v.TH_BHD[v.sample]*rdiam
v.root[-v.sample]<-NA
# ------------------------------------------------------------------

#Create the dead wood layer

s.TH@data$ID<-              c((length(s.trees)+1):((length(s.trees))+length(s.TH)))
s.TH@data$"xcoord"<-        s.TH@coords[,1]
s.TH@data$"ycoord"<-        s.TH@coords[,2]
s.TH@data$"DBH"<-           v.TH_BHD
s.TH@data$"Status"<-        2
s.TH@data$"Rootwad"<-       v.root
s.TH@data$"Length"<-        v.TH_length
s.TH@data$"Structure"<-     1 
s.TH@data$"Slope"<-         1 
s.TH@data$"TS_in"<-         1


#********************************************************************************
#Complete vegetation layer
#********************************************************************************

#Combine living and dead wood. 
#Only run this line when dead wood was generated
s.trees<-rbind(s.trees,s.TH)

#Export and save the final vegetation layer in folder "Inputdata" as "trees.shp"
writeOGR(s.trees,paste(dd,"/Inputdata",sep=""),layer="trees",driver="ESRI Shapefile",overwrite=T)

