#===============================================================================

# INPUT-DATEN

#===============================================================================

# Load Mesh Geometry

# Read Mesh as data.frame 
d.Mesh <-           read.table(paste(dd,"/Hydrodynamics/",name2dm,".2dm",sep=""),sep = "",skip=1,header=F,dec=".", fill=T)
n.Elnum<-           which(d.Mesh[,2]==1)[2]-1

# Elemente und Nods separated as data.frame 
d.Elements<-        read.table(paste(dd,"/Hydrodynamics/",name2dm,".2dm",sep=""),sep = "",skip=1,header=F,dec=".", fill=T,nrows=n.Elnum,col.names=c("Element", "ID", "Nod1", "Nod2", "Nod3", "MatID"))
d.Nods<-            read.table(paste(dd,"/Hydrodynamics/",name2dm,".2dm",sep=""),sep = "",skip=n.Elnum+1,header=F,dec=".", fill=T,col.names=c("Nod", "ID", "xcoord", "ycoord", "zcoord"))
n.nds<-             length(d.Nods[,1]) #Number of Nodes

#s.Nods<-SpatialPointsDataFrame(cbind(d.Nods$xcoord,d.Nods$ycoord),d.Nods,proj4string = kbs)
#writeOGR(s.Nods,"D:/Dokumente/MobiLab/Aare/Project/Inputdata","Nods",driver="ESRI Shapefile")

# Define lower system boundary LSB
v.outnodes<-        f.read_outnodes(x)$outnodes

# Read C-Values for mobilisation
d.C_value<-         read.csv(paste(dd,"/Inputdata/C_Werte.csv",sep=""), dec=".",sep=";",header=T)


# Load Obstacles for entrapment and log jam
if(file.exists(paste(dd,"/Inputdata/Obstacle.shp",sep=""))) {
  l.Obstacle<-        f.obstacle(x)
  m.Obstacles_Nodes<- l.Obstacle$Nodes
  m.Obstacles_Info<-  l.Obstacle$Info
  kbs<-               l.Obstacle$kbs
}
# Input und Output-hydrograph from Basement
if(file.exists(paste(dd,"/Hydrodynamics/",name,"_bnd_Outflow_th.dat",sep=""))) {
  v.Output<-          read.table(paste(dd,"/Hydrodynamics/",name,"_bnd_Outflow_th.dat",sep=""),sep="\t",dec=".")
  v.Input<-           read.table(paste(dd,"/Hydrodynamics/",name,"_bnd_Inflow_th.dat",sep=""),sep="\t",dec=".")
}
# Load vegetation layer

s.Wood<-            readOGR(paste(dd,"/Inputdata",sep=""),layer="trees")
d.Wood<-            s.Wood@data
colnames(d.Wood)<-  c("ID","xcoord","ycoord","DBH","Status", "Rootwad","Length","Structure","Slope","TS_in")

d.Wood[which(d.Wood$ID>300000),"ID"]<-which(d.Wood$ID>300000) #redefine ID's from ZULG-trees (only for this Aare-simulation)

# ------------------------------------------------------------------

# Load hydraulic results from Basement
#Skip the first timestep since its 0. Depending on the data-strucute the first 3 or 9 Rows must be skipped
if(bvers==2.5){
  dform<-".dat"
  nskp<-3
  nskp2<-3
  vel<-"vel"
} else {
  dform<-".sol"
  nskp<-9
  nskp2<-10
  vel<-"velocity"
}

# Flowdepth
d.nds_depth<-       read.table(paste(dd,"/Hydrodynamics/",name,"_nds_depth",dform,sep=""),sep=",",skip=nskp+n.nds)
d.nds_depth[grep("TS",d.nds_depth[,1]),]<-NA #Replace "TS" with NA
d.nds_depth[,1]<-   suppressWarnings(as.numeric(paste(d.nds_depth[,1])))
v.depth<-           d.nds_depth[,1] #Write flow depth in vector
rm(d.nds_depth)

# Flow velocity
d.nds_vel<-         read.table(paste(dd,"/Hydrodynamics/",name,"_nds_velocity",dform,sep=""),sep="\t",skip=nskp2,fill=T)
d.nds_vel[grep("TS",d.nds_vel[,1]),]<-NA #Replace "TS" with NA
v.velx<-            append(suppressWarnings(as.numeric(paste(d.nds_vel[,1]))),NA,after=0)
v.vely<-            append(d.nds_vel[,2],NA,after=0)
v.velx<-            v.velx[c((n.nds+2):length(v.velx))]
v.vely<-            v.vely[c((n.nds+2):length(v.vely))]
rm(d.nds_vel)

# Number of Timesteps
n.ts<-              sum(is.na(v.depth))-2 

# ------------------------------------------------------------------
#Identify trees in flooded area

v.node_maxdepth<- sapply(1:(dim(d.Nods)[1]),function(x) max(v.depth[seq(x+1,length(v.depth),dim(d.Nods)[1]+1)]))
m.NN_ID<-         get.knnx(d.Nods[,c("xcoord","ycoord")],d.Wood[,c("xcoord","ycoord")],k=3,algorithm = "kd_tree")$nn.index
m.NN_depth<-      cbind(v.node_maxdepth[m.NN_ID[,1]],v.node_maxdepth[m.NN_ID[,2]],v.node_maxdepth[m.NN_ID[,3]])
v.inwater<-       which(apply(m.NN_depth,1,sum)>0)
d.Wood<-          d.Wood[v.inwater,]



#Determine time during trees are exposed to hydrodynamics

v.timeinwater<-sapply(1:dim(d.Wood)[1], function (x) max(cbind(sum(v.depth[seq((m.NN_ID[d.Wood$ID,][x,1]+1),length(v.depth),dim(d.Nods)[1]+1)]>0),
                                                               sum(v.depth[seq((m.NN_ID[d.Wood$ID,][x,2]+1),length(v.depth),dim(d.Nods)[1]+1)]>0),
                                                               sum(v.depth[seq((m.NN_ID[d.Wood$ID,][x,3]+1),length(v.depth),dim(d.Nods)[1]+1)]>0))))



# ------------------------------------------------------------------

# Create array for saving results.
# for each layer one dimension is added.
m.tdefault<-matrix(nrow=length(d.Wood[,1]),ncol=8,data=NA)
colnames(m.tdefault)<-c('xcoord','ycoord','Status','TS_out','Obs_Nr','TS_in','TS_mobi','TS_jam')
m.tdefault[,"TS_in"]<-d.Wood[,"TS_in"]

m.t0<-m.tdefault
m.t0[c(1:length(d.Wood[,1])),1]<-d.Wood[,c("xcoord")]
m.t0[c(1:length(d.Wood[,1])),2]<-d.Wood[,c("ycoord")]
m.t0[c(1:length(d.Wood[,1])),3]<-d.Wood[,c("Status")]
m.t0[c(1:length(d.Wood[,1])),4]<-0
m.t0[c(1:length(d.Wood[,1])),5]<-0
m.t0[c(1:length(d.Wood[,1])),7]<-0
m.t0[c(1:length(d.Wood[,1])),8]<-0


m.t1<-m.t0
m.t2<-m.tdefault


# ------------------------------------------------------------------
#Define vector with time steps to save the matrices into an array
if(nsave>0){
  Tsave<-nsave/TSW
  v.save<-seq(0,TS*TSWT,Tsave)
  if(v.save[2]>1){
    v.save[1]<-1
  } else{v.save<-v.save[-1]}
} else {v.save<-1}

v.write<-seq(0,TS,nwrite)[-1]

a.out<-array(data=NA,dim=c(length(d.Wood[,1]),8,length(v.save)))
dimnames(a.out)<-list(1:dim(a.out)[1],c('xcoord','ycoord','Status','TS_out','Obs_Nr','TS_in','TS_mobi','TS_jam'),1:dim(a.out)[3])  





