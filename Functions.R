#===============================================================================

# FUNKTIONEN

#===============================================================================

#*******************************************************************************
# Functions to run the model
#*******************************************************************************

# 1.1. Function to determine Status 1 (Mobilisation). Adapted C-Value from Mazzorana(2011)

f.mobilisation<-function(x,v.velx_TS,v.vely_TS,v.depth_TS,v.timeinwater){
  n.c_value<-((v.velx_TS[x]^2+v.vely_TS[x]^2)/20)+v.depth_TS[x]
  n.tw<-v.timeinwater[x]
  if(n.c_value==0){
    1
  } else{if(d.Wood[x,"Structure"]==1){
    if(d.Wood[x,"Slope"]==1){
      if(n.c_value<0.4){
        sample(1:2,1,prob=c(1-(d.C_value[1,"rb_0.4"]/TSWT/n.tw),(d.C_value[1,"rb_0.4"]/TSWT/n.tw)))
      } else {if(n.c_value<1.5){
        sample(1:2,1,prob=c(1-(d.C_value[1,"rb_se1.5"]/TSWT/n.tw),(d.C_value[1,"rb_se1.5"]/TSWT/n.tw)))
      } else {sample(1:2,1,prob=c(1-(d.C_value[1,"rb_be1.5"]/TSWT/n.tw),(d.C_value[1,"rb_be1.5"]/TSWT/n.tw)))}}
    } else{if(n.c_value<0.4){
      sample(1:2,1,prob=c(1-(d.C_value[1,"fp_0.4"]/TSWT/n.tw),(d.C_value[1,"fp_0.4"]/TSWT/n.tw)))
    } else {if(n.c_value<1.5){
      sample(1:2,1,prob=c(1-(d.C_value[1,"fp_se1.5"]/TSWT/n.tw),(d.C_value[1,"fp_se1.5"]/TSWT/n.tw)))
    } else {sample(1:2,1,prob=c(1-(d.C_value[1,"fp_be1.5"]/TSWT/n.tw),(d.C_value[1,"fp_be1.5"]/TSWT/n.tw)))}}}
      } else {  if(d.Wood[x,"Structure"]==2){
    if(d.Wood[x,"Slope"]==1){
      if(n.c_value<0.4){
        sample(1:2,1,prob=c(1-(d.C_value[2,"rb_0.4"]/TSWT/n.tw),(d.C_value[2,"rb_0.4"]/TSWT/n.tw)))
      } else {if(n.c_value<1.5){
        sample(1:2,1,prob=c(1-(d.C_value[2,"rb_se1.5"]/TSWT/n.tw),(d.C_value[2,"rb_se1.5"]/TSWT/n.tw)))
      } else {sample(1:2,1,prob=c(1-(d.C_value[2,"rb_be1.5"]/TSWT/n.tw),(d.C_value[2,"rb_be1.5"]/TSWT/n.tw)))}}
    } else{if(n.c_value<0.4){
      sample(1:2,1,prob=c(1-(d.C_value[2,"fp_0.4"]/TSWT/n.tw),(d.C_value[2,"fp_0.4"]/TSWT/n.tw)))
    } else {if(n.c_value<1.5){
      sample(1:2,1,prob=c(1-(d.C_value[2,"fp_se1.5"]/TSWT/n.tw),(d.C_value[2,"fp_se1.5"]/TSWT/n.tw)))
    } else {sample(1:2,1,prob=c(1-(d.C_value[2,"fp_be1.5"]/TSWT/n.tw),(d.C_value[2,"fp_be1.5"]/TSWT/n.tw)))}}}
      } else { if(d.Wood[x,"Structure"]==3){
    if(d.Wood[x,"Slope"]==1){
      if(n.c_value<0.4){
        sample(1:2,1,prob=c(1-(d.C_value[3,"rb_0.4"]/TSWT/n.tw),(d.C_value[3,"rb_0.4"]/TSWT/n.tw)))
      } else {if(n.c_value<1.5){
        sample(1:2,1,prob=c(1-(d.C_value[3,"rb_se1.5"]/TSWT/n.tw),(d.C_value[3,"rb_se1.5"]/TSWT/n.tw)))
      } else {sample(1:2,1,prob=c(1-(d.C_value[3,"rb_be1.5"]/TSWT/n.tw),(d.C_value[3,"rb_be1.5"]/TSWT/n.tw)))}}
    } else{if(n.c_value<0.4){
      sample(1:2,1,prob=c(1-(d.C_value[3,"fp_0.4"]/TSWT/n.tw),(d.C_value[3,"fp_0.4"]/TSWT/n.tw)))
    } else {if(n.c_value<1.5){
      sample(1:2,1,prob=c(1-(d.C_value[3,"fp_se1.5"]/TSWT/n.tw),(d.C_value[3,"fp_se1.5"]/TSWT/n.tw)))
    } else {sample(1:2,1,prob=c(1-(d.C_value[3,"fp_be1.5"]/TSWT/n.tw),(d.C_value[3,"fp_be1.5"]/TSWT/n.tw)))}}}
      } else{if(d.Wood[x,"Slope"]==1){
    if(n.c_value<0.4){
      sample(1:2,1,prob=c(1-(d.C_value[4,"rb_0.4"]/TSWT/n.tw),(d.C_value[4,"rb_0.4"]/TSWT/n.tw)))
    } else {if(n.c_value<1.5){
      sample(1:2,1,prob=c(1-(d.C_value[4,"rb_se1.5"]/TSWT/n.tw),(d.C_value[4,"rb_se1.5"]/TSWT/n.tw)))
    } else {sample(1:2,1,prob=c(1-(d.C_value[4,"rb_be1.5"]/TSWT/n.tw),(d.C_value[4,"rb_be1.5"]/TSWT/n.tw)))}}
  } else{if(n.c_value<0.4){
    sample(1:2,1,prob=c(1-(d.C_value[4,"fp_0.4"]/TSWT/n.tw),(d.C_value[4,"fp_0.4"]/TSWT/n.tw)))
  } else {if(n.c_value<1.5){
    sample(1:2,1,prob=c(1-(d.C_value[4,"fp_se1.5"]/TSWT/n.tw),(d.C_value[4,"fp_se1.5"]/TSWT/n.tw)))
  } else {sample(1:2,1,prob=c(1-(d.C_value[4,"fp_be1.5"]/TSWT/n.tw),(d.C_value[4,"fp_be1.5"]/TSWT/n.tw)))}}}}}}
  }
}

# 1.2. Function to determine Status 2 and 3 (Transport and Deposition) 
#      It is checked wether there is a rootwad and if the log is floating or sliding/rolling

f.trans_x<-function(x,v.depth_TS,v.velx_TS,v.vely_TS,m.t1,i,j){
  if(is.na(d.Wood[x,"Rootwad"])){
    if(v.depth_TS[x]/d.Wood[x,"DBH"]<dnr){
      if(sqrt((g*mu*d.Wood[x,"DBH"]^2/(2*Cd*v.depth_TS[x]))*(pi-acos(1-(2*v.depth_TS[x]/d.Wood[x,"DBH"]))+0.5*sin(2*acos(1-(2*v.depth_TS[x]/d.Wood[x,"DBH"])))))<sqrt(v.velx_TS[x]^2+v.vely_TS[x]^2)){
        m.t1[x,"xcoord"]+(1-(1-v.depth_TS[x]/d.Wood[x,"DBH"]))*(v.velx_TS[x]*TSW)
      } else {m.t1[x,"xcoord"]}
    } else {m.t1[x,"xcoord"]+(v.velx_TS[x]*TSW)}
  }  else{if(v.depth_TS[x]/d.Wood[x,"DBH"]<dwr){
    m.t1[x,"xcoord"]
  } else {m.t1[x,"xcoord"]+(v.velx_TS[x]*TSW)}
  }
}

f.trans_y<-function(x,v.depth_TS,v.velx_TS,v.vely_TS,m.t1,i,j){
  if(is.na(d.Wood[x,"Rootwad"])){
    if(v.depth_TS[x]/d.Wood[x,"DBH"]<dnr){
      if(sqrt((g*mu*d.Wood[x,"DBH"]^2/(2*Cd*v.depth_TS[x]))*(pi-acos(1-(2*v.depth_TS[x]/d.Wood[x,"DBH"]))+0.5*sin(2*acos(1-(2*v.depth_TS[x]/d.Wood[x,"DBH"])))))<sqrt(v.velx_TS[x]^2+v.vely_TS[x]^2)){
        m.t1[x,"ycoord"]+(1-(1-v.depth_TS[x]/d.Wood[x,"DBH"]))*(v.vely_TS[x]*TSW)
      } else {m.t1[x,"ycoord"]}
    } else {m.t1[x,"ycoord"]+(v.vely_TS[x]*TSW)}
  } else{if(v.depth_TS[x]/d.Wood[x,"DBH"]<dwr){
    m.t1[x,"ycoord"]
  } else {m.t1[x,"ycoord"]+(v.vely_TS[x]*TSW)}
  }
}
  

# 1.3. Function to determine Status 4 (entrapment)

f.verklausung<-function(x,m.NN_ID,v.depth_TS,m.t1,i,j){
  if(any(m.NN_ID[x,]%in%m.Obstacles_Nodes[,"Node"]==T) & vkl==TRUE & v.depth_TS[x]>0 & m.t1[x,"TS_in"]<=j){
    nodenr<-                    m.NN_ID[x,which(m.NN_ID[x,]%in%m.Obstacles_Nodes[,"Node"]==T)[1]]
    obsnr<-                     m.Obstacles_Nodes[which(m.Obstacles_Nodes[,"Node"]==nodenr),"ID"]
    length<-                    d.Wood[x,"Length"]
    Bmin<-                      m.Obstacles_Info[obsnr,"Width"]
    H<-                         m.Obstacles_Info[obsnr,"Height"]
    probability_pier<-          (-1/15)+((2/15)*(length/Bmin))
    if(probability_pier<0) {
      probability_pier<-0
    } 
    
    if(is.na(d.Wood[x,"Rootwad"])){
      lambda<-                  length/Bmin*((v.depth_TS[x]+(d.Wood[x,"DBH"]/2))/H)
      probability_deck<-        (-0.074+0.88*lambda)
      if(probability_deck<0) {
        probability_deck        <-0
      } 
    } else {mu2<-               ((v.depth_TS[x]+(d.Wood[x,"Rootwad"]/2))/H)
            probability_deck<-  (-3.5+2.56*mu2)
            if(probability_deck<0) {
              probability_deck  <-0
            }
    }
    probability<-               probability_pier+probability_deck 
    if(probability<1){
      sample(3:4,1,prob=c(1-probability,probability))
    } else {sample(3:4,1,prob=c(0,1))}
  } else {m.t1[x,"Status"]}
}


# 1.4. Funktion to adapt the Status of a tree 

f.status<-function(x,i,j,m.t1,m.t2){
  if( m.t1[x,"Status"]==2){
    if(m.t1[x,"xcoord"]==m.t2[x,"xcoord"]){
      2
    } else {3} 
  } else{if( m.t1[x,"Status"]==3){
    if(m.t1[x,"xcoord"]==m.t2[x,"xcoord"]){
      2
    } else {3} 
  } else{if(m.t1[x,"Status"]==4){
    4
  } else {1}}}
}  


# 1.5. Function to determine if a log reached the lower system boundary LSB
f.outflow<-function(x,m.NN_ID,m.t0,m.t1,i,j){
  if(j+i>2 & m.t1[x,"TS_in"]<=j){
    if(any(m.NN_ID[x,]%in%v.outnodes==T)){
      if(m.t0[x,"TS_out"]>0){
        m.t0[x,"TS_out"]
      } else{j}
    } else {0}
  } else{0}
}

# 1.6 Function to determine the timepoint of mobilisation
f.mobitime<-function(x,m.t0,m.t1,j){
  if(m.t0[x,"Status"]<3 & m.t1[x,"Status"]==3 & m.t0[x,"TS_mobi"]==0){
    j
  } else {m.t0[x,"TS_mobi"]}
}

# 1.7 Function to write the ID of the obstacle if a log gets jammed
f.vkl_obstacle<-function(x,m.t1,m.NN_ID,i,j){
  if(m.t1[x,"Obs_Nr"]==0){
    if(m.t1[x,"Status"]==4){
      nodenr<-m.NN_ID[x,which(m.NN_ID[x,]%in%m.Obstacles_Nodes[,"Node"]==T)[1]]
      m.Obstacles_Nodes[which(m.Obstacles_Nodes[,"Node"]==nodenr),"ID"]
    } else{0}
  }else{m.t1[x,"Obs_Nr"]}
}

#1.8 Function to write the timepoint of the entrapment
f.vkltime<-function(x,m.t0,m.t1,j){
  if(m.t1[x,"Status"]==4 & m.t0[x,"Status"]!=4 & m.t0[x,"TS_jam"]==0){
    j
  } else{m.t0[x,"TS_jam"]}
}

#----------------------------------------------------------

# 2. Functions for interpolation

# 2.1. Function to calculate the distance to the 3 nearest neighbors

f.dist<-function(x,m.NN_ID,m.t1,i,j){
  d1<-sqrt((d.Nods[m.NN_ID[x,1],"xcoord"]-m.t1[x,1])^2+(d.Nods[m.NN_ID[x,1],"ycoord"]-m.t1[x,2])^2)
  d2<-sqrt((d.Nods[m.NN_ID[x,2],"xcoord"]-m.t1[x,1])^2+(d.Nods[m.NN_ID[x,2],"ycoord"]-m.t1[x,2])^2)
  d3<-sqrt((d.Nods[m.NN_ID[x,3],"xcoord"]-m.t1[x,1])^2+(d.Nods[m.NN_ID[x,3],"ycoord"]-m.t1[x,2])^2)
  return(c(d1,d2,d3))
}
# 2.2. #IDW Interpolation

f.IDW<-function(x,attribute,v.dist,i,j){
  if(any(v.dist[,x]==0)){
    attribute[x,which(v.dist[,x]==0)]
  } else{attribute[x,1]*(v.dist[1,x]^(-r)/sum(v.dist[,x]^(-r)))+ attribute[x,2]*(v.dist[2,x]^(-r)/sum(v.dist[,x]^(-r)))+ attribute[x,3]*(v.dist[3,x]^(-r)/sum(v.dist[,x]^(-r)))}
}
  
    
 
#*******************************************************************************
# Functions for data input
#*******************************************************************************

# Function to determine the lower system boundary LSB
# The Nodes must be saved as STRINGDEF from Basement under "stringdef.txt"
f.read_outnodes<-function(x){
  v.outnodes<-   read.table(paste(dd,"/Inputdata/stringdef.txt",sep=""),skip=2,sep=" ",nrows=1)
  v.outnodes<-   v.outnodes[c(3:length(v.outnodes))]
  posf<-         as.numeric(strsplit(as.character(v.outnodes[1,1]),fixed=T,split="(")[[1]][2])
  posl<-         as.numeric(strsplit(as.character(v.outnodes[1,length(v.outnodes)]),fixed=T,split=")")[[1]])
  v.outnodes<-   c(posf,as.numeric(v.outnodes[1,c(3:length(v.outnodes)-1)]),posl)
  return(list(outnodes=v.outnodes))
}


#Function to load the obstacles
#The ID of the shp must equal their order of creation (top entry must be ID=1, second ID=2,...)

f.obstacle<-function(x){
  s.obstacles<-        readOGR(paste(dd,"/Inputdata",sep=""),layer="Obstacle")
  kbs<-CRS(proj4string(s.obstacles))
  
  m.obstacles_Nodes<-c()
  for(i in 1:dim(s.obstacles)[1]){
    m.coord<-            s.obstacles@polygons[[i]]@Polygons[[1]]@coords
    v.obstacle_nods<-    d.Nods[which(point.in.polygon(d.Nods[,3],d.Nods[,4],m.coord[,1],m.coord[,2])>0),"ID"]
    v.obstacle_ID<-      rep(i,length(v.obstacle_nods))
    m.obstacle<-         cbind(v.obstacle_ID,v.obstacle_nods)
    m.obstacles_Nodes<-  rbind(m.obstacles_Nodes,m.obstacle)
  }
  colnames(m.obstacles_Nodes)<-c("ID","Node")
  
  
  m.obstacles_Info<-   matrix(nrow=length(unique(m.obstacles_Nodes[,"ID"])),ncol=6,data=NA)
  colnames(m.obstacles_Info)<-c("ID","Height","Width","Parameter","Amount","Volume")
  for(i in 1:dim(s.obstacles)[1]){
    m.obstacles_Info[i,]<-as.numeric(s.obstacles@data[i,])
  }
  return(list(Nodes=m.obstacles_Nodes,Info=m.obstacles_Info,kbs=kbs))
  
}



#*******************************************************************************
# Functions for Postprocessing
#*******************************************************************************
#Function for time variation curve of deposition

f.voldepos<-function(x){
  treeID0<-which(a.res[,"Status",x-1]==2 & 
                   a.res[,"TS_mobi",x-1]>0  & 
                   a.res[,"xcoord",x-1]==a.res[,"xcoord",dim(a.res)[3]-1] & 
                   a.res[,"TS_out",x-1]==0)
  
  treeID<-which(a.res[,"Status",x]==2 & 
                  a.res[,"TS_mobi",x]>0  &
                  a.res[,"xcoord",x]==a.res[,"xcoord",dim(a.res)[3]-1] & 
                  a.res[,"TS_out",x]==0)
  sum((d.Wood[treeID,"DBH"]^2)*(d.Wood[treeID,"Length"]*0.4))-sum((d.Wood[treeID0,"DBH"]^2)*(d.Wood[treeID0,"Length"]*0.4))
}






