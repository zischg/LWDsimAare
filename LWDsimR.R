LWDsimR<-function(){
  
  #check input parameter-----------------------------------------------
  if(!(TS==floor(TS))){
    stop("TS must me an integer.Execute Inputdaten.R again after adapting the Parameter")
  }
  if(TSW>TSB){
    stop("TSW must be smaller than TSB (TSW<TSB).Execute Inputdaten.R again after adapting the Parameter")
  }
  if(!(TSWT==floor(TSWT))){
    stop("TSWT must be an integer.Execute Inputdaten.R again after adapting the Parameter")
  }
  if(nsave<TSW){
    stop("nsave must be bigger than TSW.Execute Inputdaten.R again after adapting the Parameter")
  }  
  if(!(nsave/TSW==floor(nsave/TSW))){
    stop("nsave must be a multiple of TSW.Execute Inputdaten.R again after adapting the Parameter")
  } 
  if(dnr<1){
    stop("dnr must be bigger or equal 1.Execute Inputdaten.R again after adapting the Parameter")
  }
  if(dwr<1){
    stop("dwr must be bigger or equal 1.Execute Inputdaten.R again after adapting the Parameter")
  }
  
  #start of the two neested loops------------------------------------
    
  for(j in 1:TS){             # Outer loop: Time step of BASEMENT output (TSB)
    for(i in 1:TSWT){         # Inner loop: Time step of woody debris dynamics (TSW)
      start<-          Sys.time()
      
      #1. Localization-----------------------------------------------
      
      #Get information about hydrodynamics from the next 3 neighbours of each tree at every time step
      m.NN_ID<-        get.knnx(d.Nods[,c("xcoord","ycoord")],m.t1[,c("xcoord","ycoord")],k=3,algorithm = "kd_tree")$nn.index
      m.depth_TS<-     cbind(v.depth[m.NN_ID[,1]+(j-1)*n.nds+j],v.depth[m.NN_ID[,2]+(j-1)*n.nds+j],v.depth[m.NN_ID[,3]+(j-1)*n.nds+j])
      m.velx_TS<-      cbind(v.velx[m.NN_ID[,1]+(j-1)*n.nds+j],v.velx[m.NN_ID[,2]+(j-1)*n.nds+j],v.velx[m.NN_ID[,3]+(j-1)*n.nds+j])
      m.vely_TS<-      cbind(v.vely[m.NN_ID[,1]+(j-1)*n.nds+j],v.vely[m.NN_ID[,2]+(j-1)*n.nds+j],v.vely[m.NN_ID[,3]+(j-1)*n.nds+j])
      
      #Calculate distances to the nodes
      
      v.dist<-         sapply(1:length(d.Wood[,1]),f.dist,m.NN_ID=m.NN_ID,m.t1=m.t1,i=i,j=j)
      
      #IDW-Interpolation at location of every tree
      v.depth_TS<-     sapply(1:length(d.Wood[,1]),f.IDW,attribute=m.depth_TS,v.dist=v.dist,i=i,j=j)
      v.velx_TS<-      sapply(1:length(d.Wood[,1]),f.IDW,attribute=m.velx_TS,v.dist=v.dist,i=i,j=j)
      v.vely_TS<-      sapply(1:length(d.Wood[,1]),f.IDW,attribute=m.vely_TS,v.dist=v.dist,i=i,j=j)    
      
      #Check whether a log passed the LSB
      m.t1[,"TS_out"]<- sapply(c(1:length(d.Wood[,1])),f.outflow,m.NN_ID=m.NN_ID,m.t0=m.t0,m.t1=m.t1,i=i,j=j)      
      
      #2. Mobilization-----------------------------------------------
      
      m.t1[,"Status"]<- ifelse(m.t1[,"Status"]==1,
                               sapply(c(1:length(d.Wood[,1])),f.mobilisation,v.velx_TS=v.velx_TS,v.vely_TS=v.vely_TS,v.depth_TS=v.depth_TS,v.timeinwater=v.timeinwater),
                               m.t1[,"Status"])
      
      #3. Entrapment-----------------------------------------------
      
      m.t1[,"Status"]<- ifelse(m.t1[,"Status"]==3 & vkl==T,
                               sapply(c(1:length(d.Wood[,1])),f.verklausung,m.NN_ID=m.NN_ID,v.depth_TS=v.depth_TS,m.t1=m.t1,i=i,j=j),
                               m.t1[,"Status"])
      
      #4.Transport conditions---------------------------------------
      
      #Adapt x coordinates
      m.t2[,"xcoord"]<- ifelse(m.t1[,"TS_out"]==0 & m.t1[,"TS_in"]<=j,
                               ifelse(m.t1[,"Status"]==1,
                                      m.t1[,"xcoord"],
                                      ifelse(m.t1[,"Status"]==2 | m.t1[,"Status"]==3,
                                             sapply(c(1:length(d.Wood[,1])),f.trans_x,v.depth_TS=v.depth_TS,v.velx_TS=v.velx_TS,v.vely_TS=v.vely_TS,m.t1=m.t1,i=i,j=j),
                                             m.t1[,"xcoord"])
                               ), m.t1[,"xcoord"])
      
      
      #Adapt y coordinates
      m.t2[,"ycoord"]<- ifelse(m.t1[,"TS_out"]==0 & m.t1[,"TS_in"]<=j,
                               ifelse(m.t1[,"Status"]==1,
                                      m.t1[,"ycoord"],
                                      ifelse(m.t1[,"Status"]==2 | m.t1[,"Status"]==3,
                                             sapply(c(1:length(d.Wood[,1])),f.trans_y,v.depth_TS=v.depth_TS,v.velx_TS=v.velx_TS,v.vely_TS=v.vely_TS,m.t1=m.t1,i=i,j=j),
                                             m.t1[,"ycoord"])
                               ), m.t1[,"ycoord"])
      
      
      
      
      
      #5.Adapt Status and time points---------------------------------------
      
      #Define Status for next time step
      m.t2[,"Status"]<- sapply(c(1:length(d.Wood[,1])),f.status,i=i,j=j,m.t1=m.t1,m.t2=m.t2)
      
      #Write number of obstacle when bridge clogging
      m.t2[,"Obs_Nr"]<- sapply(c(1:length(d.Wood[,1])),f.vkl_obstacle,m.t1=m.t1,m.NN_ID=m.NN_ID,i=i,j=j)
      
      #Write time point of mobilization 
      m.t1[,"TS_mobi"]<-sapply(c(1:length(d.Wood[,1])),f.mobitime,m.t0=m.t0,m.t1=m.t1,j=j)
      
      #Write time point of entrapment
      m.t1[,"TS_jam"]<-sapply(c(1:length(d.Wood[,1])),f.vkltime,m.t0=m.t0,m.t1=m.t1,j=j)
      
      #6.Write Matrix into Array---------------------------------------
      
      if(any(i+(j-1)*TSWT==v.save)){
        a.out[,,which(v.save==i+(j-1)*TSWT)]<-m.t1
      }
      
      #7.Write safety copy on hard drive
      
      if(any(j==v.write & i==1)){
        save(list=ls(all.names = TRUE),file=paste(dd,"/Results/wscopy.RData",sep=""))
      }
      
      #8.Overwrite Matrices---------------------------------------
      
      m.t0<-m.t1
      m.t1<-m.t2
      m.t2[,c(1:5)]<-   m.tdefault[,c(1:5)]
      
      
      print(Sys.time()-start)
      print(c(j,i))
      print(paste(c(round(100/(TS*TSWT)*(i+(j-1)*TSWT),digits=2),"%"),collapse = ""))
    }
  }
  a.out<-abind(a.out,m.t1,along=3)
  
  
  return(a.out)
}

