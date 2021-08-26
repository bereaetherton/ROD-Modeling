library(readr)
library(readxl)
library(rgdal)
library(terra)
library(geosphere)

ana.mat<- read_excel("Anohola.xlsx")
#ana.mat<- read_excel("/blue/garrett/betherton/ROD/Anohola.xlsx")
ana.mat<-as.data.frame(ana.mat)
dim(ana.mat)
ana.mat<-ana.mat[,c(1,3:20)]
ana.mat<-as.matrix(ana.mat)
#View(ana.mat)
ana.mat[which(ana.mat==9)]<-4
ana.mat[which(ana.mat!=4)]<-1
ana.mat[which(ana.mat==4)]<-0
ana.mat[,1]<-1:788
#View(ana)

sum(ana.mat[,2])
sum(ana.mat[,18])

colnames(ana.mat)<-c("Node","04/24/2020",	"05/06/2020",	"05/20/2020", "06/05/2020",
                     "06/26/2020", "07/10/2020", "07/23/2020", "08/07/2020", "08/19/2020",
                     "09/09/2020", "10/21/2020","11/05/2020", "11/20/2020", "12/30/2020",
                     "02/05/2021", "02/23/2021", "03/23/2021", "05/12/2021")
ana.utm<-read_excel("Anohola.xlsx")
#ana.utm<-read_excel("/blue/garrett/betherton/ROD/Anohola.xlsx")
ana.utm<-data.matrix(ana.utm)
ana.utm<-ana.utm[,c(21:22)]
ana.ll<-vect(ana.utm,crs="+proj=utm +zone=5n +datum=WGS84  +units=m")
ana<- project(ana.ll, "+proj=longlat +datum=WGS84")
ana<-geom(ana)[,c("x","y")]
colnames(ana)<-c("Lon","Lat")
rm(ana.utm)
rm(ana.ll)

#-------------------------------------------------------------------------------#

args<-commandArgs(TRUE)
ARG<-as.numeric(strsplit(args[1],",")[[1]])

ARG=c(5.0, 1.0)

print(ARG)

beta<-ARG[1]
alpha<-ARG[2]

dates<-list("Node","04/24/2020",	"05/06/2020",	"05/20/2020", "06/05/2020",
            "06/26/2020", "07/10/2020", "07/23/2020", "08/07/2020", "08/19/2020",
            "09/09/2020", "10/21/2020","11/05/2020", "11/20/2020", "12/30/2020",
            "02/05/2021", "02/23/2021", "03/23/2021", "05/12/2021")

days<-383
size<-dim(ana.mat)[1]

#-------------------------------------------------------------------------------#

ana.dist.mat<-matrix(0,nrow=size,ncol=size)
for(i in 1:size){
  for(j in 1:size){
    if(ana.dist.mat[j,i]==0){
      ana.dist.mat[i,j]<-distGeo(ana[(i),],ana[(j),])
    }
  }
}

ivp.ana.mat<-matrix(0,nrow=size,ncol=size)

#using distance matrix to creating the inverse power law probabilities

for(i in 1:size){
  for(j in 1:size){
    if(ana.dist.mat[i,j]==0){
      ivp.ana.mat[i,j]<-0
    }
    else{
      ivp.ana.mat[i,j]<-ana.dist.mat[i,j]^(-beta)
    }
  }
}


#using distance matrix to creating the negative exponential

#for(i in 1:size){
#  for(j in 1:size){
#    if(ana.dist.mat[i,j]==0){
#      ivp.ana.mat[i,j]<-0
#    }
#    else{
#      ivp.ana.mat[i,j]<-exp(ana.dist.mat[i,j]*(-beta))
#    }
#  }
#}

min<-min(ivp.ana.mat[which(ivp.ana.mat>0)])
max<-max(ivp.ana.mat)
mm<-max-min

#rescaling between 0 and 1

rs.ivp.ana<-matrix(0,nrow=size,ncol=size)

for(i in 1:size){
  for(j in 1:size){
    if(ivp.ana.mat[i,j]==0){
      rs.ivp.ana[i,j]<-0
    }
    else{
      rs.ivp.ana[i,j]<-((ivp.ana.mat[i,j]-min)/mm)
    }
  }
}

rm(ivp.ana.mat)
rm(max)
rm(min)
rm(mm)

rs.ivp.ana<-alpha*rs.ivp.ana

#-------------------------------------------------------------------------------#

init.dat.ana<-(ana.mat[,2])

ana.spread.mat<-matrix(0,nrow=size,ncol=days)
ana.spread.mat[,1]<-init.dat.ana

bp.spread<-function(k,w){
  s<-sample(c(0,1),1,replace = TRUE,prob=c(1-rs.ivp.ana[k,w],rs.ivp.ana[k,w]))
  return(s)
}

for(t in 2:days){
  for(i in 1:size){
    if(ana.spread.mat[i,t-1]==1){
      ana.spread.mat[i,t]<-1
    }
    else{
      for(j in 1:size){  
        if(rs.ivp.ana[i,j]==0){
          if(bp.spread(j,i)==1 && ana.spread.mat[j,t-1]==1){
            ana.spread.mat[i,t]<-1
            break
          }
        }
        else{
          if(bp.spread(i,j)==1 && ana.spread.mat[j,t-1]==1){
            ana.spread.mat[i,t]<-1
            break
          }
        }
      }
    }
  }
  print(t)
}

size<-dim(ana.mat)[1]
time<-dim(ana.mat)[2]-1

days.list<-matrix(0,nrow=1,ncol=time)
days.list[1]<-0
for(i in 2:time){
  days.list[i]<-as.numeric(difftime(as.Date(dates[[i+1]],"%m/%d/%Y"),as.Date(dates[[i]],"%m/%d/%Y"),units ="days"))
  days.list[i]<-days.list[i-1]+days.list[i]
}#will say out of bounds - NOT TO WORRY
days.list[1]<-1

#-------------------------------------------------------------------------------#

#NAME <- paste(ARG[1], ARG[2], sep="_")  
#save.image(paste0("/blue/garrett/betherton/ROD/","ANAIVP.test", NAME, ".RData"))
