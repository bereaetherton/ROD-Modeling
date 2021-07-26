#-------------------------------------------------------------#

library(readxl)
library(rgdal)
library(terra)
library(geosphere)

#-------------------------------------------------------------#

#kapap data read in 

kapap<-read_excel("Kapapala_2021_02_05_for_Karen (1).xlsx",col_types = c("numeric"),
                  sheet = "FINAL")
#kapap<-read_excel("/blue/garrett/betherton/ROD/Kapapala_2021_02_05_for_Karen (1).xlsx",col_types = c("numeric"),
#                  sheet = "FINAL")
kapap<-as.matrix(kapap)
dim(kapap)
kapap<-kapap[c(1:4134,4136:4205),c(1:54)]
#View(kapap)
colnames(kapap)<-c("Node","06/10/2010",	"01/28/2011",	"01/14/2012",	"10/30/2014",	"09/25/2015", "10/26/2015",	"01/09/2016","02/19/2016","02/20/2016","11/29/2016","03/16/2017",	
                   "04/20/2017",	"05/30/2017",	"06/17/2017",	"08/08/2017",	"10/09/2017", "11/02/2017",	"12/11/2017",	"01/10/2018",	"02/09/2018",	"03/16/2018",	
                   "04/11/2018",	"05/14/2018",	"06/12/2018", "07/19/2018",	"08/20/2018",	"10/03/2018",	"11/08/2018",	"12/13/2018",	"01/11/2019",	"01/12/2019",	
                   "02/26/2019", "03/22/2019",	"04/30/2019",	"05/23/2019",	"07/02/2019","07/25/2019","08/22/2019",	"09/30/2019",	"10/28/2019", "11/27/2019",	"01/09/2020",
                   "02/12/2020",	"03/27/2020",	"04/30/2020",	"06/10/2020","07/10/2020","08/21/2020",	"09/10/2020",	"10/19/2020",	"11/19/2020",	"12/17/2020",	"01/20/2021")
kapap[,1]<-1:4204

#kapap UTM to lat and lon points


kp.utm<-read_excel("Kapapala_2021_02_05_for_Karen (1).xlsx",col_types = c("numeric"),
                  sheet = "FINAL")
#kp.utm<-read_excel("/blue/garrett/betherton/ROD/Kapapala_2021_02_05_for_Karen (1).xlsx",col_types = c("numeric"),
#                          sheet = "FINAL")
kp.utm<-data.matrix(kp.utm)
kp.utm<-kp.utm[c(1:4134,4136:4205),c(55,56)]
kp.ll<-vect(kp.utm,crs="+proj=utm +zone=5n +datum=WGS84  +units=m")
kp<- project(kp.ll, "+proj=longlat +datum=WGS84")
kp<-geom(kp)[,c("x","y")]
colnames(kp)<-c("Lon","Lat")
rm(kp.utm)
rm(kp.ll)

#-----------------------------------------------------------------------------#

#read in the arguments for beta and alpha

args<-commandArgs(TRUE)
ARG<-as.numeric(strsplit(args[1],",")[[1]])

ARG=c(5.0, 0.9)

print(ARG)

beta<-ARG[1]
alpha<-ARG[2]

size=dim(kp)[1]
time=dim(kapap)[2]-1

#-------------------------------------------------------------#

#creating a distance matrix 

kapap<-kapap[1:300,]
kp<-kp[1:300,]
size<-dim(kapap)[1]

kp.dist.mat<-matrix(0,nrow=size,ncol=size)
for(i in 1:size){
 for(j in 1:size){
    if(kp.dist.mat[j,i]==0){
      kp.dist.mat[i,j]<-distGeo(kp[(i),],kp[(j),])
    }
  }
}

ivp.kp.mat<-matrix(0,nrow=size,ncol=size)

#using distance matrix to creating the inverse power law probabilities

for(i in 1:size){
  for(j in 1:size){
    if(kp.dist.mat[i,j]==0){
      ivp.kp.mat[i,j]<-0
    }
    else{
      ivp.kp.mat[i,j]<-kp.dist.mat[i,j]^(-beta)
    }
  }
}

min<-min(ivp.kp.mat[which(ivp.kp.mat>0)])
max<-max(ivp.kp.mat)
mm<-max-min

#rescaling between 0 and 1

rs.ivp.kp<-matrix(0,nrow=size,ncol=size)

for(i in 1:size){
  for(j in 1:size){
    if(ivp.kp.mat[i,j]==0){
      rs.ivp.kp[i,j]<-0
    }
    else{
      rs.ivp.kp[i,j]<-((ivp.kp.mat[i,j]-min)/mm)
    }
  }
}

rm(ivp.kp.mat)
rm(max)
rm(min)
rm(mm)

rs.ivp.kp<-alpha*rs.ivp.kp

#-----------------------------------------------------------------#

#remove the -3 values and change the node values (cleaning up)

kapap[which(kapap<4)]<-1 
kapap[which(kapap==4)]<-0 
kapap[,1]<-1:4204
#View(kapap)

dates<-list("Node","06/10/2010",	"01/28/2011",	"01/14/2012",	"10/30/2014",	"09/25/2015", "10/26/2015",	"01/09/2016","02/19/2016","02/20/2016","11/29/2016","03/16/2017",	
         "04/20/2017",	"05/30/2017",	"06/17/2017",	"08/08/2017",	"10/09/2017", "11/02/2017",	"12/11/2017",	"01/10/2018",	"02/09/2018",	"03/16/2018",	
         "04/11/2018",	"05/14/2018",	"06/12/2018", "07/19/2018",	"08/20/2018",	"10/03/2018",	"11/08/2018",	"12/13/2018",	"01/11/2019",	"01/12/2019",	
         "02/26/2019", "03/22/2019",	"04/30/2019",	"05/23/2019",	"07/02/2019","07/25/2019","08/22/2019",	"09/30/2019",	"10/28/2019", "11/27/2019",	"01/09/2020",
         "02/12/2020",	"03/27/2020",	"04/30/2020",	"06/10/2020","07/10/2020","08/21/2020",	"09/10/2020",	"10/19/2020",	"11/19/2020",	"12/17/2020",	"01/20/2021")

#simulating disease spread

init.dat.kp<-(kapap[,10])

days<-1796 #988 days between 5/4/2017 and 1/17/2020

#1796 days between kapap[,10] and kapap[,53]

kp.spread.mat<-matrix(0,nrow=size,ncol=days)
kp.spread.mat[,1]<-init.dat.kp

bp.spread<-function(k,w){
  s<-sample(c(0,1),1,replace = TRUE,prob=c(1-rs.ivp.kp[k,w],rs.ivp.kp[k,w]))
  return(s)
}

for(t in 2:days){
  for(i in 1:size){
    if(kp.spread.mat[i,t-1]==1){
      kp.spread.mat[i,t]<-1
    }
    else{
      for(j in 1:size){  
        if(rs.ivp.kp[i,j]==0){
            if(bp.spread(j,i)==1 && kp.spread.mat[j,t-1]==1){
              kp.spread.mat[i,t]<-1
              break
            }
          }
        else{
            if(bp.spread(i,j)==1 && kp.spread.mat[j,t-1]==1){
              kp.spread.mat[i,t]<-1
              break
            }
         }
       }
     }
   }
  print(t)
}
#View(kp.spread.mat)

#------------------------------------------------------------------------------#
days.list<-matrix(0,nrow=1,ncol=time)
days.list[1]<-0
for(i in 2:time){
  days.list[i]<-as.numeric(difftime(as.Date(dates[[i+1]],"%m/%d/%Y"),as.Date(dates[[i]],"%m/%d/%Y"),units ="days"))
  days.list[i]<-days.list[i-1]+days.list[i]
}

dl<-days.list[9:53]-2081
dl[1]<-1

#filtering out specific days

kpivp<-kp.spread.mat[,dl]
beep(sound=2)

#------------------------------------------------------------------------------#


#saving the image to the working directory

NAME <- paste(ARG[1], ARG[2], sep="_")  
save.image("KPROD.test.RData")

